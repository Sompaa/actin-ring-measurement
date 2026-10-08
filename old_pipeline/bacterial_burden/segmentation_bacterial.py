"""
Segmentation Module for Bacterial Burden Quantification
Handles 3-channel 4D images [Z, C=3, Y, X] with Best Focus selection.
- Host cells: Cellpose cyto3
- Nuclei: Cellpose nuclei
- Bacteria: Sato filter (per-cell normalized) for area measurement
"""

import numpy as np
import tifffile
from pathlib import Path
from scipy.ndimage import laplace
from cellpose import models, core
import gc


class BacterialSegmenter:
    """
    Three-channel segmentation for host-pathogen interaction analysis.

    Handles 4D [Z, C=3, Y, X] images with:
    - Best Focus selection (Variance of Laplacian on host channel)
    - Host cell segmentation (Cellpose cyto3)
    - Nuclei segmentation (Cellpose nuclei)
    - Bacteria detection (Sato filter with per-cell normalization)
    """

    def __init__(self, config):
        self.config = config
        gpu_available = core.use_gpu()
        self.use_gpu = config['processing']['use_gpu'] and gpu_available

        # Channel mapping
        self.ch_host = config['channels']['host_cytoplasm']
        self.ch_nuclei = config['channels']['nuclei']
        self.ch_bacteria = config['channels']['bacteria']

        # QC thresholds
        self.qc_config = config.get('qc', {})

        # Z-projection method
        self.z_method = config.get('z_projection', {}).get('method', 'best_focus')

        if config['processing']['use_gpu'] and not gpu_available:
            print('GPU requested but not available, using CPU')

        print(f'BacterialSegmenter ready ({"GPU" if self.use_gpu else "CPU"})')

    def _adjust_brightness_contrast(self, image, params):
        """
        Adjust brightness/contrast with gamma correction.
        Algorithm from segmentation_v3.py CellSegmenter.adjust_brightness_contrast.
        """
        img = image.astype(np.float32)
        img = (img - img.min()) / (img.max() - img.min() + 1e-10)
        img = img * params['brightness']
        mean = img.mean()
        img = (img - mean) * params['contrast'] + mean
        img = np.power(np.clip(img, 1e-10, None), 1.0 / params['gamma'])
        img = np.clip(img, 0, 1)
        img = (img * 65535).astype(np.uint16)
        return img

    def _normalize_cell_region(self, region, method='percentile', percentile_low=2, percentile_high=98):
        """Delegate to shared detection module."""
        from bacterial_burden.detection import normalize_cell_region
        return normalize_cell_region(region, method, percentile_low, percentile_high)

    def _detect_bacteria_sato(self, bacteria_region, params):
        """Delegate to shared detection module."""
        from bacterial_burden.detection import detect_bacteria_in_cell
        return detect_bacteria_in_cell(bacteria_region, params)

    def _quantify_bacteria_per_cell(self, bacteria_channel, host_masks, valid_cell_ids, config_params, bacteria_raw=None):
        """
        Quantify bacterial area per cell using Sato filter + thresholding.

        NO individual bacteria segmentation - just measures total bacterial area.

        Args:
            bacteria_channel: Full 2D bacteria channel image (B/C adjusted)
            bacteria_raw: Raw bacteria channel (before B/C) for signal range check
            host_masks: 2D label array of host cells
            valid_cell_ids: List of cell IDs to process
            config_params: Dict with bacteria detection parameters

        Returns:
            Dict mapping cell_id -> {
                'bacteria_mask': binary mask of bacteria in this cell,
                'bacterial_area': total bacterial pixels,
                'host_area': host cell area in pixels,
                'area_ratio': bacterial_area / host_area
            }
        """
        from skimage.measure import regionprops

        results = {}

        # Get normalization params
        norm_params = self.config.get('per_cell_normalization', {
            'enabled': True,
            'method': 'percentile',
            'percentile_low': 2,
            'percentile_high': 98
        })

        # Sato filter params are directly in config_params (bacteria_detection)
        sato_params = config_params
        padding = config_params.get('padding', 15)

        min_range = config_params.get('min_signal_range', 100)
        print(f"Quantifying bacterial area for {len(valid_cell_ids)} cells (min_signal_range={min_range})")

        skipped_low_signal = 0
        signal_ranges = []

        for cell_id in valid_cell_ids:
            # Extract cell region
            cell_mask = (host_masks == cell_id)
            props = regionprops(cell_mask.astype(int))

            if len(props) == 0:
                continue

            bbox = props[0].bbox
            host_area = props[0].area

            # Add padding
            r0 = max(0, bbox[0] - padding)
            c0 = max(0, bbox[1] - padding)
            r1 = min(bacteria_channel.shape[0], bbox[2] + padding)
            c1 = min(bacteria_channel.shape[1], bbox[3] + padding)

            cell_mask_roi = cell_mask[r0:r1, c0:c1]

            # Check if bacteria signal inside the cell is above noise
            # Uses RAW image (before brightness/contrast) so noise isn't inflated
            # On control images (no bacteria), the raw signal is flat background
            raw_check = bacteria_raw if bacteria_raw is not None else bacteria_channel
            raw_roi = raw_check[r0:r1, c0:c1]
            cell_bacteria_raw = raw_roi[cell_mask_roi]
            if cell_bacteria_raw.size > 0 and min_range > 0:
                p2 = float(np.percentile(cell_bacteria_raw, 2))
                p98 = float(np.percentile(cell_bacteria_raw, 98))
                signal_range = p98 - p2
                signal_ranges.append(signal_range)
                if signal_range < min_range:
                    skipped_low_signal += 1
                    results[cell_id] = {
                        'bacteria_mask': np.zeros_like(cell_mask, dtype=bool),
                        'bacterial_area': 0,
                        'host_area': host_area,
                        'area_ratio': 0.0
                    }
                    continue

            # Extract ROI (adjusted image for actual detection)
            bacteria_roi = bacteria_channel[r0:r1, c0:c1].copy()

            if bacteria_roi.size == 0:
                results[cell_id] = {
                    'bacteria_mask': np.zeros_like(cell_mask, dtype=bool),
                    'bacterial_area': 0,
                    'host_area': host_area,
                    'area_ratio': 0.0
                }
                continue

            # PER-CELL NORMALIZATION
            if norm_params.get('enabled', True):
                bacteria_roi = self._normalize_cell_region(
                    bacteria_roi,
                    method=norm_params.get('method', 'percentile'),
                    percentile_low=norm_params.get('percentile_low', 2),
                    percentile_high=norm_params.get('percentile_high', 98)
                )

            # SATO FILTER + THRESHOLDING
            bacteria_mask_roi = self._detect_bacteria_sato(bacteria_roi, sato_params)

            # Restrict to cell interior only
            bacteria_mask_roi = bacteria_mask_roi & cell_mask_roi

            # Create full-size mask
            bacteria_mask_full = np.zeros_like(cell_mask, dtype=bool)
            bacteria_mask_full[r0:r1, c0:c1] = bacteria_mask_roi

            # Calculate area
            bacterial_area = np.sum(bacteria_mask_full)
            area_ratio = bacterial_area / float(host_area) if host_area > 0 else 0.0

            results[cell_id] = {
                'bacteria_mask': bacteria_mask_full,
                'bacterial_area': bacterial_area,
                'host_area': host_area,
                'area_ratio': area_ratio
            }

        if signal_ranges:
            sr = np.array(signal_ranges)
            print(f"  Raw signal ranges: min={sr.min():.0f}, median={np.median(sr):.0f}, max={sr.max():.0f} (threshold={min_range})")
        if skipped_low_signal > 0:
            print(f"  Skipped {skipped_low_signal}/{len(valid_cell_ids)} cells (low signal)")
        detected = sum(1 for r in results.values() if r['bacterial_area'] > 0)
        print(f"Quantified {len(results)} cells ({detected} with bacteria, {skipped_low_signal} skipped low signal)")

        return results

    def _compute_variance_of_laplacian(self, image_2d):
        """
        Compute Variance of Laplacian for a single 2D plane.
        Higher value = sharper focus.
        """
        lap = laplace(image_2d.astype(np.float64))
        return np.var(lap)

    def select_best_focus(self, z_stack):
        """
        Select best-focus Z-plane using Variance of Laplacian.

        Parameters
        ----------
        z_stack : ndarray, shape [Z, Y, X]
            Single-channel Z-stack (typically host channel)

        Returns
        -------
        int
            Index of best-focus Z-plane
        """
        n_z = z_stack.shape[0]
        scores = np.zeros(n_z)
        for z in range(n_z):
            scores[z] = self._compute_variance_of_laplacian(z_stack[z])
        z_best = int(np.argmax(scores))
        return z_best

    def _apply_z_projection(self, z_stack):
        """Apply MIP or middle-slice projection as fallback methods."""
        if self.z_method == 'mip':
            return np.max(z_stack, axis=0)
        elif self.z_method == 'single_middle':
            mid = z_stack.shape[0] // 2
            return z_stack[mid]
        else:
            return np.max(z_stack, axis=0)

    def qc_check_channel(self, image_2d, channel_name):
        """
        Quality control check for a single channel image.

        Only empty channels (no signal) cause a failure. Saturation is logged
        as a warning but does NOT fail the image — saturated host/nuclei
        channels are still usable for segmentation, and control images often
        have bright bacteria channels with no actual bacteria.

        Returns
        -------
        dict
            Keys: pass (bool), max_intensity, mean_intensity,
                  reason (str or None), saturated (bool)
        """
        max_val = float(image_2d.max())
        mean_val = float(image_2d.mean())

        min_thresh = self.qc_config.get('min_intensity_threshold', 50)
        if max_val < min_thresh:
            return {
                'pass': False, 'max_intensity': max_val,
                'mean_intensity': mean_val, 'reason': f'{channel_name}_empty',
                'saturated': False
            }

        # Check saturation — log but don't fail
        saturated = False
        if image_2d.dtype == np.uint16:
            sat_frac = float(np.mean(image_2d >= 65534))
        elif image_2d.dtype == np.uint8:
            sat_frac = float(np.mean(image_2d >= 254))
        else:
            sat_frac = 0.0

        sat_max = self.qc_config.get('saturation_fraction_max', 0.01)
        if sat_frac > sat_max:
            saturated = True

        return {
            'pass': True, 'max_intensity': max_val,
            'mean_intensity': mean_val, 'reason': None,
            'saturated': saturated
        }

    def load_and_preprocess(self, image_path):
        """
        Load 3-channel image, apply best-focus selection, run QC.

        Handles:
        - 4D [Z, C=3, Y, X]: best-focus or MIP, then extract 3 channels
        - 3D [C=3, Y, X]: direct channel split (no Z-projection)

        Returns
        -------
        dict or None
            Keys: host_2d, nuclei_2d, bacteria_2d, z_best, qc_results
        """
        try:
            img_stack = tifffile.imread(str(image_path))
        except Exception as e:
            print(f'Error reading {image_path}: {e}')
            return None

        z_best = -1

        if img_stack.ndim == 4:
            # [Z, C, Y, X]
            n_z, n_c, h, w = img_stack.shape

            if n_c < 3:
                print(f'Expected 3 channels, got {n_c} in {Path(image_path).name}')
                return None

            # Best focus selection on host channel
            if self.z_method == 'best_focus':
                host_z_stack = img_stack[:, self.ch_host, :, :]
                z_best = self.select_best_focus(host_z_stack)
                host_2d = img_stack[z_best, self.ch_host, :, :]
                nuclei_2d = img_stack[z_best, self.ch_nuclei, :, :]
                bacteria_2d = img_stack[z_best, self.ch_bacteria, :, :]
            else:
                # MIP or single_middle
                host_2d = self._apply_z_projection(img_stack[:, self.ch_host, :, :])
                nuclei_2d = self._apply_z_projection(img_stack[:, self.ch_nuclei, :, :])
                bacteria_2d = self._apply_z_projection(img_stack[:, self.ch_bacteria, :, :])

        elif img_stack.ndim == 3:
            # [C, Y, X] - no Z-stack
            n_c = img_stack.shape[0]
            if n_c < 3:
                print(f'Expected 3 channels, got {n_c} in {Path(image_path).name}')
                return None

            host_2d = img_stack[self.ch_host]
            nuclei_2d = img_stack[self.ch_nuclei]
            bacteria_2d = img_stack[self.ch_bacteria]

        else:
            print(f'Unexpected dimensions {img_stack.shape} for {Path(image_path).name}')
            return None

        # Free the full stack
        del img_stack
        gc.collect()

        # Keep raw bacteria for signal range check (before brightness/contrast)
        bacteria_2d_raw = bacteria_2d.copy()

        # Apply brightness/contrast per channel
        bc_config = self.config['brightness_contrast']
        host_2d = self._adjust_brightness_contrast(host_2d, bc_config['host'])
        nuclei_2d = self._adjust_brightness_contrast(nuclei_2d, bc_config['nuclei'])
        bacteria_2d = self._adjust_brightness_contrast(bacteria_2d, bc_config['bacteria'])

        # QC checks
        qc_results = {
            'host': self.qc_check_channel(host_2d, 'host'),
            'nuclei': self.qc_check_channel(nuclei_2d, 'nuclei'),
            'bacteria': self.qc_check_channel(bacteria_2d, 'bacteria'),
        }

        # Log QC issues
        for ch_name, qc in qc_results.items():
            if not qc['pass']:
                print(f'  QC SKIP [{Path(image_path).name}] {ch_name}: {qc["reason"]} '
                      f'(max={qc["max_intensity"]:.0f}, mean={qc["mean_intensity"]:.1f})')
            elif qc.get('saturated', False):
                print(f'  QC WARN [{Path(image_path).name}] {ch_name}: saturated '
                      f'(max={qc["max_intensity"]:.0f}, mean={qc["mean_intensity"]:.1f})')

        return {
            'host_2d': host_2d,
            'nuclei_2d': nuclei_2d,
            'bacteria_2d': bacteria_2d,
            'bacteria_2d_raw': bacteria_2d_raw,
            'z_best': z_best,
            'qc_results': qc_results
        }

    def segment_host_cells(self, host_2d, nuclei_2d=None):
        """
        Segment host cells using Cellpose cyto3 with optional nuclei seed.

        Parameters
        ----------
        host_2d : ndarray [Y, X]
        nuclei_2d : ndarray [Y, X] or None

        Returns
        -------
        ndarray
            Labeled host cell masks
        """
        cfg = self.config['segmentation']['host_cell']
        model = models.CellposeModel(gpu=self.use_gpu, pretrained_model=cfg['model'])

        if cfg.get('use_nuclei_channel', True) and nuclei_2d is not None:
            stacked = np.stack([host_2d, nuclei_2d], axis=-1)
        else:
            stacked = np.stack([host_2d, host_2d], axis=-1)

        masks, _, _ = model.eval(
            stacked,
            diameter=cfg['diameter'] if cfg['diameter'] > 0 else None,
            flow_threshold=cfg['flow_threshold'],
            cellprob_threshold=cfg['cellprob_threshold'],
            min_size=cfg.get('min_size', -1)
        )
        del model
        gc.collect()
        return masks

    def segment_nuclei(self, nuclei_2d):
        """
        Segment nuclei using Cellpose nuclei model.

        Parameters
        ----------
        nuclei_2d : ndarray [Y, X]

        Returns
        -------
        ndarray
            Labeled nucleus masks
        """
        cfg = self.config['segmentation']['nucleus']
        model = models.CellposeModel(gpu=self.use_gpu, pretrained_model=cfg['model'])

        masks, _, _ = model.eval(
            nuclei_2d,
            diameter=cfg['diameter'] if cfg['diameter'] > 0 else None,
            flow_threshold=cfg['flow_threshold'],
            cellprob_threshold=cfg['cellprob_threshold'],
            min_size=cfg.get('min_size', -1)
        )
        del model
        gc.collect()
        return masks

    def segment_image(self, image_path):
        """
        Full segmentation pipeline: load -> best focus -> QC -> segment host/nuclei
        -> detect bacteria per cell via Sato filter.

        Parameters
        ----------
        image_path : str or Path

        Returns
        -------
        dict or None
            Keys: host_masks, nucleus_masks, bacteria_masks,
                  host_2d, nuclei_2d, bacteria_2d, z_best, qc_results
        """
        preproc = self.load_and_preprocess(image_path)
        if preproc is None:
            return None

        # Early QC check: if any channel fails, return early with QC info
        # so the caller can skip without running expensive segmentation
        qc_results = preproc.get('qc_results', {})
        if self.config['qc'].get('skip_empty_images', True):
            if any(not ch_qc['pass'] for ch_qc in qc_results.values()):
                return {
                    'host_masks': np.zeros_like(preproc['host_2d'], dtype=np.int32),
                    'nucleus_masks': np.zeros_like(preproc['nuclei_2d'], dtype=np.int32),
                    'bacteria_masks': np.zeros_like(preproc['bacteria_2d'], dtype=np.uint32),
                    'host_2d': preproc['host_2d'],
                    'nuclei_2d': preproc['nuclei_2d'],
                    'bacteria_2d': preproc['bacteria_2d'],
                    'z_best': preproc.get('z_best', -1),
                    'qc_results': qc_results
                }

        host_2d = preproc['host_2d']
        nuclei_2d = preproc['nuclei_2d']
        bacteria_2d = preproc['bacteria_2d']
        bacteria_2d_raw = preproc.get('bacteria_2d_raw', bacteria_2d)

        # Segment host and nuclei
        host_masks = self.segment_host_cells(host_2d, nuclei_2d)
        nucleus_masks = self.segment_nuclei(nuclei_2d)

        # Filter host cells (need valid cell IDs for per-cell bacteria detection)
        from filtering_bacterial import HostCellFilter
        cell_filter = HostCellFilter(self.config)
        filter_result = cell_filter.filter_host_cells(host_masks.copy(), nucleus_masks)
        valid_cell_ids = [c['cell_id'] for c in filter_result['valid_cells']]

        # Per-cell bacteria detection via Sato filter
        bact_params = self.config.get('bacteria_detection', {})
        bacteria_results = self._quantify_bacteria_per_cell(
            bacteria_2d, host_masks, valid_cell_ids, bact_params,
            bacteria_raw=bacteria_2d_raw
        )

        # Create combined mask for visualization (bacteria labeled by parent cell ID)
        bacteria_masks = np.zeros_like(bacteria_2d, dtype=np.uint32)
        for cell_id, result in bacteria_results.items():
            bacteria_masks[result['bacteria_mask']] = cell_id

        # Store results for quantification module
        self._bacteria_area_results = bacteria_results

        # Save masks if requested
        if self.config['processing'].get('save_masks', False):
            output_dir = Path(self.config['output_dir']) / 'masks'
            output_dir.mkdir(parents=True, exist_ok=True)
            img_name = Path(image_path).stem
            tifffile.imwrite(str(output_dir / f'{img_name}_host_mask.tif'),
                             host_masks.astype(np.uint16))
            tifffile.imwrite(str(output_dir / f'{img_name}_nucleus_mask.tif'),
                             nucleus_masks.astype(np.uint16))
            tifffile.imwrite(str(output_dir / f'{img_name}_bacteria_mask.tif'),
                             bacteria_masks.astype(np.uint16))

        return {
            'host_masks': host_masks,
            'nucleus_masks': nucleus_masks,
            'bacteria_masks': bacteria_masks,
            'host_2d': host_2d,
            'nuclei_2d': nuclei_2d,
            'bacteria_2d': bacteria_2d,
            'z_best': preproc['z_best'],
            'qc_results': preproc['qc_results']
        }
