"""
Segmentation Module V3 - N-Channel Support with Optional Brightfield
Cellpose v4 with configurable channel indices and brightfield cell segmentation.
Supports 2D, 3D Z-stacks, and 4D multi-channel Z-stacks (2-4 channels).
"""

import numpy as np
import tifffile
from pathlib import Path
from cellpose import models, core
import gc


class CellSegmenter:
    """Configurable N-channel segmentation with optional brightfield and bacteria support."""

    def __init__(self, config):
        self.config = config
        gpu_available = core.use_gpu()
        self.use_gpu = config['processing']['use_gpu'] and gpu_available

        # Z-projection settings
        self.z_config = config.get('z_projection', {
            'method': 'mip',
            'use_for_segmentation': True,
            'keep_3d_volume': False
        })
        self.z_projection_method = self.z_config.get('method', 'mip')
        self.keep_3d = self.z_config.get('keep_3d_volume', False)

        # Channel indices (config-driven, backward compatible)
        bact_config = config.get('bacterial_burden', {})
        if bact_config.get('enabled', False):
            ch_map = bact_config.get('channels', {})
            self.ch_actin = ch_map.get('actin', 0)
            self.ch_nuclei = ch_map.get('nuclei', 1)
            self.ch_bacteria = ch_map.get('bacteria', None)
            self.ch_brightfield = ch_map.get('brightfield', None)
            self.cell_seg_channel = bact_config.get('cell_segmentation_channel', 'actin')
            self.bacterial_enabled = True
            self.bact_bc_params = bact_config.get('brightness_contrast', {})
        else:
            self.ch_actin = 0
            self.ch_nuclei = 1
            self.ch_bacteria = None
            self.ch_brightfield = None
            self.cell_seg_channel = 'actin'
            self.bacterial_enabled = False
            self.bact_bc_params = {}

        # Image-level normalization (CLAHE on full actin MIP)
        img_norm = config.get('image_normalization', {})
        self.img_clahe_enabled = img_norm.get('use_clahe', False)
        self.img_clahe_clip = img_norm.get('clahe_clip_limit', 0.03)
        self.img_clahe_tile = tuple(img_norm.get('clahe_tile_grid_size', [8, 8]))

        if config['processing']['use_gpu'] and not gpu_available:
            print('WARNING: GPU requested but not available, using CPU')

        # Load Cellpose models ONCE (not per-image)
        self.cell_cfg = config['segmentation']['cell']
        self.nucleus_cfg = config['segmentation']['nucleus']

        self.cell_model = models.CellposeModel(gpu=self.use_gpu, pretrained_model='cpsam')
        self.nucleus_model = models.CellposeModel(gpu=self.use_gpu, pretrained_model='cpsam')

        print(f'Cellpose will use: {"GPU" if self.use_gpu else "CPU"}')
        if self.bacterial_enabled:
            print(f'Bacterial burden: ENABLED (channels: actin={self.ch_actin}, nuclei={self.ch_nuclei}, '
                  f'bacteria={self.ch_bacteria}, brightfield={self.ch_brightfield})')
            print(f'Cell segmentation channel: {self.cell_seg_channel}')
        print('Cellpose models loaded (cell + nucleus)')

    def adjust_brightness_contrast(self, image, params):
        """Adjust brightness/contrast with gamma correction."""
        img = image.astype(np.float32)
        img = (img - img.min()) / (img.max() - img.min() + 1e-10)
        img = img * params['brightness']
        mean = img.mean()
        img = (img - mean) * params['contrast'] + mean
        img = np.power(np.clip(img, 1e-10, None), 1.0 / params['gamma'])
        img = np.clip(img, 0, 1)
        img = (img * 65535).astype(np.uint16)
        return img

    def apply_z_projection(self, z_stack, axis=0):
        """Apply Z-projection using configured method (MIP, mean, or single_middle)."""
        if self.z_projection_method == 'mip':
            return np.max(z_stack, axis=axis)
        elif self.z_projection_method == 'mean':
            return np.mean(z_stack, axis=axis).astype(z_stack.dtype)
        elif self.z_projection_method == 'single_middle':
            mid_slice = z_stack.shape[axis] // 2
            return np.take(z_stack, mid_slice, axis=axis)
        else:
            print(f'WARNING: Unknown projection method: {self.z_projection_method}, using MIP')
            return np.max(z_stack, axis=axis)

    def _validate_image(self, img_stack, image_path):
        """Validate loaded image array. Returns error message or None if valid."""
        if img_stack.size == 0:
            return f"Empty image file: {image_path}"
        if img_stack.ndim < 2 or img_stack.ndim > 4:
            return f"Unsupported dimensions ({img_stack.ndim}D): {image_path}"
        if not np.issubdtype(img_stack.dtype, np.number):
            return f"Non-numeric dtype ({img_stack.dtype}): {image_path}"
        if img_stack.max() == 0:
            return f"All-zero image: {image_path}"
        return None

    def _validate_channel_count(self, n_channels, image_path):
        """Validate that the image has enough channels for the config.
        Only actin and nuclei are mandatory; bacteria/brightfield are optional."""
        required = max(self.ch_actin, self.ch_nuclei) + 1
        if n_channels < required:
            print(f"WARNING: Image has {n_channels} channels but config requires {required} "
                  f"(actin+nuclei): {Path(image_path).name}")
            return False
        if self.ch_bacteria is not None and n_channels <= self.ch_bacteria:
            print(f"  INFO: bacteria channel ({self.ch_bacteria}) not available in "
                  f"{n_channels}-channel image — skipping bacteria for {Path(image_path).name}")
        if self.ch_brightfield is not None and n_channels <= self.ch_brightfield:
            print(f"  INFO: brightfield channel ({self.ch_brightfield}) not available in "
                  f"{n_channels}-channel image — skipping brightfield for {Path(image_path).name}")
        return True

    def load_and_preprocess(self, image_path):
        """
        Load and preprocess image with N-channel support.

        Returns dict with keys:
            actin_2d, nuclei_2d, actin_2d_raw: always present
            actin_3d, actin_3d_raw: present if keep_3d and Z-stack input
            bacteria_2d, bacteria_2d_raw: present if bacterial_burden enabled
            brightfield_2d: present if brightfield channel configured
            original_ndim: dimensionality of source image
        Returns None on error.
        """
        try:
            img_stack = tifffile.imread(image_path)
        except Exception as e:
            print(f'Error reading {image_path}: {e}')
            return None

        # Input validation
        error = self._validate_image(img_stack, image_path)
        if error:
            print(f"WARNING: {error}")
            return None

        original_ndim = img_stack.ndim

        # Initialize all outputs
        actin_2d = None
        nuclei_2d = None
        actin_3d = None
        bacteria_2d = None
        bacteria_2d_raw = None
        brightfield_2d = None

        # --- Handle different input dimensions ---
        if img_stack.ndim == 2:
            # Pure 2D single-channel image
            actin_2d = img_stack
            nuclei_2d = None

        elif img_stack.ndim == 3:
            n_first = img_stack.shape[0]
            if n_first <= 4:
                # Multi-channel 2D: [C, Y, X]
                actin_2d = img_stack[self.ch_actin]
                if n_first > self.ch_nuclei:
                    nuclei_2d = img_stack[self.ch_nuclei]
                # Extract bacteria if enabled
                if self.bacterial_enabled and self.ch_bacteria is not None and n_first > self.ch_bacteria:
                    bacteria_2d = img_stack[self.ch_bacteria].copy()
                    bacteria_2d_raw = bacteria_2d.copy()
                # Extract brightfield if configured
                if self.ch_brightfield is not None and n_first > self.ch_brightfield:
                    brightfield_2d = img_stack[self.ch_brightfield].copy()
            else:
                # Single-channel Z-stack: [Z, Y, X]
                if self.keep_3d:
                    actin_3d = img_stack.copy()
                actin_2d = self.apply_z_projection(img_stack, axis=0)
                nuclei_2d = None

        elif img_stack.ndim == 4:
            # 4D stack: [Z, C, Y, X]
            n_z, n_c, h, w = img_stack.shape
            print(f"  > 4D stack: {img_stack.shape} [Z={n_z}, C={n_c}, Y={h}, X={w}]")

            # Validate channel count
            if not self._validate_channel_count(n_c, image_path):
                return None

            # Keep 3D volumes if requested
            if self.keep_3d:
                actin_3d = img_stack[:, self.ch_actin, :, :].copy()

            # Project to 2D (MIP)
            actin_2d = self.apply_z_projection(img_stack[:, self.ch_actin, :, :], axis=0)
            if n_c > self.ch_nuclei:
                nuclei_2d = self.apply_z_projection(img_stack[:, self.ch_nuclei, :, :], axis=0)

            # Extract bacteria channel via MIP if enabled
            if self.bacterial_enabled and self.ch_bacteria is not None and n_c > self.ch_bacteria:
                bacteria_2d = self.apply_z_projection(img_stack[:, self.ch_bacteria, :, :], axis=0)
                bacteria_2d_raw = bacteria_2d.copy()

            # Extract brightfield channel via MIP if configured
            if self.ch_brightfield is not None and n_c > self.ch_brightfield:
                brightfield_2d = self.apply_z_projection(img_stack[:, self.ch_brightfield, :, :], axis=0)

        else:
            print(f'WARNING: Unexpected image dimensions: {img_stack.shape}')
            return None

        # Free the full stack
        del img_stack
        gc.collect()

        # Preserve RAW copies BEFORE brightness/contrast adjustment
        actin_2d_raw = actin_2d.copy()
        actin_3d_raw = actin_3d.copy() if actin_3d is not None else None

        # Apply brightness/contrast to actin
        actin_params = self.config['brightness_contrast']['actin']
        actin_2d = self.adjust_brightness_contrast(actin_2d, actin_params)

        if nuclei_2d is not None:
            nuclei_params = self.config['brightness_contrast']['nuclei']
            nuclei_2d = self.adjust_brightness_contrast(nuclei_2d, nuclei_params)

        # Apply B/C to 3D actin volume if present
        if actin_3d is not None:
            for z in range(actin_3d.shape[0]):
                actin_3d[z] = self.adjust_brightness_contrast(actin_3d[z], actin_params)

        # Image-level CLAHE on full actin MIP (after B/C, before per-cell skeletonization)
        # Affects both skeleton topology (Path B) AND weighted density metrics (Path C)
        if self.img_clahe_enabled:
            from skimage import exposure as _exposure
            actin_float = actin_2d.astype(np.float32) / 65535.0
            actin_float = _exposure.equalize_adapthist(
                actin_float,
                kernel_size=self.img_clahe_tile,
                clip_limit=self.img_clahe_clip,
            )
            actin_2d = (actin_float * 65535).astype(np.uint16)
            if actin_3d is not None:
                for z in range(actin_3d.shape[0]):
                    z_float = actin_3d[z].astype(np.float32) / 65535.0
                    z_float = _exposure.equalize_adapthist(
                        z_float,
                        kernel_size=self.img_clahe_tile,
                        clip_limit=self.img_clahe_clip,
                    )
                    actin_3d[z] = (z_float * 65535).astype(np.uint16)
            print(f"  Image-level CLAHE applied (clip={self.img_clahe_clip}, tile={self.img_clahe_tile})")

        # Apply B/C to bacteria channel if present
        if bacteria_2d is not None:
            bact_bc = self.bact_bc_params.get('bacteria', {
                'brightness': 1.8, 'contrast': 2.0, 'gamma': 0.85
            })
            bacteria_2d = self.adjust_brightness_contrast(bacteria_2d, bact_bc)

        # Apply B/C to brightfield if present
        if brightfield_2d is not None:
            bf_bc = self.bact_bc_params.get('brightfield', {
                'brightness': 1.0, 'contrast': 1.0, 'gamma': 1.0
            })
            brightfield_2d = self.adjust_brightness_contrast(brightfield_2d, bf_bc)

        return {
            'actin_2d': actin_2d,
            'nuclei_2d': nuclei_2d,
            'actin_2d_raw': actin_2d_raw,
            'actin_3d': actin_3d,
            'actin_3d_raw': actin_3d_raw,
            'bacteria_2d': bacteria_2d,
            'bacteria_2d_raw': bacteria_2d_raw,
            'brightfield_2d': brightfield_2d,
            'original_ndim': original_ndim
        }

    def segment_image(self, image_path):
        """
        Segment image using 2D projection with configurable cell segmentation channel.

        Returns dict with keys:
            cell_masks, nucleus_masks: 2D segmentation masks
            actin_image, nuclei_image: 2D projected (B/C-adjusted)
            actin_image_raw: Raw 2D projection for measurements
            actin_volume_3d, actin_volume_3d_raw: 3D volumes (if keep_3d)
            bacteria_2d, bacteria_2d_raw: Bacteria channel (if enabled)
            brightfield_2d: Brightfield channel (if configured)
            processing_mode: '2D' or '3D'
        Returns None on error.
        """
        preproc = self.load_and_preprocess(image_path)
        if preproc is None:
            return None

        actin_2d = preproc['actin_2d']
        nuclei_2d = preproc['nuclei_2d']
        actin_3d = preproc.get('actin_3d')

        # Determine processing mode
        processing_mode = '3D' if (actin_3d is not None and self.keep_3d) else '2D'

        # Determine which channel to use for cell boundary segmentation
        if self.cell_seg_channel == 'brightfield' and preproc.get('brightfield_2d') is not None:
            seg_channel = preproc['brightfield_2d']
            print("  > Using BRIGHTFIELD for cell segmentation (inhibitor mode)")
        else:
            seg_channel = actin_2d
            if self.cell_seg_channel == 'brightfield':
                print("  > WARNING: brightfield requested but not available, falling back to actin")

        # Cell segmentation
        if nuclei_2d is not None:
            stacked_image = np.stack([seg_channel, nuclei_2d], axis=-1)
        else:
            stacked_image = np.stack([seg_channel, seg_channel], axis=-1)

        cell_min_size = self.cell_cfg.get('min_size', -1)

        cell_masks, flows, styles = self.cell_model.eval(
            stacked_image,
            diameter=self.cell_cfg['diameter'] if self.cell_cfg['diameter'] > 0 else None,
            flow_threshold=self.cell_cfg['flow_threshold'],
            cellprob_threshold=self.cell_cfg['cellprob_threshold'],
            min_size=cell_min_size if cell_min_size > 0 else -1
        )

        # Nucleus segmentation
        if nuclei_2d is not None:
            nucleus_min_size = self.nucleus_cfg.get('min_size', -1)

            nucleus_masks, _, _ = self.nucleus_model.eval(
                nuclei_2d,
                diameter=self.nucleus_cfg['diameter'] if self.nucleus_cfg['diameter'] > 0 else None,
                flow_threshold=self.nucleus_cfg['flow_threshold'],
                cellprob_threshold=self.nucleus_cfg['cellprob_threshold'],
                min_size=nucleus_min_size if nucleus_min_size > 0 else -1
            )
        else:
            nucleus_masks = np.zeros_like(cell_masks)

        # Save masks if requested
        if self.config['processing']['save_masks']:
            output_dir = Path(self.config['output_dir']) / 'masks'
            image_name = Path(image_path).stem
            output_dir.mkdir(parents=True, exist_ok=True)
            tifffile.imwrite(output_dir / f'{image_name}_cell_mask.tif', cell_masks.astype(np.uint16))
            tifffile.imwrite(output_dir / f'{image_name}_nucleus_mask.tif', nucleus_masks.astype(np.uint16))

        # Build result dictionary (actin_image always = actin channel, not seg channel)
        result_dict = {
            'cell_masks': cell_masks,
            'nucleus_masks': nucleus_masks,
            'actin_image': actin_2d,
            'nuclei_image': nuclei_2d,
            'actin_image_raw': preproc['actin_2d_raw'],
            'processing_mode': processing_mode
        }

        # Add 3D volumes if available
        if processing_mode == '3D':
            result_dict['actin_volume_3d'] = actin_3d
            result_dict['actin_volume_3d_raw'] = preproc.get('actin_3d_raw')

        # Add bacteria data if available
        if preproc.get('bacteria_2d') is not None:
            result_dict['bacteria_2d'] = preproc['bacteria_2d']
            result_dict['bacteria_2d_raw'] = preproc['bacteria_2d_raw']

        # Add brightfield if available
        if preproc.get('brightfield_2d') is not None:
            result_dict['brightfield_2d'] = preproc['brightfield_2d']

        return result_dict
