"""
TopoActin V3 Main Pipeline - Dual-Path Processing with Weighted Skeleton Analysis
Handles both 2D and 3D Z-stack data with automatic detection.
INCLUDES: Position-based downsampling + Weighted skeleton analysis
INCLUDES: Optional bacterial burden quantification (config-driven)
"""

import warnings
warnings.filterwarnings('ignore')

import json
import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from tqdm import tqdm
from datetime import datetime
import gc
import argparse

# Import V3 modules - UPDATED FOR WEIGHTED SKELETON
try:
    from segmentation_v3 import CellSegmenter
    from filtering import NucleiFilter, CellFilter
    from skeletonization_v3 import CellSkeletonizer  # V3 WEIGHTED: Sato-weighted skeleton analysis
    from topological_analysis_v3 import TopologicalAnalyzer
except ImportError as e:
    print(f"❌ CRITICAL ERROR: Could not import core modules: {e}")
    print("Ensure segmentation_v3.py, filtering.py, skeletonization_v3_weighted.py, etc. are in the same folder.")
    sys.exit(1)

# Import visualization modules
try:
    from visualization import ResultsVisualizer
    VISUALIZATION_AVAILABLE = True
except ImportError:
    VISUALIZATION_AVAILABLE = False
    print('ℹ️  visualization.py not found - skipping per-cell visualizations')

try:
    from experiment_level_visualizer import ExperimentLevelVisualizer
    EXPERIMENT_VIZ_AVAILABLE = True
except ImportError:
    EXPERIMENT_VIZ_AVAILABLE = False
    print('experiment_level_visualizer.py not found - skipping experiment summaries')

# Optional bacterial burden integration
try:
    from bacterial_burden.detection import quantify_bacteria_per_cell, compute_zone_bacterial_density
    BACTERIAL_AVAILABLE = True
except ImportError:
    BACTERIAL_AVAILABLE = False


def parse_well_position(filename):
    """Parse well position from filename like B2_F01 or B10_F03"""
    # Try specific format first: B2_F01, B10_F03, etc.
    match = re.search(r'([A-H])(\d{1,2})_F\d+', filename)
    if match:
        row = match.group(1)
        col_str = match.group(2)
        col_num = int(col_str)
        return f"{row}{col_num}", col_num
    
    # Fallback: generic pattern with word boundary
    match = re.search(r'([A-H])(\d{1,2})\b', filename)
    if match:
        row = match.group(1)
        col_str = match.group(2)
        col_num = int(col_str)
        return f"{row}{col_num}", col_num
    
    return None, None


def assign_group(image_name, well_to_group):
    well_position, col_num = parse_well_position(image_name)
    if well_position is None or col_num is None:
        return 'N/A', 'Unassigned'
    group_name = well_to_group.get(str(col_num), 'Unassigned')
    return well_position, group_name


def load_config(config_path):
    config_file = Path(config_path)
    if not config_file.exists():
        print(f'❌ Config file not found: {config_file}')
        sys.exit(1)
    with open(config_file, 'r') as f:
        return json.load(f)


def get_image_files(input_dir):
    input_path = Path(input_dir)
    print(f"🔎 Searching for images in: {input_path.absolute()}")
    
    if not input_path.exists():
        print(f"❌ ERROR: Input directory does not exist!")
        return []
    
    files = sorted([f for f in input_path.rglob('*') if f.suffix.lower() in ('.tif', '.tiff')])
    if len(files) == 0:
        print(f"⚠️  Directory exists but contains no .tif/.tiff files.")
    return files


def apply_position_downsampling_if_enabled(file_list, config):
    """
    Apply position-based downsampling if enabled in config.
    
    Parameters
    ----------
    file_list : list
        Original file list
    config : dict
        Full configuration dict
        
    Returns
    -------
    list
        Downsampled file list (or original if disabled)
    """
    downsampling_config = config.get('downsampling', {'enabled': False})
    
    if not downsampling_config.get('enabled', False):
        return file_list
    
    try:
        from position_downsampling import apply_position_downsampling as downsample_func
        
        num_positions = downsampling_config.get('num_positions', 10)
        seed = downsampling_config.get('seed', 42)
        
        print(f"\n{'='*80}")
        print(f"POSITION DOWNSAMPLING ENABLED")
        print(f"{'='*80}")
        print(f"Target positions per well: {num_positions}")
        print(f"Random seed: {seed}")
        
        downsampled = downsample_func(file_list, num_positions, seed, verbose=True)
        
        if len(downsampled) == 0:
            print("Warning: Downsampling returned no files, using original list")
            return file_list
        
        return downsampled
        
    except ImportError:
        print("\nWarning: position_downsampling.py not found, skipping downsampling")
        return file_list
    except Exception as e:
        print(f"\nWarning: Downsampling failed ({str(e)}), using all files")
        return file_list


def process_experiment(config, experiment_config=None):
    # Build effective config with per-experiment overrides
    effective_config = config.copy()

    if experiment_config:
        experiment_name = experiment_config['name']
        input_dir = experiment_config['input_dir']
        output_dir = experiment_config['output_dir']
        if 'group_mapping' in experiment_config:
            well_to_group = experiment_config['group_mapping']['well_to_group']
            group_colors = experiment_config['group_mapping']['group_colors']
        else:
            well_to_group = config['global_well_mapping']['well_to_group']
            group_colors = config['global_well_mapping']['group_colors']
        effective_config['output_dir'] = output_dir
        effective_config['input_dir'] = input_dir
        # Per-experiment bacterial_burden override
        if 'bacterial_burden' in experiment_config:
            effective_config['bacterial_burden'] = {
                **config.get('bacterial_burden', {}),
                **experiment_config['bacterial_burden']
            }
    else:
        experiment_name = 'Single_Experiment'
        input_dir = config['input_dir']
        output_dir = config['output_dir']
        well_to_group = config.get('well_to_group', {})
        group_colors = config.get('group_colors', {})

    # Determine if bacterial burden is active for this experiment
    bact_config = effective_config.get('bacterial_burden', {})
    bacterial_enabled = bact_config.get('enabled', False) and BACTERIAL_AVAILABLE
    if bact_config.get('enabled', False) and not BACTERIAL_AVAILABLE:
        print("WARNING: bacterial_burden enabled but detection module not found. Disabling.")
    
    print('\n' + '='*80)
    print(f'PROCESSING EXPERIMENT: {experiment_name}')
    print(f'TopoActin V3 - 3D Dual-Path Analysis with Weighted Skeleton')
    print('='*80)
    print(f'Input:  {input_dir}')
    print(f'Output: {output_dir}')
    
    output_path = Path(output_dir)
    (output_path / 'results').mkdir(parents=True, exist_ok=True)
    (output_path / 'masks').mkdir(parents=True, exist_ok=True)
    (output_path / 'skeletons').mkdir(parents=True, exist_ok=True)
    (output_path / 'comparisons').mkdir(parents=True, exist_ok=True)
    
    image_files = get_image_files(input_dir)
    if not image_files:
        return None
    
    # Apply position downsampling if enabled
    image_files = apply_position_downsampling_if_enabled(image_files, config)
    if not image_files:
        print('No files remaining after downsampling')
        return None
    
    print(f'\n✓ Found {len(image_files)} images')
    
    segmenter = CellSegmenter(effective_config)  # Load Cellpose models ONCE
    nuclei_filter = NucleiFilter(effective_config)
    cell_filter = CellFilter(effective_config)
    skeletonizer = CellSkeletonizer(effective_config)
    visualizer = ResultsVisualizer(effective_config) if VISUALIZATION_AVAILABLE else None

    if bacterial_enabled:
        print(f'Bacterial burden: ENABLED')
        bact_detection_params = bact_config.get('bacteria_detection', {})
        bact_norm_params = bact_config.get('per_cell_normalization', {})
    
    stats = {
        'images_processed': 0, 'images_skipped': 0,
        'total_cells_detected': 0, 'cells_analyzed': 0,
        'rejected_nuclei': 0, 'rejected_cells': 0, 'rescued_cells': 0,
        'cells_3d_processed': 0, 'cells_2d_processed': 0
    }
    
    all_results = []
    pbar = tqdm(image_files, desc=f'Processing {experiment_name}', unit="img")
    
    for image_path in pbar:
        image_name = image_path.stem
        try:
            seg_result = segmenter.segment_image(image_path)
            if seg_result is None:
                stats['images_skipped'] += 1
                continue
            
            cell_masks = seg_result['cell_masks']
            nucleus_masks = seg_result['nucleus_masks']
            actin_image = seg_result['actin_image']
            nuclei_image = seg_result['nuclei_image']
            actin_volume_3d = seg_result.get('actin_volume_3d', None)
            # Raw (unadjusted) data for B/C/G-independent measurements
            actin_image_raw = seg_result.get('actin_image_raw', actin_image)
            actin_volume_3d_raw = seg_result.get('actin_volume_3d_raw', actin_volume_3d)
            
            stats['total_cells_detected'] += (len(np.unique(cell_masks)) - 1)
            
            # Filtering
            valid_nuclei, rescued_nucleus_ids = [], []
            if nuclei_image is not None:
                n_res = nuclei_filter.filter_nuclei(nucleus_masks, nuclei_image, image_name)
                valid_nuclei = n_res['valid_nuclei']
                stats['rejected_nuclei'] += len(n_res['excluded_nuclei'])
                if 'rescued_nuclei' in n_res:
                    rescued_nucleus_ids = [r['nucleus_id'] for r in n_res['rescued_nuclei']]
                    stats['rescued_cells'] += len(rescued_nucleus_ids)
            
            c_res = cell_filter.filter_cells(cell_masks, nucleus_masks, valid_nuclei, image_name)
            valid_cells = c_res['valid_cells']
            stats['rejected_cells'] += len(c_res['excluded_cells'])
            
            if not valid_cells:
                stats['images_skipped'] += 1
                continue

            stats['cells_analyzed'] += len(valid_cells)
            well_position, group_name = assign_group(image_name, well_to_group)

            # --- BACTERIAL BURDEN STEP (optional) ---
            bacteria_results = {}
            if bacterial_enabled and seg_result.get('bacteria_2d') is not None:
                valid_ids = [c['cell_id'] for c in valid_cells]
                bacteria_results = quantify_bacteria_per_cell(
                    bacteria_2d=seg_result['bacteria_2d'],
                    bacteria_raw=seg_result.get('bacteria_2d_raw', seg_result['bacteria_2d']),
                    cell_masks=cell_masks,
                    valid_cell_ids=valid_ids,
                    detection_params=bact_detection_params,
                    norm_params=bact_norm_params
                )

            for cell_info in valid_cells:
                cell_id = cell_info['cell_id']
                cell_mask = (cell_masks == cell_id)
                cell_info['rescued'] = cell_info['nucleus_id'] in rescued_nucleus_ids
                
                skel_res = skeletonizer.skeletonize_cell(
                    cell_mask=cell_mask, actin_image=actin_image,
                    image_name=image_name, cell_id=cell_id,
                    nucleus_mask=nucleus_masks, actin_volume_3d=actin_volume_3d,
                    actin_image_raw=actin_image_raw,
                    actin_volume_3d_raw=actin_volume_3d_raw
                )
                
                if skel_res:
                    feats = skel_res['features']
                    processing_mode = skel_res.get('processing_mode', '2D')
                    
                    if processing_mode == '3D': stats['cells_3d_processed'] += 1
                    else: stats['cells_2d_processed'] += 1
                    
                    # --- CALCULATE MISSING METRICS (Tortuosity) ---
                    mean_tortuosity = np.nan
                    graph_summary = skel_res.get('graph_summary')
                    if graph_summary is not None and len(graph_summary) > 0:
                        # Tortuosity = Branch Path Length / Euclidean Distance
                        euclidean = graph_summary['euclidean-distance']
                        path_len = graph_summary['branch-distance']
                        valid_branches = euclidean > 0
                        if valid_branches.any():
                            tortuosity = path_len[valid_branches] / euclidean[valid_branches]
                            mean_tortuosity = tortuosity.mean()

                    row_data = {
                        'experiment': experiment_name, 'image_name': image_name,
                        'cell_id': cell_id, 'group': group_name,
                        'processing_mode': processing_mode,
                        'measurement_source': 'raw',  # Measurements use raw (B/C/G-independent) data
                        'mean_tortuosity_edge': mean_tortuosity,
                        **cell_info['properties'],
                        **feats
                    }
                    
                    if skel_res.get('cortical_features'): row_data.update(skel_res['cortical_features'])
                    if skel_res.get('topological_features'): row_data.update(skel_res['topological_features'])

                    # Merge bacterial burden features
                    if cell_id in bacteria_results:
                        bact = bacteria_results[cell_id]
                        row_data['intracellular_bact_area'] = bact['bacterial_area']
                        row_data['area_ratio'] = round(bact['area_ratio'], 6)
                        row_data['has_bacteria'] = bact['bacterial_area'] > 0

                        # Per-zone bacterial density (perinuclear / middle / cortical)
                        bact_mask = bact.get('bacteria_mask')
                        if bact_mask is not None:
                            nuc_id = cell_info.get('nucleus_id')
                            nuc_mask_cell = (nucleus_masks == nuc_id).astype(bool) if nuc_id else None
                            zone_thresh = config.get('topological_analysis', {})
                            zone_thresholds = {
                                'inner': zone_thresh.get('inner_threshold', 0.33),
                                'mid':   zone_thresh.get('mid_threshold',   0.66)
                            }
                            zone_bact = compute_zone_bacterial_density(
                                bacteria_mask=bact_mask,
                                cell_mask=cell_mask,
                                nucleus_mask=nuc_mask_cell,
                                zone_thresholds=zone_thresholds
                            )
                            row_data.update(zone_bact)
                    elif bacterial_enabled:
                        row_data['intracellular_bact_area'] = np.nan
                        row_data['area_ratio'] = np.nan
                        row_data['has_bacteria'] = np.nan
                        for col in ['bact_area_perinuclear', 'bact_area_middle', 'bact_area_cortical',
                                    'bact_density_perinuclear', 'bact_density_middle', 'bact_density_cortical']:
                            row_data[col] = np.nan

                    all_results.append(row_data)
                    
                    if visualizer and len(all_results) <= 20:
                        visualizer.create_comparison(cell_mask, actin_image, skel_res['skeleton'], image_name, cell_id)
            
            stats['images_processed'] += 1
            acc_rate = (stats['cells_analyzed'] / stats['total_cells_detected'] * 100) if stats['total_cells_detected'] else 0
            pbar.set_postfix({"Analyzed": stats['cells_analyzed'], "Yield": f"{acc_rate:.1f}%"})
            gc.collect()
            
        except Exception as e:
            tqdm.write(f"❌ Error {image_name}: {e}")
            stats['images_skipped'] += 1
            continue

    # --- SAVE RESULTS ---
    print('\n' + '='*80)
    print('📊 PIPELINE ACCURACY REPORT')
    print('='*80)
    print(f"Total Detected: {stats['total_cells_detected']}")
    print(f"Total Analyzed: {stats['cells_analyzed']}")
    print(f"Rejection Rate: {((stats['rejected_nuclei']+stats['rejected_cells'])/stats['total_cells_detected']*100 if stats['total_cells_detected'] else 0):.1f}%")
    print(f"Rescued Cells:  {stats['rescued_cells']}")
    print('='*80)
    
    if not all_results:
        print('❌ No results generated.')
        return None
        
    df = pd.DataFrame(all_results)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    csv_file = output_path / 'results' / f'cytoskeleton_analysis_v3_{experiment_name}_{timestamp}.csv'
    df.to_csv(csv_file, index=False)
    print(f'\n✓ CSV Saved: {csv_file}')
    print(f'📊 Total columns: {len(df.columns)} (36 existing + 9 weighted skeleton metrics)')
    
    # --- VISUALIZATION COMPATIBILITY FIX ---
    if EXPERIMENT_VIZ_AVAILABLE and len(df) > 10:
        print('\nGenerating Charts...')
        
        # 1. Map V3 column names to V2 Visualizer names
        # V3 uses 'spectral_entropy', V2 Visualizer wants 'cortical_spectral_entropy'
        if 'spectral_entropy' in df.columns:
            df['cortical_spectral_entropy'] = df['spectral_entropy']
            
        # V3 uses 'zone_3_lacunarity' (outer), V2 Visualizer wants 'hole_size_edge'
        if 'zone_3_lacunarity' in df.columns:
            df['hole_size_edge'] = df['zone_3_lacunarity']
            
        # Ensure tortuosity is available (we just calculated it as 'mean_tortuosity_edge')
        
        viz_config = config.copy()
        viz_config['group_colors'] = group_colors
        exp_viz = ExperimentLevelVisualizer(output_path, viz_config)
        
        try:
            summary_path = exp_viz.generate_experiment_summary(df, experiment_name, timestamp)
            print(f'✓ Charts Saved: {summary_path}')
        except Exception as e:
            print(f'⚠️ Chart generation failed: {e}')
            import traceback
            traceback.print_exc()

    return {'df': df, 'csv': csv_file}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='config_batch.json')
    args = parser.parse_args()
    
    config = load_config(args.config)
    print("🚀 Starting TopoActin V3 with Weighted Skeleton Analysis...")
    
    if config.get('batch_mode'):
        for exp in config['experiments']:
            if exp.get('enabled', True): process_experiment(config, exp)
    else:
        process_experiment(config)
    
    print("\n✅ Done.")

if __name__ == '__main__':
    main()