"""
Bacterial Burden Quantification Pipeline
Standalone pipeline for host-pathogen interaction analysis.
Quantifies intracellular bacterial area per host cell in 96-well plate HCS assays.

Usage:
    python main_bacterial_burden.py --config config_bacterial.json
"""

import warnings
warnings.filterwarnings('ignore')

import json
import sys
import re
import os
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from datetime import datetime
import gc
import argparse

# Add parent directory to path for importing shared utilities
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Import bacterial burden modules (same directory)
from segmentation_bacterial import BacterialSegmenter
from filtering_bacterial import HostCellFilter
from quantification_bacterial import BacterialBurdenQuantifier

# Reuse utilities from existing TopoActin pipeline
from main_pipeline_v3 import load_config, get_image_files
from position_downsampling import apply_position_downsampling


# ---------------------------------------------------------------------------
# Well Parsing & Treatment Assignment
# ---------------------------------------------------------------------------

def parse_well_position_robust(filename):
    """
    Parse well position from filename with strict regex to prevent
    B1 matching inside B11 or B12.

    Uses negative lookahead (?![0-9]) to ensure the digit group is not
    followed by another digit. Handles B01 as Column 1 via int() conversion.

    Priority pattern: B2_F01, B10_F03 (anchored on _F separator)
    Fallback: B2, B10 with negative lookahead boundary

    Parameters
    ----------
    filename : str or Path
        Image filename (stem or full path)

    Returns
    -------
    tuple (well_str, col_int)
        well_str: e.g. "B2", col_int: e.g. 2
        Returns (None, None) if parsing fails
    """
    name = Path(filename).stem if not isinstance(filename, str) else str(filename)

    # Priority: explicit format with _F field separator
    match = re.search(r'([A-H])(\d{1,2})_F\d+', name)
    if match:
        row = match.group(1)
        col_num = int(match.group(2))
        return f"{row}{col_num}", col_num

    # Fallback: negative lookahead prevents B1 matching inside B11
    match = re.search(r'([A-H])(\d{1,2})(?![0-9])', name)
    if match:
        row = match.group(1)
        col_num = int(match.group(2))
        return f"{row}{col_num}", col_num

    return None, None


def assign_treatment(image_name, plate_layout):
    """
    Map well position to treatment name using the plate_layout dictionary.

    Parameters
    ----------
    image_name : str
        Image filename or stem
    plate_layout : dict
        Maps column number (as string key) to treatment name

    Returns
    -------
    tuple (well_position_str, treatment_name)
    """
    well_pos, col_num = parse_well_position_robust(image_name)
    if well_pos is None or col_num is None:
        return 'N/A', 'Unassigned'

    # Config stores keys as strings; also try int key
    treatment = plate_layout.get(str(col_num),
                plate_layout.get(col_num, 'Unassigned'))
    return well_pos, treatment


# ---------------------------------------------------------------------------
# Downsampling Integration
# ---------------------------------------------------------------------------

def apply_downsampling_if_enabled(file_list, config):
    """
    Apply position-based downsampling if enabled in config.
    Reuses the position_downsampling module from the parent project.
    """
    ds_config = config.get('downsampling', {'enabled': False})
    if not ds_config.get('enabled', False):
        return file_list

    num_positions = ds_config.get('num_positions', 10)
    seed = ds_config.get('seed', 42)

    print(f'\nPOSITION DOWNSAMPLING ENABLED')
    print(f'Target positions per well: {num_positions}, seed: {seed}')

    try:
        downsampled = apply_position_downsampling(file_list, num_positions, seed, verbose=True)
        if len(downsampled) == 0:
            print('Warning: Downsampling returned no files, using original list')
            return file_list
        return downsampled
    except Exception as e:
        print(f'Warning: Downsampling failed ({e}), using all files')
        return file_list


# ---------------------------------------------------------------------------
# Main Processing
# ---------------------------------------------------------------------------

def process_experiment(config, experiment_config=None):
    """
    Main processing loop for bacterial burden quantification.

    For each image:
    1. Segment host cells, nuclei, and bacteria
    2. QC check on all 3 channels
    3. Filter host cells (edge, cell morphology, nucleus association, nucleus morphology)
    4. Quantify intracellular bacterial area per valid cell
    5. Record per-cell measurements with well/treatment metadata

    Output:
    - Results_SingleCell_{name}_{timestamp}.csv
    - Results_WellSummary_{name}_{timestamp}.csv
    """
    # --- Setup ---
    if experiment_config:
        experiment_name = experiment_config['name']
        input_dir = experiment_config['input_dir']
        output_dir = experiment_config['output_dir']
        plate_layout = experiment_config.get('plate_layout', config.get('plate_layout', {}))
    else:
        experiment_name = 'BacterialBurden'
        input_dir = config['input_dir']
        output_dir = config['output_dir']
        plate_layout = config.get('plate_layout', {})

    print('\n' + '=' * 80)
    print(f'BACTERIAL BURDEN QUANTIFICATION: {experiment_name}')
    print('=' * 80)
    print(f'Input:  {input_dir}')
    print(f'Output: {output_dir}')

    output_path = Path(output_dir)
    (output_path / 'results').mkdir(parents=True, exist_ok=True)
    if config['processing'].get('save_masks', False):
        (output_path / 'masks').mkdir(parents=True, exist_ok=True)

    # --- Find images ---
    image_files = get_image_files(input_dir)
    if not image_files:
        print('No images found.')
        return None

    image_files = apply_downsampling_if_enabled(image_files, config)
    if not image_files:
        print('No files remaining after downsampling.')
        return None

    print(f'Images to process: {len(image_files)}')

    # --- Print plate layout ---
    print('\nPlate Layout:')
    for col, treatment in sorted(plate_layout.items(), key=lambda x: int(x[0])):
        print(f'  Column {col}: {treatment}')

    # --- Initialize modules ---
    segmenter = BacterialSegmenter(config)
    cell_filter = HostCellFilter(config)
    quantifier = BacterialBurdenQuantifier(config)

    # --- Processing loop ---
    stats = {
        'images_processed': 0,
        'images_skipped': 0,
        'qc_failed': 0,
        'total_hosts_detected': 0,
        'hosts_analyzed': 0,
        'hosts_excluded': 0,
    }

    all_results = []
    qc_log = []

    pbar = tqdm(image_files, desc=f'Bacterial Burden', unit='img')

    for image_path in pbar:
        image_name = image_path.stem if hasattr(image_path, 'stem') else Path(image_path).stem

        try:
            # 1. Segment
            seg_result = segmenter.segment_image(image_path)
            if seg_result is None:
                stats['images_skipped'] += 1
                continue

            # 2. QC check
            qc = seg_result.get('qc_results', {})
            qc_entry = {
                'image_name': image_name,
                'z_best': seg_result.get('z_best', -1),
            }
            any_fail = False
            for ch_name, ch_qc in qc.items():
                qc_entry[f'{ch_name}_pass'] = ch_qc['pass']
                qc_entry[f'{ch_name}_saturated'] = ch_qc.get('saturated', False)
                qc_entry[f'{ch_name}_max'] = ch_qc['max_intensity']
                qc_entry[f'{ch_name}_mean'] = ch_qc['mean_intensity']
                if not ch_qc['pass']:
                    any_fail = True
            qc_log.append(qc_entry)

            if any_fail:
                stats['qc_failed'] += 1
                if config['qc'].get('skip_empty_images', True):
                    stats['images_skipped'] += 1
                    continue

            host_masks = seg_result['host_masks']
            nucleus_masks = seg_result['nucleus_masks']

            n_hosts = len(np.unique(host_masks)) - 1  # exclude background (0)
            stats['total_hosts_detected'] += n_hosts

            # 3. Filter
            filter_result = cell_filter.filter_host_cells(host_masks, nucleus_masks)
            valid_cells = filter_result['valid_cells']
            stats['hosts_excluded'] += filter_result['exclusion_stats']['excluded_count']

            if not valid_cells:
                stats['images_skipped'] += 1
                continue

            # 4. Quantify
            well_pos, treatment = assign_treatment(image_name, plate_layout)

            bacteria_area_results = getattr(segmenter, '_bacteria_area_results', {})

            cell_measurements = quantifier.quantify_image(
                filter_result['filtered_host_masks'],
                valid_cells,
                bacteria_area_results
            )

            # 5. Collect results
            for measurement in cell_measurements:
                measurement['experiment'] = experiment_name
                measurement['image_name'] = image_name
                measurement['well'] = well_pos
                measurement['treatment'] = treatment
                measurement['z_best'] = seg_result.get('z_best', -1)
                all_results.append(measurement)

            stats['hosts_analyzed'] += len(cell_measurements)
            stats['images_processed'] += 1

            yield_pct = (stats['hosts_analyzed'] / max(stats['total_hosts_detected'], 1)) * 100
            pbar.set_postfix({
                'Hosts': stats['hosts_analyzed'],
                'Yield': f'{yield_pct:.0f}%'
            })
            gc.collect()

        except Exception as e:
            tqdm.write(f'Error processing {image_name}: {e}')
            stats['images_skipped'] += 1
            continue

    # --- Report ---
    print('\n' + '=' * 80)
    print('PIPELINE REPORT')
    print('=' * 80)
    print(f"Images processed:            {stats['images_processed']}")
    print(f"Images skipped:              {stats['images_skipped']}")
    print(f"QC failures:                 {stats['qc_failed']}")
    print(f"Host cells detected:         {stats['total_hosts_detected']}")
    print(f"Host cells analyzed:         {stats['hosts_analyzed']}")
    print(f"Host cells excluded:         {stats['hosts_excluded']}")
    exc_rate = stats['hosts_excluded'] / max(stats['total_hosts_detected'], 1) * 100
    print(f"Exclusion rate:              {exc_rate:.1f}%")
    print('=' * 80)

    if not all_results:
        print('No results generated.')
        return None

    # --- Save results ---
    df = pd.DataFrame(all_results)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    # Reorder columns for readability
    col_order = [
        'experiment', 'image_name', 'well', 'treatment',
        'cell_id', 'nucleus_id', 'z_best',
        'host_area', 'cell_area', 'cell_solidity', 'cell_eccentricity', 'cell_perimeter',
        'nucleus_area', 'nucleus_solidity', 'nucleus_eccentricity', 'nuc_cell_ratio',
        'intracellular_bact_area', 'area_ratio'
    ]
    # Only include columns that exist
    col_order = [c for c in col_order if c in df.columns]
    remaining = [c for c in df.columns if c not in col_order]
    df = df[col_order + remaining]

    single_cell_csv = output_path / 'results' / f'Results_SingleCell_{experiment_name}_{timestamp}.csv'
    df.to_csv(single_cell_csv, index=False)
    print(f'\nSingle-cell CSV: {single_cell_csv}')
    print(f'  Rows: {len(df)}, Columns: {len(df.columns)}')

    # --- Well summary ---
    summary_df = df.groupby(['well', 'treatment']).agg(
        total_host_count=('cell_id', 'count'),
        mean_area_ratio=('area_ratio', 'mean'),
        std_area_ratio=('area_ratio', 'std'),
        median_area_ratio=('area_ratio', 'median'),
        mean_host_area=('host_area', 'mean'),
        total_intracellular_bact_area=('intracellular_bact_area', 'sum'),
    ).reset_index()

    well_csv = output_path / 'results' / f'Results_WellSummary_{experiment_name}_{timestamp}.csv'
    summary_df.to_csv(well_csv, index=False)
    print(f'Well summary CSV: {well_csv}')
    print(f'  Rows: {len(summary_df)}, Columns: {len(summary_df.columns)}')

    # --- QC log ---
    if qc_log:
        qc_df = pd.DataFrame(qc_log)
        qc_csv = output_path / 'results' / f'QC_Log_{experiment_name}_{timestamp}.csv'
        qc_df.to_csv(qc_csv, index=False)
        print(f'QC log CSV: {qc_csv}')

    return {
        'df': df,
        'summary_df': summary_df,
        'single_cell_csv': single_cell_csv,
        'well_csv': well_csv,
        'stats': stats
    }


def main():
    parser = argparse.ArgumentParser(
        description='Bacterial Burden Quantification Pipeline'
    )
    parser.add_argument(
        '--config', default='config_bacterial.json',
        help='Path to bacterial burden config JSON (default: config_bacterial.json)'
    )
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.is_absolute():
        # Look in the script's directory first
        script_dir = Path(__file__).resolve().parent
        if (script_dir / config_path).exists():
            config_path = script_dir / config_path

    config = load_config(str(config_path))
    print('Starting Bacterial Burden Quantification Pipeline...')
    print(f'Config: {config_path}')

    if config.get('batch_mode'):
        for exp in config['experiments']:
            if exp.get('enabled', True):
                process_experiment(config, exp)
    else:
        process_experiment(config)

    print('\nDone.')


if __name__ == '__main__':
    main()
