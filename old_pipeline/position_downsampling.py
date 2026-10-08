"""
Position-Based Downsampling Module
Randomly selects N positions from each well for faster analysis
while maintaining statistical validity.

Example:
    If each well has 26 positions (r01c01 - r13c02), we can select
    10 random positions and analyze only those across ALL wells.
    This maintains comparability while reducing processing time.
"""

import re
import random
from pathlib import Path
from collections import defaultdict
import numpy as np


def parse_filename_position(filename):
    """
    Extract position information from filename.
    
    Handles various naming conventions, including:
    - B2_F01_Z_stack_merged.tif (well=B2, position=F01)
    - B10_F03_Z_stack_merged.tif (well=B10, position=F03)
    - r01c01, r02c03 (row/column grid)
    - s01, s02, s03 (site/position number)
    - f001, f002 (field number)
    - p01, p02 (position number)
    
    Parameters
    ----------
    filename : str or Path
        Image filename
        
    Returns
    -------
    tuple
        (well_id, position_id) or (None, None) if cannot parse
        
    Examples
    --------
    >>> parse_filename_position('B2_F01_Z_stack_merged.tif')
    ('B2', 'F01')
    
    >>> parse_filename_position('B10_F03_Z_stack_merged.tif')
    ('B10', 'F03')
    """
    stem = Path(filename).stem
    
    # PRIORITY PATTERN: User's specific format (B2_F01_Z_stack_merged)
    # Pattern: [Letter][Number]_[F][Number]_
    user_pattern = r'^([A-H]\d{1,2})_([F]\d{2,3})_'
    match = re.search(user_pattern, stem)
    if match:
        well_id = match.group(1)  # B2, B10, C5, etc.
        position_id = match.group(2)  # F01, F02, F03, etc.
        return well_id, position_id
    
    # Try common patterns for position
    position_patterns = [
        # Pattern: r01c01 (row/column)
        r'r(\d{2})c(\d{2})',
        # Pattern: s01, s001 (site)
        r's(\d{1,3})',
        # Pattern: f001, f01 (field) - lowercase
        r'f(\d{1,3})',
        # Pattern: F001, F01 (field) - uppercase
        r'F(\d{1,3})',
        # Pattern: p01, p001 (position)
        r'p(\d{1,3})',
    ]
    
    # Try common patterns for well
    well_patterns = [
        # Pattern: B2, B10, C5 (single letter + 1-2 digits)
        r'\b([A-H]\d{1,2})\b',
        # Pattern: well_A01, well_B03
        r'well[_-]?([A-H]\d{2})',
        # Pattern: f01, fld01 (field/well ID)
        r'f(?:ld)?(\d{1,3})',
        # Pattern: w01, w1 (well number)
        r'w(\d{1,2})',
    ]
    
    # Extract position
    position_id = None
    for pattern in position_patterns:
        match = re.search(pattern, stem, re.IGNORECASE)
        if match:
            if len(match.groups()) == 2:  # r01c01 format
                position_id = f"r{match.group(1)}c{match.group(2)}"
            else:  # s01, f01, p01 format
                position_id = f"{match.group(0)}"
            break
    
    # Extract well
    well_id = None
    for pattern in well_patterns:
        match = re.search(pattern, stem, re.IGNORECASE)
        if match:
            well_id = match.group(1)
            break
    
    # If no explicit well found, try to infer from beginning of filename
    if well_id is None:
        # Sometimes well is just the first part of the filename
        parts = stem.split('_')
        if len(parts) > 0 and len(parts[0]) <= 4:
            well_id = parts[0]
    
    return well_id, position_id


def group_files_by_well_and_position(file_list):
    """
    Group files by well and position.
    
    Parameters
    ----------
    file_list : list of str/Path
        List of image file paths
        
    Returns
    -------
    dict
        Nested dict: {well_id: {position_id: [files]}}
        
    Example
    -------
    >>> files = ['well_A01_s01_z001.tif', 'well_A01_s01_z002.tif', 
    ...          'well_A01_s02_z001.tif']
    >>> grouped = group_files_by_well_and_position(files)
    >>> grouped['A01']['s01']
    ['well_A01_s01_z001.tif', 'well_A01_s01_z002.tif']
    """
    grouped = defaultdict(lambda: defaultdict(list))
    unparsed = []
    
    for file_path in file_list:
        well_id, position_id = parse_filename_position(file_path)
        
        if well_id and position_id:
            grouped[well_id][position_id].append(file_path)
        else:
            unparsed.append(file_path)
    
    if unparsed:
        print(f"Warning: Could not parse position info for {len(unparsed)} files")
        print(f"Example unparsed: {unparsed[0]}")
    
    return dict(grouped)


def select_random_positions(grouped_files, num_positions, seed=42):
    """
    Randomly select N positions to use across ALL wells.
    
    Strategy:
    1. Find ALL unique positions across all wells (UNION, not intersection)
    2. Randomly select N positions from that set
    3. Apply SAME N positions to all wells (each well uses the positions it has)
    4. ALL WELLS ARE PROCESSED (not filtered out!)
    
    Parameters
    ----------
    grouped_files : dict
        Nested dict from group_files_by_well_and_position
        Format: {well_id: {position_id: [file_paths]}}
    num_positions : int
        Number of positions to select
    seed : int
        Random seed for reproducibility
        
    Returns
    -------
    list
        Selected position IDs (e.g., ['F01', 'F03', 'F07', ...])
        
    Example
    -------
    >>> grouped = {
    ...     'B2': {'F01': [...], 'F02': [...], 'F03': [...]},
    ...     'B3': {'F01': [...], 'F02': [...], 'F04': [...]}
    ... }
    >>> select_random_positions(grouped, num_positions=2, seed=42)
    ['F01', 'F03']
    # Then both B2 and B3 will be processed, using only F01 and F03 positions
    """
    random.seed(seed)
    
    if not grouped_files:
        print("Warning: No files to downsample")
        return []
    
    # Get ALL unique positions across all wells (UNION, not intersection!)
    all_positions = set()
    for well_id, positions in grouped_files.items():
        all_positions.update(positions.keys())
    
    all_positions = sorted(list(all_positions))
    
    print(f"Info: Found {len(all_positions)} unique positions across {len(grouped_files)} wells")
    
    # Show position distribution
    wells_per_position = {}
    for pos in all_positions:
        count = sum(1 for well_positions in grouped_files.values() if pos in well_positions)
        wells_per_position[pos] = count
    
    # Select random positions
    if len(all_positions) <= num_positions:
        selected = all_positions
        print(f"Warning: Requested {num_positions} positions but only {len(all_positions)} available")
        print(f"Info: Using all {len(all_positions)} positions")
    else:
        selected = random.sample(all_positions, num_positions)
        print(f"Info: Randomly selected {num_positions} positions from {len(all_positions)} available")
    
    # Show how many wells have each selected position
    coverage = [wells_per_position[pos] for pos in selected]
    print(f"Info: Selected positions appear in {min(coverage)}-{max(coverage)} wells (min-max)")
    
    return sorted(selected)


def apply_position_downsampling(file_list, num_positions, seed=42, verbose=True):
    """
    Apply position-based downsampling to file list.
    
    This function:
    1. Groups files by well and position
    2. Randomly selects N positions (same for all wells)
    3. Returns only files from those positions
    
    Parameters
    ----------
    file_list : list of str/Path
        Input file paths
    num_positions : int
        Number of positions to keep
    seed : int
        Random seed for reproducibility
    verbose : bool
        Print detailed statistics
        
    Returns
    -------
    list
        Filtered file paths (only selected positions)
        
    Example
    -------
    >>> files = get_all_tiff_files('input_dir')  # 26 positions × 10 wells = 260 files
    >>> downsampled = apply_position_downsampling(files, num_positions=10, seed=42)
    >>> len(downsampled)  # 10 positions × 10 wells = 100 files
    100
    """
    if num_positions <= 0:
        print("Warning: num_positions <= 0, returning empty list")
        return []
    
    # Group files
    grouped = group_files_by_well_and_position(file_list)
    
    if not grouped:
        print("Warning: Could not group any files by position")
        return file_list
    
    # Select positions
    selected_positions = select_random_positions(grouped, num_positions, seed)
    
    if not selected_positions:
        print("Warning: No positions selected")
        return []
    
    # Filter files
    filtered_files = []
    for well_id, positions in grouped.items():
        for pos_id in selected_positions:
            if pos_id in positions:
                filtered_files.extend(positions[pos_id])
    
    # Statistics
    if verbose:
        # Count wells and positions in result
        result_grouped = {}
        for f in filtered_files:
            well, pos = parse_filename_position(f)
            if well:
                if well not in result_grouped:
                    result_grouped[well] = set()
                if pos:
                    result_grouped[well].add(pos)
        
        print(f"\n{'='*60}")
        print(f"POSITION DOWNSAMPLING STATISTICS")
        print(f"{'='*60}")
        print(f"Original files:        {len(file_list)}")
        print(f"Downsampled files:     {len(filtered_files)}")
        print(f"Reduction:             {(1 - len(filtered_files)/len(file_list))*100:.1f}%")
        print(f"\nWells:")
        print(f"  Original wells:      {len(grouped)}")
        print(f"  Processed wells:     {len(result_grouped)} (ALL WELLS)")
        print(f"\nPositions:")
        print(f"  Requested per well:  {num_positions}")
        print(f"  Selected positions:  {', '.join(selected_positions[:10])}")
        if len(selected_positions) > 10:
            print(f"                       ... and {len(selected_positions)-10} more")
        print(f"  Avg per well:        {len(filtered_files)/len(result_grouped):.1f} files")
        print(f"\nRandom seed:           {seed} (for reproducibility)")
        print(f"{'='*60}\n")
    
    return filtered_files


def get_position_statistics(file_list):
    """
    Get statistics about positions in the dataset.
    
    Useful for deciding how many positions to select.
    
    Parameters
    ----------
    file_list : list of str/Path
        Input file paths
        
    Returns
    -------
    dict
        Statistics including:
        - num_wells: Number of unique wells
        - positions_per_well: Dict of well -> num positions
        - common_positions: Number of positions common to all wells
        - total_positions: Total unique positions across all wells
        
    Example
    -------
    >>> stats = get_position_statistics(file_list)
    >>> print(f"Wells: {stats['num_wells']}, Avg positions: {stats['avg_positions']:.1f}")
    Wells: 10, Avg positions: 26.3
    """
    grouped = group_files_by_well_and_position(file_list)
    
    if not grouped:
        return {
            'num_wells': 0,
            'positions_per_well': {},
            'common_positions': 0,
            'total_positions': 0,
            'avg_positions': 0,
            'min_positions': 0,
            'max_positions': 0
        }
    
    positions_per_well = {well: len(positions) for well, positions in grouped.items()}
    
    # Common positions
    well_positions = {well: set(positions.keys()) for well, positions in grouped.items()}
    common_positions = len(set.intersection(*well_positions.values())) if well_positions else 0
    
    # Total unique positions
    all_positions = set()
    for positions in well_positions.values():
        all_positions.update(positions)
    
    stats = {
        'num_wells': len(grouped),
        'positions_per_well': positions_per_well,
        'common_positions': common_positions,
        'total_positions': len(all_positions),
        'avg_positions': np.mean(list(positions_per_well.values())),
        'min_positions': min(positions_per_well.values()),
        'max_positions': max(positions_per_well.values()),
    }
    
    return stats


# Example usage and testing
if __name__ == "__main__":
    # Test filename parsing
    test_files = [
        'r01c01f01p01-ch1sk1fk1fl1.tiff',
        'well_B03_s05_z001.tif',
        'f001_s01_z001.tif',
        'A01_p03_t001.tif',
    ]
    
    print("Testing filename parsing:")
    print("="*60)
    for f in test_files:
        well, pos = parse_filename_position(f)
        print(f"{f:40s} -> well={well:6s}, pos={pos}")
    
    # Test downsampling
    print("\n\nTesting downsampling:")
    print("="*60)
    
    # Simulate file list
    simulated_files = []
    for well in ['A01', 'A02', 'A03']:
        for pos in range(1, 27):  # 26 positions per well
            for z in range(1, 6):  # 5 z-slices
                simulated_files.append(f'well_{well}_s{pos:02d}_z{z:03d}.tif')
    
    print(f"Created {len(simulated_files)} simulated files")
    
    # Get statistics
    stats = get_position_statistics(simulated_files)
    print(f"\nDataset statistics:")
    print(f"  Wells: {stats['num_wells']}")
    print(f"  Common positions: {stats['common_positions']}")
    print(f"  Avg positions per well: {stats['avg_positions']:.1f}")
    
    # Apply downsampling
    downsampled = apply_position_downsampling(
        simulated_files, 
        num_positions=10, 
        seed=42,
        verbose=True
    )
    
    print(f"\nDownsampled to {len(downsampled)} files")
    print(f"Reduction: {(1 - len(downsampled)/len(simulated_files))*100:.1f}%")