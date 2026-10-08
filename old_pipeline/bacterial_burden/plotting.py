#!/usr/bin/env python3
"""
Bacterial Burden Analysis Plotting Script
Generates violin+box plots for 1h and 3h treatments
Windows-compatible version
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import os

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# Get the directory where the script is located
script_dir = Path(__file__).parent
data_dir = script_dir  # Assuming data files are in the same directory

# Load data - update these paths if your CSV files are elsewhere
data_1h_path = data_dir / 'Results_SingleCell_BacterialBurden_1htreatment_20260210_171315.csv'
data_3h_path = data_dir / 'Results_SingleCell_BacterialBurden_3htreatment_20260210_211203.csv'

# Check if files exist
if not data_1h_path.exists():
    print(f"Error: Cannot find {data_1h_path}")
    print(f"Current directory: {os.getcwd()}")
    exit(1)
if not data_3h_path.exists():
    print(f"Error: Cannot find {data_3h_path}")
    print(f"Current directory: {os.getcwd()}")
    exit(1)

# Load data
data_1h = pd.read_csv(data_1h_path)
data_3h = pd.read_csv(data_3h_path)

# Add timepoint column
data_1h['timepoint'] = '1h'
data_3h['timepoint'] = '3h'

# Get unique strains (excluding header row if present)
strains = sorted([s for s in data_1h['treatment'].unique() if s != 'treatment'])
print(f"Strains found: {strains}")

# Define color palettes
colors_1h = sns.color_palette("Set2", n_colors=len(strains))
colors_3h = sns.color_palette("Set1", n_colors=len(strains))

# Create output directory in the same folder as the script
output_dir = script_dir / 'bacterial_burden_plots'
output_dir.mkdir(exist_ok=True)
print(f"Output directory: {output_dir}")

# ============================================================================
# PART 1: SEPARATE VIOLIN+BOX PLOTS FOR 1h AND 3h
# ============================================================================

metrics = [
    ('area_ratio', 'Bacterial Burden (Area Ratio)', 'Area Ratio'),
    ('intracellular_bact_area', 'Intracellular Bacterial Area (pixels)', 'Bacterial Area (pixels)'),
    ('host_area', 'Host Cell Area (pixels)', 'Cell Area (pixels)')
]

for metric, title_suffix, ylabel in metrics:
    # 1h plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Create violin plot
    parts = ax.violinplot(
        [data_1h[data_1h['treatment'] == strain][metric].dropna().values 
         for strain in strains],
        positions=range(len(strains)),
        widths=0.7,
        showmeans=False,
        showmedians=False,
        showextrema=False
    )
    
    # Color the violins
    for pc, color in zip(parts['bodies'], colors_1h):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
    
    # Add boxplots on top
    bp = ax.boxplot(
        [data_1h[data_1h['treatment'] == strain][metric].dropna().values 
         for strain in strains],
        positions=range(len(strains)),
        widths=0.3,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color='red', linewidth=2),
        boxprops=dict(facecolor='white', alpha=0.8, linewidth=1.5),
        whiskerprops=dict(linewidth=1.5),
        capprops=dict(linewidth=1.5)
    )
    
    ax.set_xlabel('Strain', fontsize=12, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(f'{title_suffix} - 1h Treatment', fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(strains)))
    ax.set_xticklabels(strains, rotation=45, ha='right')
    ax.grid(True, alpha=0.3)
    
    # Add mean values above each violin
    y_max = ax.get_ylim()[1]
    for i, strain in enumerate(strains):
        mean_val = data_1h[data_1h['treatment'] == strain][metric].mean()
        ax.text(i, y_max * 0.98, f'μ={mean_val:.4f}', 
                ha='center', va='top', fontsize=8, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_dir / f'{metric}_1h_violin_box.png', bbox_inches='tight')
    plt.close()
    print(f"Saved: {metric}_1h_violin_box.png")
    
    # 3h plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Create violin plot
    parts = ax.violinplot(
        [data_3h[data_3h['treatment'] == strain][metric].dropna().values 
         for strain in strains],
        positions=range(len(strains)),
        widths=0.7,
        showmeans=False,
        showmedians=False,
        showextrema=False
    )
    
    # Color the violins
    for pc, color in zip(parts['bodies'], colors_3h):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
    
    # Add boxplots on top
    bp = ax.boxplot(
        [data_3h[data_3h['treatment'] == strain][metric].dropna().values 
         for strain in strains],
        positions=range(len(strains)),
        widths=0.3,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color='red', linewidth=2),
        boxprops=dict(facecolor='white', alpha=0.8, linewidth=1.5),
        whiskerprops=dict(linewidth=1.5),
        capprops=dict(linewidth=1.5)
    )
    
    ax.set_xlabel('Strain', fontsize=12, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(f'{title_suffix} - 3h Treatment', fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(strains)))
    ax.set_xticklabels(strains, rotation=45, ha='right')
    ax.grid(True, alpha=0.3)
    
    # Add mean values above each violin
    y_max = ax.get_ylim()[1]
    for i, strain in enumerate(strains):
        mean_val = data_3h[data_3h['treatment'] == strain][metric].mean()
        ax.text(i, y_max * 0.98, f'μ={mean_val:.4f}', 
                ha='center', va='top', fontsize=8, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_dir / f'{metric}_3h_violin_box.png', bbox_inches='tight')
    plt.close()
    print(f"Saved: {metric}_3h_violin_box.png")

# ============================================================================
# PART 2: SIDE-BY-SIDE COMPARISON (1h vs 3h for all strains)
# ============================================================================

# Prepare combined data
data_combined = pd.concat([data_1h, data_3h], ignore_index=True)

for metric, title_suffix, ylabel in metrics:
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Prepare positions for side-by-side violins
    x_positions = []
    labels = []
    for i, strain in enumerate(strains):
        x_positions.extend([i*3, i*3 + 1])
        labels.extend([f'{strain}\n1h', f'{strain}\n3h'])
    
    # Create violin plots
    data_for_violins = []
    colors_for_violins = []
    
    for i, strain in enumerate(strains):
        data_for_violins.append(data_1h[data_1h['treatment'] == strain][metric].dropna().values)
        data_for_violins.append(data_3h[data_3h['treatment'] == strain][metric].dropna().values)
        colors_for_violins.extend([colors_1h[i], colors_3h[i]])
    
    parts = ax.violinplot(
        data_for_violins,
        positions=x_positions,
        widths=0.8,
        showmeans=False,
        showmedians=False,
        showextrema=False
    )
    
    # Color the violins
    for pc, color in zip(parts['bodies'], colors_for_violins):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
    
    # Add boxplots
    bp = ax.boxplot(
        data_for_violins,
        positions=x_positions,
        widths=0.4,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color='red', linewidth=2),
        boxprops=dict(facecolor='white', alpha=0.8, linewidth=1.5),
        whiskerprops=dict(linewidth=1.5),
        capprops=dict(linewidth=1.5)
    )
    
    ax.set_xlabel('Strain and Timepoint', fontsize=12, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(f'{title_suffix} - 1h vs 3h Comparison', fontsize=14, fontweight='bold')
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    ax.grid(True, alpha=0.3)
    
    # Add mean values above each violin
    y_max = ax.get_ylim()[1]
    for i, strain in enumerate(strains):
        mean_1h = data_1h[data_1h['treatment'] == strain][metric].mean()
        mean_3h = data_3h[data_3h['treatment'] == strain][metric].mean()
        ax.text(i*3, y_max * 0.98, f'μ={mean_1h:.3f}', 
                ha='center', va='top', fontsize=7, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
        ax.text(i*3 + 1, y_max * 0.98, f'μ={mean_3h:.3f}', 
                ha='center', va='top', fontsize=7, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=colors_1h[0], alpha=0.7, edgecolor='black', label='1h'),
        Patch(facecolor=colors_3h[0], alpha=0.7, edgecolor='black', label='3h')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_dir / f'{metric}_comparison_1h_vs_3h.png', bbox_inches='tight')
    plt.close()
    print(f"Saved: {metric}_comparison_1h_vs_3h.png")

# ============================================================================
# PART 3: INDIVIDUAL STRAIN COMPARISONS (Each strain separately)
# ============================================================================

for strain in strains:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f'Strain: {strain} - 1h vs 3h Comparison', fontsize=16, fontweight='bold')
    
    for ax, (metric, title_suffix, ylabel) in zip(axes, metrics):
        # Get data for this strain
        data_1h_strain = data_1h[data_1h['treatment'] == strain][metric].dropna().values
        data_3h_strain = data_3h[data_3h['treatment'] == strain][metric].dropna().values
        
        # Create violin plots
        parts = ax.violinplot(
            [data_1h_strain, data_3h_strain],
            positions=[0, 1],
            widths=0.7,
            showmeans=False,
            showmedians=False,
            showextrema=False
        )
        
        # Color the violins
        strain_idx = strains.index(strain)
        colors = [colors_1h[strain_idx], colors_3h[strain_idx]]
        for pc, color in zip(parts['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
            pc.set_edgecolor('black')
            pc.set_linewidth(1.5)
        
        # Add boxplots
        bp = ax.boxplot(
            [data_1h_strain, data_3h_strain],
            positions=[0, 1],
            widths=0.35,
            patch_artist=True,
            showfliers=False,
            medianprops=dict(color='red', linewidth=2.5),
            boxprops=dict(facecolor='white', alpha=0.9, linewidth=1.5),
            whiskerprops=dict(linewidth=1.5),
            capprops=dict(linewidth=1.5)
        )
        
        ax.set_ylabel(ylabel, fontsize=11, fontweight='bold')
        ax.set_title(title_suffix, fontsize=11, fontweight='bold')
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['1h', '3h'], fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Add sample sizes and mean values
        n_1h = len(data_1h_strain)
        n_3h = len(data_3h_strain)
        mean_1h = data_1h_strain.mean()
        mean_3h = data_3h_strain.mean()
        
        y_max = ax.get_ylim()[1]
        ax.text(0, y_max * 0.95, f'n={n_1h}\nμ={mean_1h:.3f}', 
                ha='center', va='top', fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        ax.text(1, y_max * 0.95, f'n={n_3h}\nμ={mean_3h:.3f}', 
                ha='center', va='top', fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    # Create safe filename
    safe_strain = strain.replace('+', 'plus').replace('/', '_')
    plt.savefig(output_dir / f'individual_{safe_strain}_1h_vs_3h.png', bbox_inches='tight')
    plt.close()
    print(f"Saved: individual_{safe_strain}_1h_vs_3h.png")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)

for strain in strains:
    print(f"\n{strain}:")
    print("-" * 40)
    
    for metric, _, _ in metrics:
        data_1h_strain = data_1h[data_1h['treatment'] == strain][metric].dropna()
        data_3h_strain = data_3h[data_3h['treatment'] == strain][metric].dropna()
        
        print(f"\n  {metric}:")
        print(f"    1h: n={len(data_1h_strain)}, mean={data_1h_strain.mean():.4f}, "
              f"median={data_1h_strain.median():.4f}, std={data_1h_strain.std():.4f}")
        print(f"    3h: n={len(data_3h_strain)}, mean={data_3h_strain.mean():.4f}, "
              f"median={data_3h_strain.median():.4f}, std={data_3h_strain.std():.4f}")

print("\n" + "="*80)
print(f"All plots saved to: {output_dir.absolute()}")
print("="*80)