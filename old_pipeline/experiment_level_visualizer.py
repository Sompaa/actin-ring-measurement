"""
Experiment-Level Visualization Module
Designed for large datasets (1000+ cells per experiment)
Generates aggregated, statistical plots instead of per-cell visualizations
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')


class ExperimentLevelVisualizer:
    """
    Create experiment-level aggregated visualizations for large datasets
    """
    
    def __init__(self, output_dir, config):
        self.output_dir = Path(output_dir)
        self.config = config
        
        # Create output directories
        (self.output_dir / 'experiment_summaries').mkdir(parents=True, exist_ok=True)
        (self.output_dir / 'group_comparisons').mkdir(parents=True, exist_ok=True)
        (self.output_dir / 'statistical_validation').mkdir(parents=True, exist_ok=True)
    
    def generate_experiment_summary(self, df, experiment_name, timestamp):
        """
        Generate comprehensive experiment-level summary
        
        Creates 3x3 grid with:
        - Sample size per group
        - Distribution plots (violin + box)
        - Effect size heatmap
        - Statistical test results
        
        Optimized for large N (>1000 cells)
        """
        
        fig = plt.figure(figsize=(20, 15))
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        fig.suptitle(f'Experiment Summary: {experiment_name}\n'
                     f'Total Cells: {len(df)} | Images: {df["image_name"].nunique()}',
                     fontsize=16, fontweight='bold')
        
        # Get group colors
        group_colors = self.config.get('group_colors', {})
        if not group_colors:
            unique_groups = df['group'].unique()
            palette = sns.color_palette('tab10', len(unique_groups))
            group_colors = dict(zip(unique_groups, palette))
        
        # Panel 1: Sample sizes
        ax = fig.add_subplot(gs[0, 0])
        group_counts = df['group'].value_counts().sort_index()
        bars = ax.barh(range(len(group_counts)), group_counts.values,
                       color=[group_colors.get(g, 'gray') for g in group_counts.index])
        ax.set_yticks(range(len(group_counts)))
        ax.set_yticklabels(group_counts.index)
        ax.set_xlabel('Number of Cells', fontweight='bold')
        ax.set_title('Sample Size per Group', fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        
        # Add count labels
        for i, (bar, count) in enumerate(zip(bars, group_counts.values)):
            ax.text(count + max(group_counts)*0.01, i, f'n={count}',
                    va='center', fontweight='bold')
        
        # Panel 2: Peripheral Gradient (Violin + Strip)
        ax = fig.add_subplot(gs[0, 1])
        self._plot_violin_strip(ax, df, 'peripheral_gradient', 'Peripheral Gradient\n(Cortical Fortification)',
                                group_colors)
        ax.axhline(y=0.001, color='red', linestyle='--', linewidth=2, alpha=0.5, label='Threshold')
        ax.legend()
        
        # Panel 3: Spectral Entropy (Violin + Strip)
        ax = fig.add_subplot(gs[0, 2])
        self._plot_violin_strip(ax, df, 'cortical_spectral_entropy', 'Spectral Entropy\n(Network Periodicity)',
                                group_colors)
        ax.axhline(y=2.0, color='orange', linestyle='--', linewidth=2, alpha=0.5, label='Periodic/Chaotic')
        ax.legend()
        
        # Panel 4: Cortical/Inner Ratio (Violin + Strip)
        ax = fig.add_subplot(gs[1, 0])
        self._plot_violin_strip(ax, df, 'cortical_inner_ratio', 'Cortical/Inner Ratio',
                                group_colors)
        ax.axhline(y=1.0, color='black', linestyle='--', linewidth=2, alpha=0.5, label='Ratio = 1.0')
        ax.legend()
        
        # Panel 5: Skeleton Density (Violin + Strip)
        ax = fig.add_subplot(gs[1, 1])
        self._plot_violin_strip(ax, df, 'skeleton_density', 'Skeleton Density',
                                group_colors)
        
        # Panel 6: Hole Size Edge (Violin + Strip)
        ax = fig.add_subplot(gs[1, 2])
        self._plot_violin_strip(ax, df, 'hole_size_edge', 'Hole Size (Edge)\n(Mesh Tightness)',
                                group_colors)
        
        # Panel 7: Effect Size Heatmap (Cohen's d vs Control)
        ax = fig.add_subplot(gs[2, 0])
        self._plot_effect_size_heatmap(ax, df, group_colors)
        
        # Panel 8: Statistical Significance Matrix
        ax = fig.add_subplot(gs[2, 1])
        self._plot_significance_matrix(ax, df, 'peripheral_gradient')
        
        # Panel 9: Group Means with Error Bars
        ax = fig.add_subplot(gs[2, 2])
        self._plot_group_means_comparison(ax, df, 'peripheral_gradient', 
                                          'Mean Peripheral Gradient\n(± SEM)', group_colors)
        
        # Save
        plt.tight_layout()
        save_path = self.output_dir / 'experiment_summaries' / f'{experiment_name}_summary_{timestamp}.png'
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        return save_path
    
    def _plot_violin_strip(self, ax, df, metric, title, group_colors):
        """Plot violin + strip for large N"""
        
        # Remove NaN
        plot_df = df[df[metric].notna()].copy()
        
        if len(plot_df) == 0:
            ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title, fontweight='bold')
            return
        
        # Violin plot (shows distribution)
        parts = ax.violinplot([plot_df[plot_df['group'] == g][metric].values 
                               for g in plot_df['group'].unique()],
                              positions=range(len(plot_df['group'].unique())),
                              widths=0.7, showmeans=False, showmedians=False)
        
        # Color violins
        for i, (pc, group) in enumerate(zip(parts['bodies'], plot_df['group'].unique())):
            pc.set_facecolor(group_colors.get(group, 'gray'))
            pc.set_alpha(0.6)
        
        # Box plot overlay (shows quartiles)
        ax.boxplot([plot_df[plot_df['group'] == g][metric].values 
                    for g in plot_df['group'].unique()],
                   positions=range(len(plot_df['group'].unique())),
                   widths=0.2, showfliers=False,
                   boxprops=dict(color='black', linewidth=1.5),
                   whiskerprops=dict(color='black', linewidth=1.5),
                   capprops=dict(color='black', linewidth=1.5),
                   medianprops=dict(color='red', linewidth=2))
        
        ax.set_xticks(range(len(plot_df['group'].unique())))
        ax.set_xticklabels(plot_df['group'].unique(), rotation=45, ha='right')
        ax.set_ylabel(metric.replace('_', ' ').title(), fontweight='bold')
        ax.set_title(title, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
    
    def _plot_effect_size_heatmap(self, ax, df, group_colors):
        """Cohen's d effect sizes vs control"""
        
        metrics = ['peripheral_gradient', 'cortical_spectral_entropy', 
                   'cortical_inner_ratio', 'skeleton_density', 'hole_size_edge']
        
        groups = [g for g in df['group'].unique() if g.lower() != 'control' and g.lower() != 'empty']
        
        # Find control group
        control_group = None
        for g in df['group'].unique():
            if 'control' in g.lower():
                control_group = g
                break
        
        if control_group is None:
            ax.text(0.5, 0.5, 'No control group found', ha='center', va='center', transform=ax.transAxes)
            return
        
        # Calculate Cohen's d
        effect_sizes = []
        for metric in metrics:
            row = []
            control_data = df[df['group'] == control_group][metric].dropna()
            
            for group in groups:
                treatment_data = df[df['group'] == group][metric].dropna()
                
                if len(control_data) > 0 and len(treatment_data) > 0:
                    # Cohen's d
                    pooled_std = np.sqrt((control_data.std()**2 + treatment_data.std()**2) / 2)
                    if pooled_std > 0:
                        d = (treatment_data.mean() - control_data.mean()) / pooled_std
                        row.append(d)
                    else:
                        row.append(0)
                else:
                    row.append(np.nan)
            
            effect_sizes.append(row)
        
        # Plot heatmap
        im = ax.imshow(effect_sizes, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
        
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=45, ha='right')
        ax.set_yticks(range(len(metrics)))
        ax.set_yticklabels([m.replace('_', ' ').title() for m in metrics])
        ax.set_title("Effect Size (Cohen's d) vs Control", fontweight='bold')
        
        # Add text annotations
        for i in range(len(metrics)):
            for j in range(len(groups)):
                if not np.isnan(effect_sizes[i][j]):
                    text = ax.text(j, i, f'{effect_sizes[i][j]:.2f}',
                                   ha='center', va='center', color='white' if abs(effect_sizes[i][j]) > 1 else 'black',
                                   fontweight='bold')
        
        plt.colorbar(im, ax=ax, label="Cohen's d")
    
    def _plot_significance_matrix(self, ax, df, metric):
        """P-value matrix for pairwise comparisons"""
        
        groups = df['group'].unique()
        n_groups = len(groups)
        
        p_matrix = np.ones((n_groups, n_groups))
        
        for i, g1 in enumerate(groups):
            for j, g2 in enumerate(groups):
                if i < j:
                    data1 = df[df['group'] == g1][metric].dropna()
                    data2 = df[df['group'] == g2][metric].dropna()
                    
                    if len(data1) > 0 and len(data2) > 0:
                        _, p = stats.mannwhitneyu(data1, data2, alternative='two-sided')
                        p_matrix[i, j] = p
                        p_matrix[j, i] = p
        
        # Plot
        im = ax.imshow(p_matrix, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=0.05)
        
        ax.set_xticks(range(n_groups))
        ax.set_xticklabels(groups, rotation=45, ha='right')
        ax.set_yticks(range(n_groups))
        ax.set_yticklabels(groups)
        ax.set_title(f'Pairwise P-values\n({metric.replace("_", " ").title()})', fontweight='bold')
        
        # Add significance markers
        for i in range(n_groups):
            for j in range(n_groups):
                if i != j:
                    p = p_matrix[i, j]
                    if p < 0.001:
                        text = '***'
                    elif p < 0.01:
                        text = '**'
                    elif p < 0.05:
                        text = '*'
                    else:
                        text = 'ns'
                    
                    ax.text(j, i, text, ha='center', va='center',
                            color='white' if p < 0.01 else 'black', fontweight='bold')
        
        plt.colorbar(im, ax=ax, label='P-value')
    
    def _plot_group_means_comparison(self, ax, df, metric, title, group_colors):
        """Bar plot with means and SEM error bars"""
        
        groups = df['group'].unique()
        means = []
        sems = []
        
        for group in groups:
            data = df[df['group'] == group][metric].dropna()
            if len(data) > 0:
                means.append(data.mean())
                sems.append(data.sem())
            else:
                means.append(0)
                sems.append(0)
        
        x = range(len(groups))
        bars = ax.bar(x, means, yerr=sems, capsize=5, 
                      color=[group_colors.get(g, 'gray') for g in groups],
                      edgecolor='black', linewidth=1.5, alpha=0.7)
        
        ax.set_xticks(x)
        ax.set_xticklabels(groups, rotation=45, ha='right')
        ax.set_ylabel(metric.replace('_', ' ').title(), fontweight='bold')
        ax.set_title(title, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        # Add n values
        for i, (bar, group) in enumerate(zip(bars, groups)):
            n = len(df[df['group'] == group][metric].dropna())
            ax.text(i, means[i] + sems[i] + max(means)*0.05, f'n={n}',
                    ha='center', fontsize=9)
    
    def generate_multi_experiment_comparison(self, all_dfs, experiment_names, timestamp):
        """
        Compare multiple experiments side-by-side
        
        Parameters:
        -----------
        all_dfs : list of DataFrames
            One DataFrame per experiment
        experiment_names : list of str
            Names of experiments
        timestamp : str
            Timestamp string
        """
        
        # Combine all data with experiment labels
        combined = []
        for df, exp_name in zip(all_dfs, experiment_names):
            df_copy = df.copy()
            df_copy['experiment'] = exp_name
            combined.append(df_copy)
        
        full_df = pd.concat(combined, ignore_index=True)
        
        # Create comparison plot
        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        fig.suptitle(f'Multi-Experiment Comparison ({len(experiment_names)} experiments)',
                     fontsize=16, fontweight='bold')
        
        metrics = [
            ('peripheral_gradient', 'Peripheral Gradient'),
            ('cortical_spectral_entropy', 'Spectral Entropy'),
            ('cortical_inner_ratio', 'Cortical/Inner Ratio'),
            ('skeleton_density', 'Skeleton Density'),
            ('hole_size_edge', 'Hole Size (Edge)'),
            ('mean_tortuosity_edge', 'Path Tortuosity (Edge)')
        ]
        
        for ax, (metric, label) in zip(axes.flat, metrics):
            # Group by experiment and group
            plot_data = full_df[full_df[metric].notna()].copy()
            
            if len(plot_data) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(label)
                continue
            
            # Create grouped box plot
            unique_experiments = plot_data['experiment'].unique()
            unique_groups = plot_data['group'].unique()
            
            x_pos = []
            data_to_plot = []
            labels = []
            colors = []
            
            pos = 0
            for exp in unique_experiments:
                for group in unique_groups:
                    subset = plot_data[(plot_data['experiment'] == exp) & (plot_data['group'] == group)]
                    if len(subset) > 0:
                        x_pos.append(pos)
                        data_to_plot.append(subset[metric].values)
                        labels.append(f'{exp}\n{group}')
                        colors.append(self.config.get('group_colors', {}).get(group, 'gray'))
                        pos += 1
                pos += 1  # Gap between experiments
            
            bp = ax.boxplot(data_to_plot, positions=x_pos, widths=0.6,
                           patch_artist=True, showfliers=False)
            
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            
            ax.set_xticks(x_pos)
            ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
            ax.set_ylabel(label, fontweight='bold')
            ax.set_title(label, fontweight='bold')
            ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        save_path = self.output_dir / 'group_comparisons' / f'multi_experiment_comparison_{timestamp}.png'
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        return save_path