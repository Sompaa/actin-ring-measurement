"""
Visualization Module - Results Visualization
Creates side-by-side comparison images and summary visualizations
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend


class ResultsVisualizer:
    """Create visualizations for results"""
    
    def __init__(self, config):
        self.config = config
        self.output_dir = Path(config['output_dir']) / 'comparisons'
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def create_comparison(self, cell_mask, actin_image, skeleton, image_name, cell_id):
        """
        Create side-by-side comparison of original segmented cell and skeletonized cell
        
        Parameters:
        -----------
        cell_mask : ndarray
            Binary mask of the cell
        actin_image : ndarray
            Full actin image
        skeleton : ndarray
            Skeleton of the cell (cropped)
        image_name : str
            Name of source image
        cell_id : int
            Cell identifier
        """
        try:
            # Extract cell region from full image
            rows, cols = np.where(cell_mask)
            if len(rows) == 0:
                return
            
            r0, r1 = rows.min(), rows.max() + 1
            c0, c1 = cols.min(), cols.max() + 1
            
            # Add padding
            pad = 10
            r0 = max(0, r0 - pad)
            c0 = max(0, c0 - pad)
            r1 = min(actin_image.shape[0], r1 + pad)
            c1 = min(actin_image.shape[1], c1 + pad)
            
            # Crop original cell
            cell_region = actin_image[r0:r1, c0:c1]
            cell_mask_crop = cell_mask[r0:r1, c0:c1]
            
            # Mask background
            cell_region_masked = cell_region * cell_mask_crop
            
            # Resize skeleton to match if needed
            if skeleton.shape != cell_region.shape:
                # Skeleton might be smaller due to cropping in skeletonization
                # Center it in a canvas of the right size
                canvas = np.zeros_like(cell_region_masked, dtype=bool)
                h_offset = (canvas.shape[0] - skeleton.shape[0]) // 2
                w_offset = (canvas.shape[1] - skeleton.shape[1]) // 2
                canvas[
                    h_offset:h_offset + skeleton.shape[0],
                    w_offset:w_offset + skeleton.shape[1]
                ] = skeleton
                skeleton_display = canvas
            else:
                skeleton_display = skeleton
            
            # Create figure
            fig, axes = plt.subplots(1, 2, figsize=(10, 5))
            
            # Left: Original segmented cell
            axes[0].imshow(cell_region_masked, cmap='gray')
            axes[0].contour(cell_mask_crop, colors='cyan', linewidths=1, levels=[0.5])
            axes[0].set_title('Original Segmented Cell', fontsize=12, fontweight='bold')
            axes[0].axis('off')
            
            # Right: Skeletonized cell
            axes[1].imshow(cell_region_masked, cmap='gray', alpha=0.7)
            axes[1].imshow(skeleton_display, cmap='hot', alpha=0.8)
            axes[1].contour(cell_mask_crop, colors='cyan', linewidths=1, levels=[0.5])
            axes[1].set_title('Skeletonized Cytoskeleton', fontsize=12, fontweight='bold')
            axes[1].axis('off')
            
            # Add overall title
            fig.suptitle(
                f'{image_name} - Cell {cell_id:03d}',
                fontsize=14,
                fontweight='bold',
                y=0.98
            )
            
            plt.tight_layout()
            
            # Save
            output_file = self.output_dir / f'{image_name}_cell{cell_id:03d}_comparison.png'
            plt.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
            plt.close(fig)
            
        except Exception as e:
            print(f'  ⚠️  Error creating comparison for cell {cell_id}: {str(e)}')
    
    def create_batch_summary(self, df, output_path):
        """
        Create summary visualizations for a batch of results
        
        Parameters:
        -----------
        df : pandas.DataFrame
            Results dataframe
        output_path : Path
            Path to save summary figure
        """
        import seaborn as sns
        
        fig, axes = plt.subplots(3, 3, figsize=(18, 15))
        fig.suptitle('Batch Analysis Summary', fontsize=16, fontweight='bold')
        
        # 1. Skeleton length distribution
        axes[0, 0].hist(df['total_length'], bins=50, edgecolor='black', alpha=0.7)
        axes[0, 0].set_xlabel('Total Skeleton Length (pixels)', fontsize=10)
        axes[0, 0].set_ylabel('Count', fontsize=10)
        axes[0, 0].set_title('Skeleton Length Distribution')
        axes[0, 0].grid(alpha=0.3)
        
        # 2. Number of branches
        axes[0, 1].hist(df['num_branches'], bins=30, edgecolor='black', alpha=0.7, color='orange')
        axes[0, 1].set_xlabel('Number of Branches', fontsize=10)
        axes[0, 1].set_ylabel('Count', fontsize=10)
        axes[0, 1].set_title('Branch Count Distribution')
        axes[0, 1].grid(alpha=0.3)
        
        # 3. Mean branch length
        axes[0, 2].hist(df['mean_branch_length'], bins=50, edgecolor='black', alpha=0.7, color='green')
        axes[0, 2].set_xlabel('Mean Branch Length (pixels)', fontsize=10)
        axes[0, 2].set_ylabel('Count', fontsize=10)
        axes[0, 2].set_title('Mean Branch Length Distribution')
        axes[0, 2].grid(alpha=0.3)
        
        # 4. Cell area
        axes[1, 0].hist(df['cell_area'], bins=50, edgecolor='black', alpha=0.7, color='purple')
        axes[1, 0].set_xlabel('Cell Area (pixels²)', fontsize=10)
        axes[1, 0].set_ylabel('Count', fontsize=10)
        axes[1, 0].set_title('Cell Area Distribution')
        axes[1, 0].grid(alpha=0.3)
        
        # 5. Skeleton density
        axes[1, 1].hist(df['skeleton_density'], bins=50, edgecolor='black', alpha=0.7, color='red')
        axes[1, 1].set_xlabel('Skeleton Density', fontsize=10)
        axes[1, 1].set_ylabel('Count', fontsize=10)
        axes[1, 1].set_title('Skeleton Density Distribution')
        axes[1, 1].grid(alpha=0.3)
        
        # 6. Cell solidity
        axes[1, 2].hist(df['cell_solidity'], bins=50, edgecolor='black', alpha=0.7, color='teal')
        axes[1, 2].set_xlabel('Cell Solidity', fontsize=10)
        axes[1, 2].set_ylabel('Count', fontsize=10)
        axes[1, 2].set_title('Cell Solidity Distribution')
        axes[1, 2].grid(alpha=0.3)
        
        # 7. Skeleton length vs cell area
        axes[2, 0].scatter(df['cell_area'], df['total_length'], alpha=0.5, s=20)
        axes[2, 0].set_xlabel('Cell Area (pixels²)', fontsize=10)
        axes[2, 0].set_ylabel('Total Skeleton Length (pixels)', fontsize=10)
        axes[2, 0].set_title('Cell Area vs Skeleton Length')
        axes[2, 0].grid(alpha=0.3)
        
        # 8. Branches vs skeleton length
        axes[2, 1].scatter(df['num_branches'], df['total_length'], alpha=0.5, s=20, color='orange')
        axes[2, 1].set_xlabel('Number of Branches', fontsize=10)
        axes[2, 1].set_ylabel('Total Skeleton Length (pixels)', fontsize=10)
        axes[2, 1].set_title('Branches vs Skeleton Length')
        axes[2, 1].grid(alpha=0.3)
        
        # 9. Cells per image
        if 'image_name' in df.columns:
            cells_per_image = df.groupby('image_name').size()
            axes[2, 2].bar(range(len(cells_per_image)), cells_per_image.values)
            axes[2, 2].set_xlabel('Image Index', fontsize=10)
            axes[2, 2].set_ylabel('Number of Cells', fontsize=10)
            axes[2, 2].set_title('Cells Per Image')
            axes[2, 2].grid(alpha=0.3)
        else:
            axes[2, 2].text(0.5, 0.5, 'No image grouping', 
                           ha='center', va='center', fontsize=12)
            axes[2, 2].axis('off')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        
        print(f'✓ Summary visualization saved to: {output_path}')
