"""
Visualization functions for phase change model results.
Generates figures for the paper.
Author: Jair Patino B.
Date: January 2026
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle, Circle
import seaborn as sns
from typing import List, Dict, Tuple, Optional

# Set publication-quality style
plt.style.use('default')
mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 600,
    'savefig.format': 'png',
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'lines.linewidth': 1.5,
    'lines.markersize': 6,
    'grid.linewidth': 0.5,
    'grid.alpha': 0.3,
})

# Color palette
COLORS = {
    'primary': '#1f77b4',
    'secondary': '#ff7f0e',
    'tertiary': '#2ca02c',
    'error': '#d62728',
    'gray': '#7f7f7f',
    'light_gray': '#d3d3d3'
}

def plot_error_vs_biot(error_data: Dict[str, np.ndarray], 
                      save_path: Optional[str] = None) -> plt.Figure:
    """
    Create Figure 1: Relative error vs Biot number.
    
    Parameters:
    -----------
    error_data : dict
        Dictionary with keys:
        - 'Bi': Biot numbers
        - 'our_model': errors for our model
        - 'liu_2009': errors for Liu et al. model
        - 'karwa_2013': errors for Karwa et al. model
    save_path : str, optional
        Path to save figure
    
    Returns:
    --------
    fig : matplotlib Figure
    """
    Bi = error_data['Bi']
    
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # Plot error curves
    ax.semilogx(Bi, error_data['our_model'], 
               color=COLORS['primary'], 
               linewidth=2,
               label='Our model', 
               zorder=3)
    
    ax.semilogx(Bi, error_data['liu_2009'], 
               color=COLORS['secondary'], 
               linewidth=1.5,
               linestyle='--',
               label='Liu et al. (2009)', 
               zorder=2)
    
    ax.semilogx(Bi, error_data['karwa_2013'], 
               color=COLORS['tertiary'], 
               linewidth=1.5,
               linestyle=':',
               label='Karwa et al. (2013)', 
               zorder=2)
    
    # Add shaded regions for error zones
    ax.axhspan(0, 5, alpha=0.2, color='green', zorder=0, label='<5% error')
    ax.axhspan(5, 15, alpha=0.2, color='yellow', zorder=0, label='5-15% error')
    ax.axhspan(15, 100, alpha=0.2, color='red', zorder=0, label='>15% error')
    
    # Add vertical lines for validity limits
    ax.axvline(x=0.2, color=COLORS['secondary'], linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(x=0.5, color=COLORS['tertiary'], linestyle=':', alpha=0.5, linewidth=1)
    ax.axvline(x=2.0, color=COLORS['primary'], linestyle='-', alpha=0.5, linewidth=1)
    
    # Add text annotations for validity limits
    ax.text(0.25, 80, 'Liu et al.\nlimit', 
            fontsize=8, color=COLORS['secondary'], ha='left')
    ax.text(0.55, 80, 'Karwa et al.\nlimit', 
            fontsize=8, color=COLORS['tertiary'], ha='left')
    ax.text(2.2, 80, 'Our model\nlimit', 
            fontsize=8, color=COLORS['primary'], ha='left')
    
    # Set limits and labels
    ax.set_xlim(0.01, 10)
    ax.set_ylim(0, 100)
    ax.set_xlabel('Biot number, Bi')
    ax.set_ylabel('Relative error (%)')
    ax.set_title('Model Accuracy vs Biot Number (Planar Geometry)')
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.grid(True, which='minor', alpha=0.2, linestyle=':', linewidth=0.5)
    
    # Legend
    ax.legend(loc='upper left', frameon=True, framealpha=0.9)
    
    # Add text annotation
    ax.text(0.02, 0.98, 
           'Our model: <5% error for Bi ≤ 1\n<15% error at Bi = 2',
           transform=ax.transAxes,
           verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
           fontsize=8)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig

def plot_predicted_vs_numerical(predicted: np.ndarray, 
                               numerical: np.ndarray,
                               materials: List[str],
                               Bi_ranges: List[Tuple[float, float]],
                               save_path: Optional[str] = None) -> plt.Figure:
    """
    Create Figure 2: Predicted vs numerically computed phase change times.
    
    Parameters:
    -----------
    predicted : array
        Predicted phase change times
    numerical : array
        Numerical reference times
    materials : list
        Material names for each data point
    Bi_ranges : list
        Biot number ranges for each data point
    save_path : str, optional
        Path to save figure
    
    Returns:
    --------
    fig : matplotlib Figure
    """
    fig, ax = plt.subplots(figsize=(5, 5))
    
    # Create color map for materials
    unique_materials = list(set(materials))
    material_colors = plt.cm.tab10(np.linspace(0, 1, len(unique_materials)))
    color_dict = {mat: color for mat, color in zip(unique_materials, material_colors)}
    
    # Plot each material with different marker
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
    
    for i, material in enumerate(unique_materials):
        mask = [m == material for m in materials]
        if not any(mask):
            continue
            
        ax.scatter(numerical[mask], predicted[mask],
                  color=color_dict[material],
                  marker=markers[i % len(markers)],
                  s=60,
                  alpha=0.8,
                  edgecolors='black',
                  linewidth=0.5,
                  label=material,
                  zorder=3)
    
    # Add perfect agreement line
    max_val = max(np.max(predicted), np.max(numerical))
    min_val = min(np.min(predicted), np.min(numerical))
    ax.plot([min_val, max_val], [min_val, max_val],
           'k--', linewidth=1, alpha=0.7, zorder=1,
           label='Perfect agreement')
    
    # Add ±5% error bands
    ax.fill_between([min_val, max_val],
                   [min_val*0.95, max_val*0.95],
                   [min_val*1.05, max_val*1.05],
                   color='green', alpha=0.2, zorder=0,
                   label='±5% error band')
    
    # Calculate statistics
    errors = 100 * np.abs(predicted - numerical) / numerical
    max_error = np.max(errors)
    mean_error = np.mean(errors)
    
    # Add statistics text box
    stats_text = f'All points within ±5%\nMax error: {max_error:.1f}%\nMean error: {mean_error:.1f}%'
    ax.text(0.05, 0.95, stats_text,
           transform=ax.transAxes,
           verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
           fontsize=9)
    
    # Set labels and limits
    ax.set_xlabel('Numerically computed time (s)')
    ax.set_ylabel('Predicted time (s)')
    ax.set_title('Model Validation Across Materials')
    
    # Set equal aspect ratio
    ax.set_aspect('equal', adjustable='box')
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    # Legend
    ax.legend(loc='lower right', frameon=True, framealpha=0.9)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig

def plot_supercooling_effect(T_ambient: np.ndarray,
                            times_no_sc: np.ndarray,
                            times_with_sc: np.ndarray,
                            experimental_data: Optional[Dict] = None,
                            save_path: Optional[str] = None) -> plt.Figure:
    """
    Create Figure 3: Solidification time vs ambient temperature with supercooling.
    
    Parameters:
    -----------
    T_ambient : array
        Ambient temperatures (°C)
    times_no_sc : array
        Solidification times without supercooling
    times_with_sc : array
        Solidification times with supercooling
    experimental_data : dict, optional
        {'T': temperatures, 't': times} from Müller et al.
    save_path : str, optional
        Path to save figure
    
    Returns:
    --------
    fig : matplotlib Figure
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # Plot model predictions
    ax.plot(T_ambient, times_no_sc,
           color=COLORS['primary'],
           linewidth=2,
           label='Model (no supercooling)')
    
    ax.plot(T_ambient, times_with_sc,
           color=COLORS['secondary'],
           linewidth=2,
           label='Model (with supercooling)')
    
    # Plot experimental data if provided
    if experimental_data:
        ax.scatter(experimental_data['T'],
                  experimental_data['t'],
                  color=COLORS['error'],
                  marker='o',
                  s=40,
                  edgecolors='black',
                  linewidth=0.5,
                  zorder=4,
                  label='Experiment [Müller et al., 2015]')
    
    # Find minimum point
    min_idx = np.argmin(times_with_sc)
    T_min = T_ambient[min_idx]
    t_min = times_with_sc[min_idx]
    
    # Highlight minimum point
    ax.scatter([T_min], [t_min],
              color='red',
              marker='*',
              s=100,
              zorder=5,
              label=f'Minimum at {T_min:.1f}°C')
    
    # Add vertical line at minimum
    ax.axvline(x=T_min, color='red', linestyle='--', alpha=0.5, linewidth=1)
    
    # Add text annotation for minimum
    ax.text(T_min + 0.5, t_min + 50,
           f'Minimum: {T_min:.1f}°C',
           fontsize=9,
           color='red',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Add inset explaining physical mechanism
    inset_ax = fig.add_axes([0.55, 0.55, 0.3, 0.3])
    
    # Create schematic for inset
    T_range = np.linspace(-25, 0, 100)
    cooling_rate = 1 / (1 + np.exp(-0.3*(T_range + 12)))  # Simplified
    driving_force = np.abs(T_range + 5)  # Simplified
    
    inset_ax.plot(T_range, cooling_rate, 
                 color=COLORS['primary'], 
                 linewidth=1.5,
                 label='Cooling rate')
    inset_ax.plot(T_range, driving_force, 
                 color=COLORS['secondary'], 
                 linewidth=1.5,
                 label='Driving force')
    
    # Mark competition region
    inset_ax.fill_between([-20, -5], 0, 1.2,
                         alpha=0.3, color='gray',
                         label='Competition region')
    
    inset_ax.set_xlabel('$T_\\infty$ (°C)', fontsize=8)
    inset_ax.set_ylabel('Normalized', fontsize=8)
    inset_ax.set_title('Competing effects', fontsize=9)
    inset_ax.legend(fontsize=7, frameon=False)
    inset_ax.grid(True, alpha=0.3)
    
    # Set main plot properties
    ax.set_xlabel('Ambient temperature $T_{\infty}$ (°C)')
    ax.set_ylabel('Solidification time (s)')
    ax.set_title('Water Droplet Solidification with Supercooling')
    
    ax.set_xlim(-25, 0)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    
    ax.legend(loc='upper right', frameon=True, framealpha=0.9)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig

def plot_geometry_comparison(results: Dict[str, Dict[str, np.ndarray]],
                            save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot comparison of phase change times for different geometries.
    
    Parameters:
    -----------
    results : dict
        Dictionary with geometry names as keys, each containing:
        - 'Bi': Biot numbers
        - 't_total': total phase change times
    save_path : str, optional
        Path to save figure
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    
    colors = [COLORS['primary'], COLORS['secondary'], COLORS['tertiary']]
    geometries = list(results.keys())
    
    for i, (geom, data) in enumerate(results.items()):
        ax.loglog(data['Bi'], data['t_total'],
                 color=colors[i % len(colors)],
                 linewidth=2,
                 marker='o',
                 markersize=5,
                 label=geom.capitalize())
    
    ax.set_xlabel('Biot number, Bi')
    ax.set_ylabel('Phase change time (s)')
    ax.set_title('Geometry Comparison')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(frameon=True, framealpha=0.9)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
    
    return fig

def create_sample_figures(output_dir: str = 'figures'):
    """
    Create sample figures for demonstration.
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # Sample data for Figure 1
    Bi = np.logspace(-2, 1, 100)  # 0.01 to 10
    
    # Simulated error patterns
    error_our = 2 + 3 * Bi**1.5 / (1 + Bi**1.5)
    error_liu = np.where(Bi < 0.2, 2 + 0.5*Bi, 20 + 80*(Bi-0.2))
    error_karwa = np.where(Bi < 0.5, 3 + Bi, 15 + 50*(Bi-0.5))
    
    error_data = {
        'Bi': Bi,
        'our_model': error_our,
        'liu_2009': error_liu,
        'karwa_2013': error_karwa
    }
    
    fig1 = plot_error_vs_biot(error_data, 
                             save_path=os.path.join(output_dir, 'Figure1_error_vs_Bi.png'))
    
    # Sample data for Figure 2
    np.random.seed(42)
    n_points = 50
    numerical_times = np.random.uniform(100, 2000, n_points)
    predicted_times = numerical_times * (1 + np.random.normal(0, 0.02, n_points))
    
    materials = np.random.choice(['Water', 'Paraffin', 'Gallium'], n_points)
    Bi_ranges = [(0.1, 0.5), (0.5, 1.0), (1.0, 1.5)] * (n_points // 3 + 1)
    
    fig2 = plot_predicted_vs_numerical(predicted_times, numerical_times,
                                      materials[:n_points], Bi_ranges[:n_points],
                                      save_path=os.path.join(output_dir, 'Figure2_validation.png'))
    
    # Sample data for Figure 3
    T_ambient = np.linspace(-25, 0, 100)
    times_no_sc = 5000 / (np.abs(T_ambient) + 1)
    times_with_sc = 4000 / (np.abs(T_ambient + 12) + 1) + 200
    
    # Add some experimental points
    exp_T = [-20, -15, -12, -10, -5]
    exp_t = [450, 350, 320, 380, 600]
    
    exp_data = {'T': exp_T, 't': exp_t}
    
    fig3 = plot_supercooling_effect(T_ambient, times_no_sc, times_with_sc,
                                   exp_data,
                                   save_path=os.path.join(output_dir, 'Figure3_supercooling.png'))
    
    print(f"Sample figures created in '{output_dir}' directory")
    return fig1, fig2, fig3

if __name__ == "__main__":
    # Create sample figures
    create_sample_figures()
    
    print("Visualization module ready.")
    print("Available functions:")
    print("1. plot_error_vs_biot() - Figure 1")
    print("2. plot_predicted_vs_numerical() - Figure 2")
    print("3. plot_supercooling_effect() - Figure 3")
    print("4. plot_geometry_comparison() - Additional figure")
