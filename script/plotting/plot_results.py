# scripts/plotting/plot_results.py
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path

from loaders import load_all_results
from summaries import compute_cv_avg, compute_max_over_mean_avg, compute_teps_harmonic, compute_times_p90   # assuming this exists from previous step

# A few constants used in graphs
VARIANT_COLORS = {
    'basic':  '#1f77b4',   
    'hybrid': "#f02525", 
}

TIME_COLORS = {
    'total': '#ffcc00',     
    'comm':  '#9932cc',     
}

BAR_WIDTH = 0.35              # width of each bar in grouped bar plot

plt.style.use("seaborn-v0_8-whitegrid")

SAVE_FOLDER = Path("./plots").resolve()
SAVE_FOLDER.mkdir(exist_ok=True, parents=True)


def save_plot(name: str, fig=None):
    if fig is None:
        fig = plt.gcf()
    path = SAVE_FOLDER / f"{name}.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    print(f"Saved → {path}")
    plt.close(fig)

# ==================================================================================================
# TEPS comparison
# ==================================================================================================
def plot_teps_comparison():
    df = load_all_results(results_root="./results")
    
    # Safety filter: only positive TEPS
    df = df[df['TEPS'] > 0].copy()
    
    summary = compute_teps_harmonic(df)
    
    # Get sorted unique graphs
    graphs = sorted(summary['graph'].unique())
    
    for graph_name in graphs:
        data = summary[summary['graph'] == graph_name].copy()
        
        if data.empty:
            continue
            
        # Sort by cores
        data = data.sort_values('cores')
        
        fig, ax = plt.subplots(figsize=(10, 5.8))
        
        # Plot both variants if present
        for variant in ['basic', 'hybrid']:
            sub = data[data['variant'] == variant]
            if sub.empty:
                continue
                
            cores = sub['cores']
            teps = sub['TEPS_harmonic']
            
            label = f"{variant.capitalize()}"
            ax.plot(cores, teps, marker='o', linestyle='-', linewidth=1.9,
                markersize=7, label=label,
                color=VARIANT_COLORS.get(variant, 'gray'))
        
        # X-axis: CORES
        ax.set_xscale('log', base=2)
        unique_cores = sorted(data['cores'].unique())
        ax.set_xticks(unique_cores)
        ax.set_xticklabels([str(c) for c in unique_cores])
        ax.tick_params(axis='x', rotation=0)
        
        # Y-axis: TEPS
        ax.set_ylabel("TEPS (harmonic mean)", fontsize=12)
        # Add scientific notation hint in title / label area if needed
        # Most people prefer just clean number + unit, but we can force sci notation:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
        # This puts e.g. ×10⁸ in the top-left-ish area automatically
        
        ax.set_xlabel("Number of cores", fontsize=12)
        title_name = "Kronecker generated" if graph_name == "weak_scaling" else graph_name
        ax.set_title(f"TEPS comparison - {title_name}", fontsize=13, pad=12)        

        ax.legend(frameon=True, fontsize=10.5)
        ax.grid(True, which="major", alpha=0.5)
        ax.grid(True, which="minor", alpha=0.15, linestyle=":")
        
        save_plot(f"TEPS_comparison_{graph_name}")


# ==================================================================================================
# TIME comparison
# ==================================================================================================
def plot_total_time_p90():
    df = load_all_results(results_root="./results")
    
    # Keep only valid positive times
    df = df[df['total_time'] > 0].copy()
    
    summary = compute_times_p90(df)
    
    graphs = sorted(summary['graph'].unique())
    
    for graph_name in graphs:
        data = summary[summary['graph'] == graph_name].copy()
        if data.empty:
            continue
            
        data = data.sort_values('cores')
        
        fig, ax = plt.subplots(figsize=(10, 5.8))
        
        for variant in ['basic', 'hybrid']:
            sub = data[data['variant'] == variant]
            if sub.empty:
                continue
                
            cores = sub['cores']
            time_p90 = sub['total_time_p90']
            
            label = f"{variant.capitalize()}"
            ax.plot(cores, time_p90, marker='o', linestyle='-', linewidth=1.9,
                    markersize=7, label=label,
                    color=VARIANT_COLORS.get(variant, 'gray'))
        
        ax.set_xscale('log', base=2)
        unique_cores = sorted(data['cores'].unique())
        ax.set_xticks(unique_cores)
        ax.set_xticklabels([str(c) for c in unique_cores])
        ax.tick_params(axis='x', rotation=0)
        
        ax.set_yscale('log')   # ← usually very useful for time plots
        ax.set_xlabel("Number of cores", fontsize=12)
        ax.set_ylabel("Total time per solution – p90 (s)", fontsize=12)
        ax.set_title(f"Runtime comparison (p90) – {graph_name}", fontsize=13, pad=12)
        
        ax.legend(frameon=True, fontsize=10.5)
        ax.grid(True, which="major", alpha=0.5)
        ax.grid(True, which="minor", alpha=0.15, linestyle=":")
        
        save_plot(f"runtime_p90_{graph_name}")


# ==================================================================================================
# TOTAL vs COMMUNICATION time comparison
# ==================================================================================================
def plot_times_p90_bars():
    df = load_all_results(results_root="./results")
    df = df[df['total_time'] > 0].copy()
    
    summary = compute_times_p90(df)
    
    graphs = sorted(summary['graph'].unique())
    
    for graph_name in graphs:
        data = summary[summary['graph'] == graph_name].copy()
        if data.empty:
            continue
        
        data = data.sort_values('cores')
        cores_list = sorted(data['cores'].unique())
        
        fig, ax = plt.subplots(figsize=(11, 6))
        
        x_positions = np.arange(len(cores_list))  # 0,1,2,... for each core count
        
        for i, variant in enumerate(['basic', 'hybrid']):
            sub = data[data['variant'] == variant]
            if sub.empty:
                continue
            
            # map cores to integer positions
            pos = [x_positions[list(cores_list).index(c)] + i*BAR_WIDTH for c in sub['cores']]
            
            # total time bars (background)
            ax.bar(pos, sub['total_time_p90'], width=BAR_WIDTH,
                   color=TIME_COLORS['total'], label=f"{variant.capitalize()} total" if i==0 else "",
                   edgecolor='black', linewidth=0.8)
            
            # comm time bars (on top / overlapping visually, but shifted)
            ax.bar(pos, sub['comm_time_p90'], width=BAR_WIDTH,
                   color=TIME_COLORS['comm'], label=f"{variant.capitalize()} comm" if i==0 else "",
                   edgecolor='black', linewidth=0.8, alpha=0.95)
        
        # clean x-ticks
        ax.set_xticks(x_positions + BAR_WIDTH/2)
        ax.set_xticklabels([str(c) for c in cores_list])
        
        ax.set_xlabel("Number of cores")
        ax.set_ylabel("90th percentile time (s)")
        title_name = "Kronecker generated" if graph_name == "weak_scaling" else graph_name
        ax.set_title(f"p90 Runtime Breakdown – {title_name}", fontsize=13, pad=12)
        
        ax.legend(frameon=True, fontsize=10, ncol=2)
        ax.grid(True, axis='y', alpha=0.35, linestyle="--")
        
        save_plot(f"runtime_p90_bars_{graph_name}")

def plot_max_over_mean():
    df = load_all_results(results_root="./results")
    
    # Filter valid positive values (safety, assuming >1 typical)
    df = df[df['max_over_mean'] > 0].copy()
    
    summary = compute_max_over_mean_avg(df)
    
    graphs = sorted(summary['graph'].unique())
    
    for graph_name in graphs:
        data = summary[summary['graph'] == graph_name].copy()
        if data.empty:
            continue
            
        data = data.sort_values('cores')
        
        fig, ax = plt.subplots(figsize=(10, 5.8))
        
        for variant in ['basic', 'hybrid']:
            sub = data[data['variant'] == variant]
            if sub.empty:
                continue
                
            cores = sub['cores']
            mom = sub['max_over_mean_avg']
            
            label = f"{variant.capitalize()}"
            ax.plot(cores, mom, marker='o', linestyle='-', linewidth=1.9,
                    markersize=7, label=label,
                    color=VARIANT_COLORS.get(variant, 'gray'))
        
        ax.set_xscale('log', base=2)
        unique_cores = sorted(data['cores'].unique())
        ax.set_xticks(unique_cores)
        ax.set_xticklabels([str(c) for c in unique_cores])
        ax.tick_params(axis='x', rotation=0)
        
        # Y linear (ratios usually start ~1)
        ax.set_xlabel("Number of cores", fontsize=12)
        ax.set_ylabel("Max/Mean Load Imbalance (mean over runs)", fontsize=12)
        
        title_name = "Kronecker generated" if graph_name == "weak_scaling" else graph_name
        ax.set_title(f"Load Imbalance (Max/Mean) – {title_name}", fontsize=13, pad=12)
        
        ax.legend(frameon=True, fontsize=10.5)
        ax.grid(True, which="major", alpha=0.5)
        ax.grid(True, which="minor", alpha=0.15, linestyle=":")
        
        save_plot(f"load_imbalance_max_over_mean_{graph_name}")


def plot_cv():
    df = load_all_results(results_root="./results")
    
    df = df[df['CV'] >= 0].copy()  # CV is non-negative
    
    summary = compute_cv_avg(df)
    
    graphs = sorted(summary['graph'].unique())
    
    for graph_name in graphs:
        data = summary[summary['graph'] == graph_name].copy()
        if data.empty:
            continue
            
        data = data.sort_values('cores')
        
        fig, ax = plt.subplots(figsize=(10, 5.8))
        
        for variant in ['basic', 'hybrid']:
            sub = data[data['variant'] == variant]
            if sub.empty:
                continue
                
            cores = sub['cores']
            cv = sub['CV_avg']
            
            label = f"{variant.capitalize()}"
            ax.plot(cores, cv, marker='o', linestyle='-', linewidth=1.9,
                    markersize=7, label=label,
                    color=VARIANT_COLORS.get(variant, 'gray'))
        
        ax.set_xscale('log', base=2)
        unique_cores = sorted(data['cores'].unique())
        ax.set_xticks(unique_cores)
        ax.set_xticklabels([str(c) for c in unique_cores])
        ax.tick_params(axis='x', rotation=0)
        
        ax.set_xlabel("Number of cores", fontsize=12)
        ax.set_ylabel("Coefficient of Variation (mean over runs)", fontsize=12)
        
        title_name = "Kronecker generated" if graph_name == "weak_scaling" else graph_name
        ax.set_title(f"Load Variation (CV) – {title_name}", fontsize=13, pad=12)
        
        ax.legend(frameon=True, fontsize=10.5)
        ax.grid(True, which="major", alpha=0.5)
        ax.grid(True, which="minor", alpha=0.15, linestyle=":")
        
        save_plot(f"load_variation_cv_{graph_name}")

def main():
    plot_teps_comparison()
    plot_total_time_p90()
    plot_times_p90_bars()
    plot_max_over_mean() 
    plot_cv()              

if __name__ == "__main__":
    main()