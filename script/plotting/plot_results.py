# scripts/plotting/plot_results.py
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path

from loaders import load_all_results
from summaries import compute_teps_summary   # assuming this exists from previous step

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


def plot_teps_comparison():
    df = load_all_results(results_root="./results")
    
    # Safety filter: only positive TEPS
    df = df[df['TEPS'] > 0].copy()
    
    summary = compute_teps_summary(df)
    
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
                    markersize=7, label=label)
        
        # X-axis: log scale but show clean powers-of-two labels
        ax.set_xscale('log', base=2)
        unique_cores = sorted(data['cores'].unique())
        ax.set_xticks(unique_cores)
        ax.set_xticklabels([str(c) for c in unique_cores])
        ax.tick_params(axis='x', rotation=0)
        
        # Y-axis improvements
        ax.set_ylabel("TEPS (harmonic mean)", fontsize=12)
        # Add scientific notation hint in title / label area if needed
        # Most people prefer just clean number + unit, but we can force sci notation:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
        # This puts e.g. ×10⁸ in the top-left-ish area automatically
        
        ax.set_xlabel("Number of cores", fontsize=12)
        ax.set_title(f"TEPS comparison – {graph_name}", fontsize=13, pad=12)
        
        ax.legend(frameon=True, fontsize=10.5)
        ax.grid(True, which="major", alpha=0.5)
        ax.grid(True, which="minor", alpha=0.15, linestyle=":")
        
        save_plot(f"TEPS_comparison_{graph_name}")


def main():
    plot_teps_comparison()


if __name__ == "__main__":
    main()