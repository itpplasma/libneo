#!/usr/bin/env python3
"""Plot Aphi gradient validation results."""

import sys
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def plot_error_profile(csv_file, output_file, title, xlabel):
    """Plot error profile from CSV file."""
    df = pd.read_csv(csv_file)

    plt.figure(figsize=(10, 6))
    plt.semilogy(df['coord'], df['error'], 'b-', linewidth=2)
    plt.axhline(y=1e-4, color='r', linestyle='--', label='0.01% tolerance')
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel('Relative error (log scale)', fontsize=12)
    plt.title(title, fontsize=14)
    plt.grid(True, alpha=0.3, which='both')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"Saved plot to {output_file}")


def main():
    if len(sys.argv) < 2:
        print("Usage: plot_aphi_gradients.py <plot_directory>")
        sys.exit(1)

    plot_dir = Path(sys.argv[1])

    plots = [
        {
            'file': plot_dir / 'aphi_dR_error_pointwise.csv',
            'title': r'Relative error in $\partial A_\phi/\partial R$ (point-by-point FD)',
            'xlabel': 'R coordinate (cm)',
            'output': plot_dir / 'aphi_dR_error_pointwise.png'
        },
        {
            'file': plot_dir / 'aphi_dZ_error_pointwise.csv',
            'title': r'Relative error in $\partial A_\phi/\partial Z$ (point-by-point FD)',
            'xlabel': 'Z coordinate (cm)',
            'output': plot_dir / 'aphi_dZ_error_pointwise.png'
        }
    ]

    for spec in plots:
        if spec['file'].exists():
            plot_error_profile(spec['file'], spec['output'], spec['title'], spec['xlabel'])

    print("All plots generated successfully")


if __name__ == '__main__':
    main()
