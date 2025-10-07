#!/usr/bin/env python3
"""Plot vector potential derivative errors from test_coil_tools_vector_potential_derivs."""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt


def load_error_data(filename):
    """Load error data from CSV file."""
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    return data[:, 0], data[:, 1]


def plot_errors(plot_dir):
    """Generate all error plots."""
    plots = [
        {
            'file': 'dAphi_dR_error.csv',
            'title': r'Relative error $|\Delta(dA_\phi/dR)|$',
            'xlabel': 'R coordinate (cm)',
            'output': 'coil_tools_dAphi_dR_error.png'
        },
        {
            'file': 'dAphi_dZ_error.csv',
            'title': r'Relative error $|\Delta(dA_\phi/dZ)|$',
            'xlabel': 'Z coordinate (cm)',
            'output': 'coil_tools_dAphi_dZ_error.png'
        }
    ]

    for plot_spec in plots:
        input_file = os.path.join(plot_dir, plot_spec['file'])
        output_file = os.path.join(plot_dir, plot_spec['output'])

        if not os.path.exists(input_file):
            print(f"Warning: {input_file} not found, skipping")
            continue

        coords, errors = load_error_data(input_file)

        plt.figure(figsize=(8, 6))
        plt.plot(coords, errors, 'b-', linewidth=1.5)
        plt.title(plot_spec['title'], fontsize=14)
        plt.xlabel(plot_spec['xlabel'], fontsize=12)
        plt.ylabel('Relative error', fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_file, dpi=150)
        plt.close()
        print(f"Saved plot to {output_file}")


if __name__ == '__main__':
    if len(sys.argv) > 1:
        plot_directory = sys.argv[1]
    else:
        plot_directory = os.getenv('COIL_TOOLS_PLOT_DIR', './plots')

    if not os.path.exists(plot_directory):
        print(f"Error: plot directory {plot_directory} does not exist")
        sys.exit(1)

    plot_errors(plot_directory)
    print("All plots generated successfully")
