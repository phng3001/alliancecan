# P=NP
import os
import argparse
import sys

def check_pandas():
    """Check if Pandas is available."""
    try:
        global pd
        import pandas as pd
    except ImportError:
        print("Pandas is not available. Please install/load it.")
        sys.exit(1)

def check_plotly():
    """Check if Plotly is available."""
    try:
        global px
        import plotly.express as px
    except ImportError:
        print("Plotly is not available. Please install/load it.")
        sys.exit(1)

def check_required_columns(df, file_path):
    """Check if the required columns exist in the input DataFrame."""
    required_columns = {'#CHROM', 'START', 'NORMALIZED_LOG2_READS_RATIO_MUTANT/WT', 'GENE'}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        print(f"Error: Missing columns {missing_columns} in file {file_path}.")
        print("Required columns are: '#CHROM', 'START', 'NORMALIZED_LOG2_READS_RATIO_MUTANT/WT' and 'GENE'")
        sys.exit(1)

def read_data(input_files):
    """Read and concatenate input files into a single dataframe."""
    data_frames = []
    for file in input_files:
        sample_name = os.path.splitext(os.path.basename(file))[0]
        df = pd.read_csv(file, sep='\t')

        # Check for required columns
        check_required_columns(df, file)

        df['Sample'] = sample_name
        data_frames.append(df)
    return pd.concat(data_frames, ignore_index=True)

def inches_to_pixels(inches, dpi=96):
    """Converts inches to pixels based on a given DPI."""
    return int(inches * dpi)

def plot_interactive(data, output_dir, ylim, gridstep, dpi, width, height):
    """Generates interactive plots for each chromosome."""
    width = inches_to_pixels(width, dpi)
    height = inches_to_pixels(height, dpi)

    chromosomes = data['#CHROM'].unique()
    for chrom in chromosomes:
        chrom_data = data[data['#CHROM'] == chrom]
        fig = px.line(
            chrom_data,
            x='START',
            y='NORMALIZED_LOG2_READS_RATIO_MUTANT/WT',
            color='Sample',
            title=f'Chromosome {chrom}',
            labels={'START': 'Position',
                    'NORMALIZED_LOG2_READS_RATIO_MUTANT/WT': 'Normalized Log2 Ratio',
                    'GENE': 'Gene'},
            hover_data={'GENE': True}            
        )

        # Default grid step
        step = gridstep if gridstep else 1

        # y-axis layout
        yaxis_config = {
            'tickmode': 'linear',
            'dtick': step,  # Custom or default step size for gridlines
            'showgrid': True
        }

        if ylim:
            yaxis_config['range'] = ylim  # Add ylim if provided

        # Update layout
        fig.update_layout(
            yaxis=yaxis_config,
            xaxis=dict(showgrid=True),
            legend=dict(
                title='Samples',
                yanchor="top",
                y=1.02,
                xanchor="left",
                x=1.05
            ),
            xaxis_title='Genomic Coordinate',
            yaxis_title='Normalized Log2 Ratio',
            template="plotly_white",
            width=width,
            height=height
        )

        output_path = os.path.join(output_dir, f'{chrom}_interactive.html')
        fig.write_html(output_path)

def main():
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Generate interactive plots of the normalized log2 mutant/wild type read ratio for each chromosome.')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Input TSV files (one per sample)')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory to save plots')
    parser.add_argument('--ylim', nargs=2, type=float, default=None, help='y-axis limits as min max (e.g., --ylim -2 2)')
    parser.add_argument('--gridstep', type=float, default=None, help='Step size for y-axis gridlines. Default is 1')
    parser.add_argument('--dpi', type=int, default=96, help='Plot resolution in dpi (dots per inch). Default is 96')
    parser.add_argument('--width', type=float, default=12, help='Figure width in inches. Default is 12')
    parser.add_argument('--height', type=float, default=3, help='Figure height in inches. Default is 3')
    
    # Parse the arguments
    args = parser.parse_args()
    input = args.input
    outdir = args.outdir
    ylim = args.ylim
    gridstep = args.gridstep
    dpi = args.dpi
    width = args.width
    height = args.height

    # Check if the output directory exists
    if os.path.exists(outdir):
        print(f"Warning: The output directory {outdir} already exists. Exiting.")
        sys.exit(1)
    else:
        os.makedirs(outdir)

    # Check if Pandas and Plotly are available
    check_pandas()
    check_plotly()

    # Read and concatenate input files
    data = read_data(input)

    # Generates interactive plots for each chromosome
    plot_interactive(data, outdir, ylim, gridstep, dpi, width, height)

if __name__ == '__main__':
    # Check if no arguments are provided
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
