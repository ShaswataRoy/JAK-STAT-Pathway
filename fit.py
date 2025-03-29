import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from glob import glob
import math
import os
from scipy.optimize import minimize
from tabulate import tabulate  # Make sure to install this: pip install tabulate

def log_choose(x, y):
    return math.lgamma(x + 1) - math.lgamma(y + 1) - math.lgamma(x - y + 1)

def log_P(w, p, r, n):
    """
    Calculate the log probability of a given weight and price.
    """
    if n == 0:
        return math.log(w + (1-w) * (1-p)**r)
    else:
        return math.log(1-w) + \
            log_choose(n+r-1, n) + \
            n * math.log(p) + \
            r * math.log(1-p)

def log_likelihood(w, p, r, count_dist):
    """
    Calculate the log likelihood of a given weight and price.
    """
    result = 0
    for n in count_dist.index:
        result += count_dist[n] * log_P(w, p, r, n)
    return result

# Define the negative log-likelihood function (to minimize)
def neg_log_likelihood(params, count_dist):
    w, p, r = params
    # Add constraints to ensure parameters are in valid ranges
    if not (0 <= w <= 1 and 0 <= p <= 1 and r > 0):
        return float('inf')
    return -log_likelihood(w, p, r, count_dist)

def fit_model(count_file):
    """Fit model to a single count file and return the results"""
    # Extract condition and time from filename
    filename = os.path.basename(count_file)
    # Extract timepoint and condition from filename (format: counts_{tp}h_{condition}.csv)
    parts = filename.replace('.csv', '').split('_')
    timepoint = parts[1].replace('h', '')
    condition = parts[2]
    
    # Read data
    counts = pd.read_csv(count_file)
    
    # Process data
    counts_spots = counts.query('type in ["Mature","Nascent"]').groupby(['cell','frame']).agg({'count': 'sum'})
    counts_no_spots = counts.query('type in ["No Spots"]')
    count_distribution = counts_spots['count'].value_counts().sort_index()
    count_distribution[0] = len(counts_no_spots)
    count_distribution = count_distribution.sort_index()
    
    # Initial guess
    initial_guess = [0.5, 0.1, 2]
    
    # Perform optimization
    result = minimize(
        neg_log_likelihood, 
        initial_guess, 
        args=(count_distribution,),
        method='Nelder-Mead',
        options={'maxiter': 10000}
    )
    
    # Extract optimized parameters
    w_opt, p_opt, r_opt = result.x
    
    # Calculate negative log-likelihood at optimum
    nll = result.fun
    
    # Calculate AIC and BIC
    n_samples = count_distribution.sum()
    n_params = 3  # w, p, r
    aic = 2 * nll + 2 * n_params
    bic = 2 * nll + n_params * np.log(n_samples)
    
    # Calculate mean and variance based on the model
    mean = (1 - w_opt) * r_opt * (1 - p_opt) / p_opt
    variance = (1 - w_opt) * r_opt * (1 - p_opt) / (p_opt**2)
    
    # Plot fit
    plot_fit(count_distribution, w_opt, p_opt, r_opt, timepoint, condition)
    
    return {
        'filename': filename,
        'timepoint': timepoint,
        'condition': condition,
        'w': w_opt,
        'p': p_opt,
        'r': r_opt,
        'nll': nll,
        'aic': aic,
        'bic': bic,
        'mean': mean,
        'variance': variance,
        'converged': result.success
    }

def plot_fit(count_dist, w, p, r, timepoint, condition, max_count=50):
    """Plot the empirical distribution vs. fitted model."""
    # Calculate theoretical probabilities
    x_range = np.arange(0, max_count+1)
    theoretical_probs = np.zeros(len(x_range))
    for i, k in enumerate(x_range):
        theoretical_probs[i] = math.exp(log_P(w, p, r, k))
    
    # Calculate empirical probabilities
    empirical_probs = count_dist / count_dist.sum()
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    # Plot empirical data as bars
    plt.bar(empirical_probs.index, empirical_probs.values, alpha=0.5, label='Empirical')
    
    # Plot theoretical curve
    plt.plot(x_range, theoretical_probs, 'r-', label='Fitted model')
    
    # Add details to plot
    plt.xlabel('Count')
    plt.ylabel('Probability')
    plt.title(f'Zero-Inflated Negative Binomial Fit - {timepoint}h {condition}')
    plt.legend()
    plt.grid(alpha=0.3)
    
    # Display parameter values
    textstr = f'w (zero-inflation) = {w:.4f}\nr (size) = {r:.4f}\np (prob) = {p:.4f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.95, 0.95, textstr, transform=plt.gca().transAxes, 
             fontsize=10, verticalalignment='top', horizontalalignment='right', bbox=props)
    
    # Save the plot
    plt.tight_layout()
    output_dir = os.path.join(os.path.dirname(__file__), "fish_quant_analysis/fits")
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f"fit_{timepoint}h_{condition}.png"))
    plt.close()

def main():
    # Get all count files
    count_files = glob("fish_quant_analysis/counts_*.csv")
    
    if not count_files:
        print("No count files found in fish_quant_analysis directory")
        return
    
    print(f"Found {len(count_files)} count files")
    
    # Create output directory for results
    results_dir = os.path.join(os.path.dirname(__file__), "fish_quant_analysis/results")
    os.makedirs(results_dir, exist_ok=True)
    
    # Fit model to each file
    results = []
    for file in count_files:
        print(f"Processing {os.path.basename(file)}...")
        try:
            result = fit_model(file)
            results.append(result)
            print(f"  Fitted parameters: w={result['w']:.4f}, p={result['p']:.4f}, r={result['r']:.4f}")
        except Exception as e:
            print(f"  Error processing {file}: {e}")
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results)
    
    # Sort by timepoint and condition
    results_df['timepoint'] = pd.to_numeric(results_df['timepoint'])
    results_df = results_df.sort_values(['condition', 'timepoint'])
    
    # Save results to CSV
    results_df.to_csv(os.path.join(results_dir, "fit_results.csv"), index=False)
    
    # Display results table
    table_data = results_df[['timepoint', 'condition', 'w', 'p', 'r', 'mean', 'variance', 'aic', 'converged']]
    print("\nFit Results:")
    print(tabulate(table_data, headers='keys', tablefmt='pipe', floatfmt='.4f'))
    
    # Create summary plots
    plot_parameter_vs_timepoint(results_df, 'w', 'Zero-inflation parameter')
    plot_parameter_vs_timepoint(results_df, 'p', 'Probability parameter')
    plot_parameter_vs_timepoint(results_df, 'r', 'Size parameter')
    plot_parameter_vs_timepoint(results_df, 'mean', 'Mean of fitted distribution')
    plot_parameter_vs_timepoint(results_df, 'variance', 'Variance of fitted distribution')

def plot_parameter_vs_timepoint(results_df, param, title):
    """Plot parameter vs timepoint, grouped by condition."""
    plt.figure(figsize=(10, 6))
    
    # Get unique conditions
    conditions = results_df['condition'].unique()
    
    for condition in conditions:
        subset = results_df[results_df['condition'] == condition]
        plt.plot(subset['timepoint'], subset[param], marker='o', label=condition)
    
    plt.xlabel('Timepoint (h)')
    plt.ylabel(param)
    plt.title(title)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    output_dir = os.path.join(os.path.dirname(__file__), "fish_quant_analysis/summary_plots")
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f"{param}_vs_timepoint.png"))
    plt.close()

if __name__ == "__main__":
    main()