using DataFrames
using CSV
using Plots
using Statistics
using Optim
using SpecialFunctions

"""
Calculate the logarithm of the binomial coefficient.
"""
function log_choose(x, y)
    return lgamma(x + 1) - lgamma(y + 1) - lgamma(x - y + 1)
end

"""
Calculate the log probability of a given weight and price.
"""
function log_P(w, p, r, n)
    if n == 0
        return log(w + (1-w) * (1-p)^r)
    else
        return log(1-w) + 
            log_choose(n+r-1, n) + 
            n * log(p) + 
            r * log(1-p)
    end
end

"""
Calculate the log likelihood of a given weight and price.
"""
function log_likelihood(w, p, r, count_dist)
    result = 0.0
    for n in keys(count_dist)
        result += count_dist[n] * log_P(w, p, r, n)
    end
    return result
end

"""
Define the negative log-likelihood function (to minimize).
"""
function neg_log_likelihood(params, count_dist)
    w, p, r = params
    # Add constraints to ensure parameters are in valid ranges
    if !(0 <= w <= 1 && 0 <= p <= 1 && r > 0)
        return Inf
    end
    return -log_likelihood(w, p, r, count_dist)
end

"""
Format a number to 4 decimal places as a string.
"""
function format_decimal(x)
    return string(round(x, digits=4))
end

"""
Plot the empirical distribution vs. fitted model.
"""
function plot_fit(count_dist, w, p, r, timepoint, condition; max_count=50)
    # Calculate theoretical probabilities
    x_range = 0:max_count
    theoretical_probs = zeros(length(x_range))
    for (i, k) in enumerate(x_range)
        theoretical_probs[i] = exp(log_P(w, p, r, k))
    end
    
    # Calculate empirical probabilities
    total_counts = sum(values(count_dist))
    empirical_probs = Dict(k => v/total_counts for (k, v) in count_dist)
    
    # Create plot (use plt instead of p to avoid parameter name collision)
    plt = plot(size=(800, 500), legend=:topright)
    
    # Plot empirical data as bars
    bar!(plt, collect(keys(empirical_probs)), collect(values(empirical_probs)), 
         alpha=0.5, label="Empirical", color=:blue)
    
    # Plot theoretical curve
    plot!(plt, x_range, theoretical_probs, line=:solid, marker=:circle, 
          markersize=3, color=:red, label="Fitted model")
    
    # Add details to plot
    xlabel!(plt, "Count")
    ylabel!(plt, "Probability")
    title!(plt, "Zero-Inflated Negative Binomial Fit - $(timepoint)h $(condition)")
    
    # Display parameter values using string interpolation
    w_formatted = format_decimal(w)
    r_formatted = format_decimal(r)
    p_formatted = format_decimal(p)
    
    textstr = "w (zero-inflation) = $(w_formatted)\n" *
              "r (size) = $(r_formatted)\n" *
              "p (prob) = $(p_formatted)"
    
    # Calculate annotation position
    x_pos = max_count * 0.95
    y_pos = maximum(theoretical_probs) * 0.95
    
    annotate!(plt, x_pos, y_pos, text(textstr, :right, 8))
    
    # Save the plot
    output_dir = joinpath(dirname(@__FILE__), "fish_quant_analysis/fits")
    mkpath(output_dir)
    savefig(plt, joinpath(output_dir, "fit_$(timepoint)h_$(condition).png"))
    
    return plt
end

"""
Fit model to a single count file and return the results.
"""
function fit_model(count_file)
    # Extract condition and time from filename
    filename = basename(count_file)
    # Extract timepoint and condition from filename (format: counts_{tp}h_{condition}.csv)
    parts = split(replace(filename, ".csv" => ""), "_")
    timepoint = replace(parts[2], "h" => "")
    condition = parts[3]
    
    # Read data
    counts = CSV.read(count_file, DataFrame)
    
    # Check if file exists and has data
    if isempty(counts)
        println("  Skipping empty file: $(filename)")
        return nothing
    end
    
    # Process data
    counts_spots = filter(row -> row.type in ["Mature", "Nascent"], counts)
    counts_no_spots = filter(row -> row.type == "No Spots", counts)
    
    # Skip if no spot data
    if isempty(counts_spots)
        println("  Skipping file with no spot data: $(filename)")
        return nothing
    end
    
    # Group by cell and frame
    grouped = combine(groupby(counts_spots, [:cell, :frame]), :count => sum => :count)
    
    # Create count distribution dictionary
    counts_dict = Dict{Int, Int}()
    for row in eachrow(grouped)
        count = row.count
        counts_dict[count] = get(counts_dict, count, 0) + 1
    end
    
    # Add zero counts
    counts_dict[0] = nrow(counts_no_spots)
    
    # If we have too few data points, skip fitting
    if length(counts_dict) < 4
        println("  Skipping file with insufficient data points: $(filename)")
        return nothing
    end
    
    # Sort keys for consistent processing
    count_keys = sort(collect(keys(counts_dict)))
    count_distribution = Dict(k => counts_dict[k] for k in count_keys)
    
    # Initial guess
    initial_guess = [0.5, 0.1, 2.0]
    
    # Perform optimization
    result = optimize(params -> neg_log_likelihood(params, count_distribution), 
                      initial_guess, 
                      NelderMead(), 
                      Optim.Options(iterations=10000))
    
    # Check if optimization converged
    if !Optim.converged(result)
        println("  Warning: Optimization did not converge for $(filename)")
    end
    
    # Extract optimized parameters
    w_opt, p_opt, r_opt = Optim.minimizer(result)
    
    # Check for unreasonable parameters
    if p_opt < 0.00001 || p_opt > 0.99999 || r_opt < 0.1 || r_opt > 1000
        println("  Skipping file with unreasonable fit parameters: $(filename)")
        return nothing
    end
    
    # Calculate negative log-likelihood at optimum
    nll = Optim.minimum(result)
    
    # Calculate AIC and BIC
    n_samples = sum(values(count_distribution))
    n_params = 3  # w, p, r
    aic = 2 * nll + 2 * n_params
    bic = 2 * nll + n_params * log(n_samples)
    
    # Calculate mean and variance based on the model
    mean_val = (1 - w_opt) * r_opt * (1 - p_opt) / p_opt
    variance = (1 - w_opt) * r_opt * (1 - p_opt) / (p_opt^2)
    
    # Plot fit
    plot_fit(count_distribution, w_opt, p_opt, r_opt, timepoint, condition)
    
    return Dict(
        "filename" => filename,
        "timepoint" => parse(Int, timepoint),
        "condition" => condition,
        "w" => w_opt,
        "p" => p_opt,
        "r" => r_opt,
        "nll" => nll,
        "aic" => aic,
        "bic" => bic,
        "mean" => mean_val,
        "variance" => variance,
        "converged" => Optim.converged(result),
        "n_samples" => n_samples,
        "n_datapoints" => length(counts_dict)
    )
end

"""
Plot parameter vs timepoint, grouped by condition.
"""
function plot_parameter_vs_timepoint(results_df, param, title)
    plt = plot(size=(800, 500), legend=:topright)
    
    # Get unique conditions
    conditions = unique(results_df.condition)
    
    # Plot each condition
    for condition in conditions
        subset = filter(row -> row.condition == condition, results_df)
        # Sort by timepoint
        sort!(subset, :timepoint)
        
        # Access columns directly by name as strings
        timepoints = subset[!, "timepoint"]
        param_values = subset[!, param]  # Use direct column access
        
        plot!(plt, timepoints, param_values, marker=:circle, label=condition)
    end
    
    xlabel!(plt, "Timepoint (h)")
    ylabel!(plt, param)
    title!(plt, title)
    
    # Save the plot
    output_dir = joinpath(dirname(@__FILE__), "fish_quant_analysis/summary_plots")
    mkpath(output_dir)
    savefig(plt, joinpath(output_dir, "$(param)_vs_timepoint.png"))
    
    return plt
end

"""
Main function to run the data processing pipeline.
"""
function main()
    # Get all count files
    fish_quant_dir = joinpath(dirname(@__FILE__), "fish_quant_analysis")
    count_files = filter(f -> startswith(basename(f), "counts_") && endswith(f, ".csv"), 
                         readdir(fish_quant_dir, join=true))
    
    if isempty(count_files)
        println("No count files found in fish_quant_analysis directory")
        return
    end
    
    println("Found $(length(count_files)) count files")
    
    # Create output directory for results
    results_dir = joinpath(dirname(@__FILE__), "fish_quant_analysis/results")
    mkpath(results_dir)
    
    # Fit model to each file
    results = []
    for file in count_files
        println("Processing $(basename(file))...")
        try
            result = fit_model(file)
            if result !== nothing
                push!(results, result)
                println("  Fitted parameters: w=$(format_decimal(result["w"])), p=$(format_decimal(result["p"])), r=$(format_decimal(result["r"]))")
            else
                println("  Skipped due to data issues")
            end
        catch e
            println("  Error processing $file: $e")
            # Print stack trace for debugging
            println(stacktrace(catch_backtrace()))
        end
    end
    
    # Check if we have any results
    if isempty(results)
        println("No successful fits were performed. Check the errors above.")
        return
    end
    
    # Create DataFrame from results
    results_df = DataFrame(results)
    
    # Sort by condition and timepoint
    sort!(results_df, [:condition, :timepoint])
    
    # Save results to CSV
    CSV.write(joinpath(results_dir, "fit_results.csv"), results_df)
    
    # Display results table
    println("\nFit Results:")
    pretty_table = "| Timepoint | Condition | w | p | r | Mean | Variance | AIC | N | Converged |\n" *
                   "|-----------|-----------|---|---|---|------|----------|-----|---|----------|\n"
    
    for row in eachrow(results_df)
        w_str = format_decimal(row.w)
        p_str = format_decimal(row.p)
        r_str = format_decimal(row.r)
        mean_str = format_decimal(row.mean)
        var_str = format_decimal(row.variance)
        aic_str = format_decimal(row.aic)
        
        pretty_table *= "| $(row.timepoint) | $(row.condition) | $(w_str) | $(p_str) | $(r_str) | " *
                        "$(mean_str) | $(var_str) | $(aic_str) | $(row.n_samples) | $(row.converged) |\n"
    end
    println(pretty_table)
    
    # Create summary plots
    plot_parameter_vs_timepoint(results_df, "w", "Zero-inflation parameter")
    plot_parameter_vs_timepoint(results_df, "p", "Probability parameter")
    plot_parameter_vs_timepoint(results_df, "r", "Size parameter")
    plot_parameter_vs_timepoint(results_df, "mean", "Mean of fitted distribution")
    plot_parameter_vs_timepoint(results_df, "variance", "Variance of fitted distribution")
end

# Run the main function
main()