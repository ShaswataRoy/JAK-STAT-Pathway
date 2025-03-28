# Import necessary packages
using Catalyst
using DifferentialEquations
using Plots
using Statistics  # For mean, var functions
using StatsPlots

theme(:dao)

# Define the two-state gene expression model with Catalyst
gene_model = @reaction_network begin
    k_on, DNA_off --> DNA_on     # DNA activation
    k_off, DNA_on --> DNA_off    # DNA inactivation
    k_tx, DNA_on --> DNA_on + mRNA   # Transcription (only in ON state)
    d_mrna, mRNA --> ∅           # mRNA degradation
end

# Function to simulate the gene expression model
function simulate_gene_expression(k_on, k_off, k_tx, d_mrna, initial_state, tspan)
    # Parameter values
    p = [k_on, k_off, k_tx, d_mrna]
    
    # Initial state [DNA_off, DNA_on, mRNA]
    u0 = initial_state
    
    # Define the Jump Problem (for discrete stochastic simulation)
    discrete_prob = DiscreteProblem(gene_model, u0, tspan, p)
    jump_prob = JumpProblem(gene_model, discrete_prob, Direct())
    
    # Solve the jump problem (simulate the process)
    sol = solve(jump_prob, SSAStepper())
    
    return sol
end

# Function to properly plot a step function for discrete state changes
function step_plot(t, y; kwargs...)
    # Create vectors for proper step drawing
    n = length(t)
    t_step = Vector{Float64}(undef, 2n)
    y_step = Vector{Float64}(undef, 2n)
    
    # Handle first point
    t_step[1] = t[1]
    y_step[1] = y[1]
    
    # Create step pattern
    for i in 1:n-1
        t_step[2i] = t[i+1]
        t_step[2i+1] = t[i+1]
        y_step[2i] = y[i]
        y_step[2i+1] = y[i+1]
    end
    
    # Handle last point (remove extra point)
    t_step = t_step[1:end-1]
    y_step = y_step[1:end-1]
    
    return plot(t_step, y_step; kwargs...)
end

# Function to analyze results
function analyze_results(sol, regime_name)
    # Extract data
    dna_on = sol[2,:]
    mrna = sol[3,:]
    
    # Calculate statistics
    mean_dna_on = mean(dna_on)
    fraction_time_on = mean_dna_on
    mean_mrna = mean(mrna)
    var_mrna = var(mrna)
    fano_factor = var_mrna / mean_mrna
    
    println("\nAnalysis Results for $regime_name regime:")
    println("Fraction of time DNA is ON: $(round(fraction_time_on, digits=3))")
    println("Mean mRNA level: $(round(mean_mrna, digits=2))")
    println("Variance of mRNA level: $(round(var_mrna, digits=2))")
    println("Fano factor (variance/mean): $(round(fano_factor, digits=2))")
        
    # Return analysis results as dictionary
    results = Dict(
        "fraction_on" => fraction_time_on,
        "mean_mrna" => mean_mrna,
        "var_mrna" => var_mrna,
        "fano_factor" => fano_factor
    )
    
    return results
end

# Simple animation function with minimal GKS operations
function create_simple_animation(sol, title, filename, max_time=200.0)
    # This function creates an animation with fewer GKS operations
    # to avoid the "too many open files" error

    # Find time points up to max_time
    time_idx = findall(t -> t <= max_time, sol.t)
    
    # Maximum mRNA level for consistent y-axis
    max_mrna = maximum(sol[3, time_idx]) * 1.1
    max_mrna = max(5, max_mrna)  # Ensure minimum height
    
    # Create a reasonable number of frames (fewer = less likely to have errors)
    frames = 90
    
    # Create the animation
    println("Creating $title animation...")
    anim = @animate for i in 1:frames
        # Find the appropriate index
        frame_time = i * (max_time / frames)
        idx = findlast(t -> t <= frame_time, sol.t)
        
        if isnothing(idx)
            idx = 1
        end
        
        # Current state
        current_dna_state = sol[2, idx]
        
        # Create a simple plot layout
        p = plot(layout=(2,1), size=(800, 600), 
                title="$title (Time: $(round(sol.t[idx], digits=1)))")
        
        # DNA state plot (top)
        plot!(p[1], sol.t[1:idx], sol[2, 1:idx], 
              label="DNA State", 
              color=:green, 
              lw=2,
              xlabel="",
              ylabel="State (0/1)",
              ylims=(-0.1, 1.1),
              xlims=(0, max_time),
              seriestype=:steppost)  # Use built-in step plot
        
        # Add DNA state indicator
        if current_dna_state > 0
            annotate!(p[1], 10, 0.5, text("DNA: ON", :green, 12))
        else
            annotate!(p[1], 10, 0.5, text("DNA: OFF", :red, 12))
        end
        
        # mRNA level plot (bottom)
        plot!(p[2], sol.t[1:idx], sol[3, 1:idx], 
              label="mRNA", 
              color=:blue, 
              lw=2,
              xlabel="Time",
              ylabel="mRNA Count",
              ylims=(0, max_mrna),
              xlims=(0, max_time),
              seriestype=:steppost)  # Use built-in step plot
        
        # Current frame
        p
    end
    
    # Save the animation
    try
        gif(anim, filename, fps=15)
        println("Animation saved to $filename")
    catch e
        println("Error saving animation: $e")
        println("Trying with fewer frames...")
        
        # Create a reduced animation with even fewer frames as fallback
        anim_small = @animate for i in 1:30
            frame_time = i * (max_time / 30)
            idx = findlast(t -> t <= frame_time, sol.t)
            
            if isnothing(idx)
                idx = 1
            end
            
            # Create a simple plot
            plot(layout=(2,1), size=(800, 600), 
                 title="$title (Time: $(round(sol.t[idx], digits=1)))")
            
            # DNA state plot
            plot!(sol.t[1:idx], sol[2, 1:idx], 
                  subplot=1, 
                  label="DNA State", 
                  color=:green, 
                  lw=2,
                  ylims=(-0.1, 1.1),
                  xlims=(0, max_time),
                  seriestype=:steppost)
            
            # mRNA level plot
            plot!(sol.t[1:idx], sol[3, 1:idx], 
                  subplot=2, 
                  label="mRNA", 
                  color=:blue, 
                  lw=2,
                  ylims=(0, max_mrna),
                  xlims=(0, max_time),
                  seriestype=:steppost)
        end
        
        gif(anim_small, "small_$filename", fps=10)
        println("Reduced animation saved to small_$filename")
    end
    
    return filename
end

# Function to create a smoother animation with consistent time sampling
function create_smooth_animation(sol, title, filename, max_time=400.0)
    # This function creates a smoother animation by using interpolation
    # for equally spaced time points
    
    println("Creating smooth $title animation...")
    
    # Maximum mRNA level for consistent y-axis
    time_idx = findall(t -> t <= max_time, sol.t)
    max_mrna = maximum(sol[3, time_idx]) * 1.1
    max_mrna = max(5, max_mrna)  # Ensure minimum height
    
    # Number of frames for animation
    frames = 240
    
    # Create evenly spaced time points for smooth animation
    even_times = range(0, max_time, length=frames)
    
    # Create the animation
    anim = @animate for i in 1:frames
        current_time = even_times[i]
        
        # Find the two indices that bracket this time
        idx_before = findlast(t -> t <= current_time, sol.t)
        if isnothing(idx_before)
            idx_before = 1
        end
        
        # Get the current state (use exact simulation points to preserve step function)
        # For DNA state, use the exact value at idx_before
        current_dna_state = sol[2, idx_before]
        
        # Create the plot
        p = plot(layout=(2,1), size=(1800, 600))
        
        # Plot all data points up to the current time
        # First find all indices up to current time
        plot_idx = findall(t -> t <= current_time, sol.t)
        
        # DNA state plot (top)
        plot!(p[1], sol.t[plot_idx], sol[2, plot_idx], 
              color=:green, 
              legend = false,
              lw=2,
              xlabel="",
              ylabel="DNA State",
              ylims=(-0.1, 1.3),
              xlims=(0, max_time),
              seriestype=:steppost)  # Use built-in step plot

        # title!("$title (Time: $(round(current_time, digits=1)))")
        title!("")
        
        # Add DNA state indicator
        if current_dna_state > 0
            annotate!(p[1], max_time * 0.1, 1.2, text("DNA: ON", :green, 14))
        else
            annotate!(p[1], max_time * 0.1, 1.2, text("DNA: OFF", :red, 14))
        end
        
        # mRNA level plot (bottom)
        plot!(p[2], sol.t[plot_idx], sol[3, plot_idx],
              color=:blue, 
              lw=2,
              legend = false,
              xlabel="Time",
              ylabel="mRNA Count",
              ylims=(0, max_mrna),
              xlims=(0, max_time),
              seriestype=:steppost)  # Use built-in step plot
        
        # Add current time marker (vertical line)
        vline!(p[1], [current_time], label="", color=:black, linestyle=:dash, alpha=0.5)
        vline!(p[2], [current_time], label="", color=:black, linestyle=:dash, alpha=0.5)
        
        # Current frame
        p
    end
    
    # Save the animation with higher frame rate for smoother playback
    try
        gif(anim, filename, fps=24)
        println("Animation saved to $filename")
    catch e
        println("Error saving animation: $e")
        println("Trying with fewer frames...")
        
        # Create a reduced version as fallback
        reduced_frames = 60
        reduced_times = range(0, max_time, length=reduced_frames)
        
        anim_small = @animate for i in 1:reduced_frames
            current_time = reduced_times[i]
            
            # Find index before this time
            idx_before = findlast(t -> t <= current_time, sol.t)
            if isnothing(idx_before)
                idx_before = 1
            end
            
            # Get current state
            current_dna_state = sol[2, idx_before]
            
            # Find all points up to current time
            plot_idx = findall(t -> t <= current_time, sol.t)
            
            # Create simple plot
            p = plot(layout=(2,1), size=(800, 600), 
                   title="$title (Time: $(round(current_time, digits=1)))")
            
            # DNA plot
            plot!(p[1], sol.t[plot_idx], sol[2, plot_idx], 
                  label="DNA", color=:green, lw=2,
                  ylims=(-0.1, 1.1), xlims=(0, max_time),
                  seriestype=:steppost)
            
            # mRNA plot  
            plot!(p[2], sol.t[plot_idx], sol[3, plot_idx], 
                  label="mRNA", color=:blue, lw=2,
                  ylims=(0, max_mrna), xlims=(0, max_time),
                  seriestype=:steppost)
            
            p
        end
        
        gif(anim_small, "small_$filename", fps=15)
        println("Reduced animation saved to small_$filename")
    end
    
    return filename
end

# Function to run simulations and create animations
function run_simulations_and_animations()
    # Reset the plotting backend to try to avoid GKS errors
    gr(fmt=:png)  # Use PNG format which tends to be more stable
    
    tspan = (0.0, 1000.0)
    initial_state = [0, 1, 0]  # Start with DNA ON and no mRNA
    
    # Common parameters
    d_mrna = 0.05    # mRNA degradation rate
    
    # Fast switching parameters
    k_on_fast = 0.05
    k_off_fast = 0.05
    k_tx_fast = 1.0
    
    # Slow switching parameters
    k_on_slow = 0.05
    k_off_slow = 0.1
    k_tx_slow = 2.0  # Higher to keep similar mean mRNA level
    
    # Theoretical steady-state mRNA levels
    theoretical_mean_fast = (k_on_fast/(k_on_fast+k_off_fast)) * (k_tx_fast/d_mrna)
    theoretical_mean_slow = (k_on_slow/(k_on_slow+k_off_slow)) * (k_tx_slow/d_mrna)
    
    println("\n--- Theoretical Predictions ---")
    println("Fast switching theoretical mean mRNA: $(round(theoretical_mean_fast, digits=2))")
    println("Slow switching theoretical mean mRNA: $(round(theoretical_mean_slow, digits=2))")
    
    # Run simulations
    println("\nRunning fast switching simulation...")
    sol_fast = simulate_gene_expression(k_on_fast, k_off_fast, k_tx_fast, d_mrna, initial_state, tspan)
    
    println("\nRunning slow switching simulation...")
    sol_slow = simulate_gene_expression(k_on_slow, k_off_slow, k_tx_slow, d_mrna, initial_state, tspan)
    
    # Analyze results
    results_fast = analyze_results(sol_fast, "Faster Switching")
    results_slow = analyze_results(sol_slow, "Higher Transcription Rate")
    
    # Create smooth animations one at a time
    create_smooth_animation(sol_fast, "Faster Switching", "fast_switching_smooth.gif")
    GC.gc()  # Force garbage collection
    
    create_smooth_animation(sol_slow, "Higher Transcription Rate", "slow_switching_smooth.gif")
    GC.gc()  # Force garbage collection
    
    # Also create static images as a backup in case animations fail
    create_static_image_comparison(sol_fast, sol_slow)
    
    return sol_fast, sol_slow, results_fast, results_slow
end

# Function to create static image comparison as a backup
function create_static_image_comparison(sol_fast, sol_slow)
    # Time window for comparison
    max_time = 200.0
    
    # Find indices up to max_time
    fast_idx = findall(t -> t <= max_time, sol_fast.t)
    slow_idx = findall(t -> t <= max_time, sol_slow.t)
    
    # Plot DNA states
    p1 = plot(sol_fast.t[fast_idx], sol_fast[2, fast_idx],
              title="Fast Switching: DNA State",
              label="DNA ON/OFF",
              color=:green,
              lw=2,
              ylims=(-0.1, 1.1),
              seriestype=:steppost)
    
    p2 = plot(sol_slow.t[slow_idx], sol_slow[2, slow_idx],
              title="Slow Switching: DNA State",
              label="DNA ON/OFF",
              color=:green,
              lw=2,
              ylims=(-0.1, 1.1),
              seriestype=:steppost)
    
    # Plot mRNA levels
    p3 = plot(sol_fast.t[fast_idx], sol_fast[3, fast_idx],
              title="Fast Switching: mRNA",
              label="mRNA",
              color=:blue,
              lw=2,
              seriestype=:steppost)
    
    p4 = plot(sol_slow.t[slow_idx], sol_slow[3, slow_idx],
              title="Slow Switching: mRNA",
              label="mRNA",
              color=:blue,
              lw=2,
              seriestype=:steppost)
    
    # Create a 2x2 comparison plot
    full_comparison = plot(p1, p2, p3, p4, 
                          layout=(2,2), 
                          size=(1000, 800),
                          title="Fast vs. Slow Promoter Switching")
    
    # Save images
    savefig(full_comparison, "switching_comparison.png")
    
    println("Static comparison image saved to switching_comparison.png")
    return full_comparison
end

# Function to create histograms of mRNA levels
function plot_mrna_histograms(sol_fast, sol_slow)
    # Use the second half of the simulation time to ensure steady state
    halfway_fast = div(length(sol_fast.t), 2)
    halfway_slow = div(length(sol_slow.t), 2)
    
    # Extract mRNA counts after reaching steady state
    mrna_fast = sol_fast[3, halfway_fast:end]
    mrna_slow = sol_slow[3, halfway_slow:end]
    
    # Calculate theoretical Poisson for fast switching
    mean_mrna_fast = mean(mrna_fast)
    max_count = max(maximum(mrna_fast), maximum(mrna_slow)) + 5
    x_range = 0:1:max_count
    poisson_pmf = [exp(-mean_mrna_fast) * mean_mrna_fast^k / factorial(big(k)) for k in x_range]
    
    # Create histogram for fast switching
    p1 = histogram(mrna_fast, 
                  bins=0:1:max_count,
                  normalize=:probability,
                  title="Fast Switching mRNA Distribution",
                  xlabel="mRNA Count",
                  ylabel="Probability",
                  label="Simulation",
                  alpha=0.7,
                  color=:blue)
    
    # Add theoretical Poisson distribution for comparison
    plot!(p1, x_range, poisson_pmf, 
          line=:stem, marker=:circle, 
          label="Poisson (λ=$(round(mean_mrna_fast, digits=1)))",
          color=:red)
    
    # Create histogram for slow switching
    p2 = histogram(mrna_slow, 
                  bins=0:1:max_count,
                  normalize=:probability,
                  title="Slow Switching mRNA Distribution",
                  xlabel="mRNA Count",
                  ylabel="Probability",
                  label="Simulation",
                  alpha=0.7,
                  color=:green)
    
    # Create combined plot for direct comparison
    p3 = histogram(mrna_fast,
                  bins=0:1:max_count, 
                  normalize=:probability,
                  title="Distribution Comparison",
                  xlabel="mRNA Count",
                  ylabel="Probability",
                  label="Fast Switching",
                  alpha=0.5,
                  color=:blue)
                  
    histogram!(p3, mrna_slow, 
              bins=0:1:max_count, 
              normalize=:probability,
              label="Slow Switching",
              alpha=0.5,
              color=:green)
    
    # Create a 2x2 layout with the histograms and a text box
    layout = @layout [a b; c d]
    
    # Statistical info
    fast_fano = var(mrna_fast) / mean(mrna_fast)
    slow_fano = var(mrna_slow) / mean(mrna_slow)
    
    stats_text = """
    Fast Switching:
        Mean: $(round(mean(mrna_fast), digits=2))
        Variance: $(round(var(mrna_fast), digits=2))
        Fano factor: $(round(fast_fano, digits=2))
    
    Slow Switching:
        Mean: $(round(mean(mrna_slow), digits=2))
        Variance: $(round(var(mrna_slow), digits=2))
        Fano factor: $(round(slow_fano, digits=2))
    
    Noise amplification: $(round(slow_fano/fast_fano, digits=2))x
    """
    
    # Create a plot with only text
    p4 = plot(title="Statistical Comparison", grid=false, showaxis=false, ticks=false)
    annotate!(p4, 0.5, 0.5, text(stats_text, 10))
    
    # Combine all plots
    histogram_plot = plot(p1, p2, p3, p4, 
                          layout=layout,
                          size=(1000, 800))
    
    savefig(histogram_plot, "mrna_histograms.png")
    println("Histograms saved to mrna_histograms.png")
    
    return histogram_plot
end

# Main function
function main()
    # Try to reset GKS state
    try
        Plots.reset_defaults()
        gr(fmt=:png)  # Use PNG format which tends to be more stable
    catch
        # Continue even if this fails
    end
    
    # Run the simulations and create animations
    sol_fast, sol_slow, results_fast, results_slow = run_simulations_and_animations()
    
    # Add histogram generation
    println("\nGenerating mRNA distribution histograms...")
    plot_mrna_histograms(sol_fast, sol_slow)
    
    println("\n--- Key Comparison Results ---")
    println("Fast Switching Fano Factor: $(round(results_fast["fano_factor"], digits=2))")
    println("Slow Switching Fano Factor: $(round(results_slow["fano_factor"], digits=2))")
    println("Noise Amplification Factor: $(round(results_slow["fano_factor"] / results_fast["fano_factor"], digits=2))x")
    
    # Print key concepts of comparison
    println("\n--- Key Insights ---")
    println("1. Fast switching: DNA state changes rapidly, mRNA shows less variability")
    println("2. Slow switching: DNA switches rarely, creating mRNA bursts and higher variability")
    println("3. Note how mRNA increases when DNA is ON and decreases when DNA is OFF")
    println("4. In fast switching, frequent ON/OFF cycles create more consistent mRNA levels")
    println("5. In slow switching, long ON periods create bursts of mRNA followed by decay during OFF periods")
    
    println("\nDone! Smooth animations and static images created.")
end

# Run the main function
main()