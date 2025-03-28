using CSV
using DataFrames
using DifferentialEquations
using JumpProcesses
using Plots
using Interpolations
using Statistics
using Random
using StatsPlots
using ProgressMeter

# 1. Load the IRF1 protein trajectory
function load_irf1_protein_trajectory(filepath)
    data = CSV.read(filepath, DataFrame)
    irf1_protein = data[:, "IRF1_protein_sim1"]
    time_min = data[:, "time_min"]
    return irf1_protein, time_min
end

# Path to the CSV file
filepath = joinpath(@__DIR__, "irf1_trajectories.csv")

# Load the IRF1 protein trajectory
irf1_protein, time_min = load_irf1_protein_trajectory(filepath)

# 2. Create an interpolation function for the protein data
protein_interpolation = linear_interpolation(time_min, irf1_protein)

# 3. Define the stochastic model parameters
# Parameter values (adjust these as needed for your model)
p_values = Dict(
    :V_max => .3,       # V_max
    :K_m => 1.3e3,       # K_m (chosen relative to the protein level range)
    :k_deg => 5*6.9e-04,   # k_deg
    :basal_rate => 0.0   # basal_rate
)

# 4. Set up the stochastic simulation
tspan = (0.0, time_min[end])  # Simulate over the same time range as our data

# Initial condition: starting with no mRNA
u0 = [0.0]  # mRNA count as float for ODE solver

# Create time in hours for plotting
time_hrs = time_min ./ 60.0

# Number of stochastic ensemble simulations to run
num_simulations = 1000  # Increased for better histogram distribution

# Create a function to build the SDE problem
function create_sde_problem(u0, tspan, p_values)
    # Define the deterministic part (mRNA dynamics) - using in-place convention
    function mrna_dynamics!(du, u, p, t)
        # Get protein level at current time
        prot_level = t <= time_min[end] ? protein_interpolation(t) : protein_interpolation(time_min[end])
        
        # Michaelis-Menten kinetics for synthesis
        synthesis_rate = p[:basal_rate] + p[:V_max] * prot_level / (p[:K_m] + prot_level)
        
        # First-order degradation
        degradation_rate = p[:k_deg] * u[1]
        
        # Rate of change
        du[1] = synthesis_rate - degradation_rate
    end
    
    # Define the noise function - using in-place convention to match mrna_dynamics!
    function noise_func!(du, u, p, t)
        # Noise proportional to square root of state (Poisson-like noise)
        du[1] = sqrt(max(0.1, u[1]))
    end
    
    # Create SDE problem
    sde_prob = SDEProblem(mrna_dynamics!, noise_func!, u0, tspan, p_values)
    
    return sde_prob
end

# 5. Run the ensemble of simulations
println("Running ensemble of $num_simulations stochastic simulations...")

# Run the ensemble simulations
ensemble_solutions = []

# Create a progress meter
p = Progress(num_simulations, desc="Running simulations: ", dt=0.5)

for i in 1:num_simulations
    # Set different random seed for each simulation
    Random.seed!(42 + i)
    
    # Create a new SDE problem for each simulation
    sde_prob = create_sde_problem(u0, tspan, p_values)
    
    # Solve the stochastic model with regular sampling
    sol = solve(sde_prob, SRIW1(), saveat=time_min)
    
    # Store the solution
    push!(ensemble_solutions, sol)
    
    # Update progress meter
    next!(p)
end

println("Ensemble of $num_simulations stochastic simulations completed!")

# 6. Calculate ensemble statistics
# Extract mRNA trajectories from all simulations for easier calculation
mRNA_trajectories = zeros(num_simulations, length(time_min))
for (i, sol) in enumerate(ensemble_solutions)
    # Ensure no negative values in simulation
    mRNA_trajectories[i, :] = max.(0, sol[1, :])
end

# Calculate statistics
mean_mRNA = mean(mRNA_trajectories, dims=1)[1, :]
std_mRNA = std(mRNA_trajectories, dims=1)[1, :]

# 7. Plot the results
plt = plot(layout=(2,1), size=(800, 600))

# Top plot: Protein input
plot!(plt[1], time_hrs, irf1_protein, 
      label="IRF1 Protein (Input)", 
      title="IRF1 Protein Level",
      xlabel="",
      ylabel="Protein Level",
      linewidth=2,
      color=:blue,
      legend=:topleft)

# Bottom plot: mRNA stochastic simulations
plot!(plt[2], time_hrs, mean_mRNA, 
      label="mRNA (Mean)", 
      title="mRNA Dynamics (Stochastic Model, n=$num_simulations)",
      xlabel="Time (hours)",
      ylabel="mRNA Count",
      linewidth=2,
      color=:red,
      legend=:topleft)

# Save the figures
savefig(plt, joinpath(@__DIR__,"mrna/irf1_stochastic_simulation.png"))

# 8. Create a second plot showing just the individual trajectories
plt_trajectories = plot(time_hrs, irf1_protein,
      label="IRF1 Protein (Input)",
      title="IRF1 mRNA Stochastic Trajectories (n=$num_simulations)",
      xlabel="Time (hours)",
      ylabel="Molecule Count",
      linewidth=2,
      color=:blue,
      legend=:topleft)

# Plot each trajectory with a different color
for i in 1:min(10, num_simulations)  # Only plot first 10 trajectories to avoid overcrowding
    plot!(plt_trajectories, time_hrs, mRNA_trajectories[i, :],
          linewidth=1.5,
          alpha=0.7,
          label="mRNA Trajectory $i")
end

savefig(plt_trajectories, joinpath(@__DIR__,"mrna/irf1_stochastic_trajectories.png"))

# 9. NEW: Create histograms of mRNA counts at 6 different time points
# Select 6 time points distributed across the simulation
time_indices = round.(Int, range(1, length(time_min), length=6))
selected_times = time_min[time_indices]

# Create a 2x3 plot layout for histograms
plt_histograms = plot(layout=(2, 3), size=(1200, 800), 
                      legend=false, 
                      title="mRNA Count Distribution at Different Time Points")

# Color palette for the histograms
colors = [:blue, :green, :red, :purple, :orange, :brown]

# Create histograms for each selected time point
for (i, time_idx) in enumerate(time_indices)
    # Extract mRNA counts at this time point from all simulations
    counts_at_time = mRNA_trajectories[:, time_idx]
    
    # Get time in hours for nicer display
    time_hrs = round(selected_times[i] / 60, digits=1)
    
    # Determine appropriate bin count based on data range
    bin_count = max(10, min(30, round(Int, maximum(counts_at_time) - minimum(counts_at_time) + 1)))
    
    # Create histogram
    histogram!(plt_histograms[i], counts_at_time, 
            #    bins=bin_count,
               bins = 0:1:100,  # Set bins to 0-100
               color=colors[i],
               alpha=0.7,
               title="t = $(time_hrs) hrs",
               xlabel="mRNA Count",
               ylabel="Frequency",
               xlims = (0, 100),
               normalize=:pdf)  # Use probability density for better comparison
end

# Adjust overall title
plot!(plt_histograms, plot_title="mRNA Count Distribution (n=$num_simulations simulations)")

# Save the histogram plot
savefig(plt_histograms, joinpath(@__DIR__,"mrna/mrna_histograms.png"))

println("Stochastic simulation complete! Results saved.")
println("Parameters used:")
println("V_max: $(p_values[:V_max]), K_m: $(p_values[:K_m]), k_deg: $(p_values[:k_deg]), basal_rate: $(p_values[:basal_rate])")
println("Number of simulations: $num_simulations")
println("Histograms created for time points (hrs): $(round.(selected_times ./ 60, digits=1))")