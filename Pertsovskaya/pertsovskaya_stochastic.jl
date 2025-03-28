using DifferentialEquations
using JumpProcesses
using Plots
using DataFrames
using CSV
using Random
using ProgressMeter
using Statistics
using Interpolations

# Define parameter values
p_values = Dict(
    :b_S => 0.0,           # Receptor activation rate
    :b_exp => 0.08,        # Nuclear export rate
    :b_imp => 0.013,       # Nuclear import rate 
    :b_ph => 1300.0,       # Phosphorylation rate
    :b_deph => 0.036,      # Dephosphorylation rate 
    :k_A => 4680.0,        # Scaling factor for STAT1
    :k_I => 82680.0,       # Scaling factor for SOCS1 inhibition
    :k_r => 23400.0,       # Scaling factor for SOCS1 transcription
    :k_f => 7.3e+03,       # Scaling factor for IRF1 transcription
    :k_F => 130e+03,       # Scaling factor for STAT1 transcription
    :b_A => 65.0,          # STAT1 translation rate
    :b_r => 12.8,          # SOCS1 transcription rate
    :b_R => 1.0e+02,       # SOCS1 translation rate
    :b_f => 2.7,           # IRF1 transcription rate
    :b_F => 10.0,          # IRF1 translation rate
    :b_a => 1e-01,         # STAT1 transcription rate
    :lambda_S => 0.02,     # Receptor degradation rate
    :lambda_r => 0.03,     # SOCS1 mRNA degradation rate
    :lambda_R => 0.02,     # SOCS1 protein degradation rate
    :lambda_f => 0.017,    # IRF1 mRNA degradation rate
    :lambda_F => 0.01,     # IRF1 protein degradation rate
    :lambda_a => 0.006,    # STAT1 mRNA degradation rate
    :q => 4,               # Hill coefficient for SOCS inhibition
    :n => 3,               # Hill coefficient for SOCS1 transcription
    :m => 2,               # Hill coefficient for IRF1 transcription
    :u => 1,               # Hill coefficient for STAT1 transcription
    :lambda_stat => 6.9e-04, # STAT1 degradation rate
    :basal_stat => 0.006    # Basal STAT1 transcription rate
)

# Define our custom safe rate functions
function phospho_rate(S, A, R, p)
    S_safe = max(S, 0.0)
    A_safe = max(A, 0.0)
    R_safe = max(R, 0.0)
    
    # Calculate inhibition term safely
    inhibition = 1.0 + (A_safe/p[:k_A]) + (R_safe/p[:k_I])^p[:q]
    
    return (p[:b_ph] * S_safe * (A_safe/p[:k_A])) / inhibition
end

function safe_hill(x, vmax, n, K)
    # Ensure x is non-negative 
    x_safe = max(x, 0.0)
    return vmax * (x_safe^n) / (K^n + x_safe^n)
end

function socs1_transcription(Ap_n, p)
    return safe_hill(Ap_n/p[:k_r], p[:b_r], p[:n], 1)
end

function irf1_transcription(Ap_n, p)
    return safe_hill(Ap_n/p[:k_f], p[:b_f], p[:m], 1)
end

function stat1_transcription(F, p)
    return safe_hill(F/p[:k_F], p[:b_a], p[:u], 1) + p[:basal_stat]
end

# Define species and their indices for the stochastic simulation
species_names = [:S, :A, :Ap_c, :Ap_n, :r, :R, :f, :F, :a]
S_idx, A_idx, Ap_c_idx, Ap_n_idx, r_idx, R_idx, f_idx, F_idx, a_idx = 1:9

# Initial conditions
u0 = [
    1000.0,        # S: Active receptors
    1e5,         # A: Non-phosphorylated STAT1
    100.0,         # Ap_c: Phosphorylated cytoplasmic STAT1
    1.0,          # Ap_n: Phosphorylated nuclear STAT1
    1.0,          # r: SOCS1 mRNA
    1.0,          # R: SOCS1 protein
    1.0,          # f: IRF1 mRNA
    1.0,          # F: IRF1 protein
    1.0           # a: STAT1 mRNA
]

# Define the time span (24 hours = 1440 minutes)
tspan = (0.0, 1440.0)

# Define the number of simulations for the ensemble
num_simulations = 1

# Create a discrete problem
function create_jump_problem(u0, tspan, p_values)
    # Create a discrete problem
    dprob = DiscreteProblem(u0, tspan, p_values)
    
    # Create all jump callbacks for each reaction
    
    # Reaction 1: Receptor production (0 -> S)
    rate1(u, p, t) = p[:b_S]
    affect1!(integrator) = (integrator.u[S_idx] += 1)
    jump1 = ConstantRateJump(rate1, affect1!)
    
    # Reaction 2: Receptor degradation (S -> 0)
    rate2(u, p, t) = p[:lambda_S] * u[S_idx]
    affect2!(integrator) = (integrator.u[S_idx] -= 1)
    jump2 = ConstantRateJump(rate2, affect2!)
    
    # Reaction 3: STAT1 phosphorylation (A -> Ap_c)
    rate3(u, p, t) = phospho_rate(u[S_idx], u[A_idx], u[R_idx], p)
    affect3!(integrator) = (integrator.u[A_idx] -= 1; integrator.u[Ap_c_idx] += 1)
    jump3 = ConstantRateJump(rate3, affect3!)
    
    # Reaction 4: STAT1 dephosphorylation (Ap_c -> A)
    rate4(u, p, t) = p[:b_deph] * u[Ap_c_idx]
    affect4!(integrator) = (integrator.u[Ap_c_idx] -= 1; integrator.u[A_idx] += 1)
    jump4 = ConstantRateJump(rate4, affect4!)
    
    # Reaction 5: Nuclear import (Ap_c -> Ap_n)
    rate5(u, p, t) = p[:b_imp] * u[Ap_c_idx]
    affect5!(integrator) = (integrator.u[Ap_c_idx] -= 1; integrator.u[Ap_n_idx] += 1)
    jump5 = ConstantRateJump(rate5, affect5!)
    
    # Reaction 6: Nuclear export (Ap_n -> Ap_c)
    rate6(u, p, t) = p[:b_exp] * u[Ap_n_idx]
    affect6!(integrator) = (integrator.u[Ap_n_idx] -= 1; integrator.u[Ap_c_idx] += 1)
    jump6 = ConstantRateJump(rate6, affect6!)
    
    # Reaction 7: SOCS1 mRNA transcription (0 -> r)
    rate7(u, p, t) = socs1_transcription(u[Ap_n_idx], p)
    affect7!(integrator) = (integrator.u[r_idx] += 1)
    jump7 = ConstantRateJump(rate7, affect7!)
    
    # Reaction 8: IRF1 mRNA transcription (0 -> f)
    rate8(u, p, t) = irf1_transcription(u[Ap_n_idx], p)
    affect8!(integrator) = (integrator.u[f_idx] += 1)
    jump8 = ConstantRateJump(rate8, affect8!)
    
    # Reaction 9: STAT1 mRNA transcription (0 -> a)
    rate9(u, p, t) = stat1_transcription(u[F_idx], p)
    affect9!(integrator) = (integrator.u[a_idx] += 1)
    jump9 = ConstantRateJump(rate9, affect9!)
    
    # Reaction 10: STAT1 translation (a -> a + A)
    rate10(u, p, t) = p[:b_A] * u[a_idx]
    affect10!(integrator) = (integrator.u[A_idx] += 1)
    jump10 = ConstantRateJump(rate10, affect10!)
    
    # Reaction 11: SOCS1 translation (r -> r + R)
    rate11(u, p, t) = p[:b_R] * u[r_idx]
    affect11!(integrator) = (integrator.u[R_idx] += 1)
    jump11 = ConstantRateJump(rate11, affect11!)
    
    # Reaction 12: IRF1 translation (f -> f + F)
    rate12(u, p, t) = p[:b_F] * u[f_idx]
    affect12!(integrator) = (integrator.u[F_idx] += 1)
    jump12 = ConstantRateJump(rate12, affect12!)
    
    # Reaction 13: SOCS1 mRNA degradation (r -> 0)
    rate13(u, p, t) = p[:lambda_r] * u[r_idx]
    affect13!(integrator) = (integrator.u[r_idx] -= 1)
    jump13 = ConstantRateJump(rate13, affect13!)
    
    # Reaction 14: SOCS1 protein degradation (R -> 0)
    rate14(u, p, t) = p[:lambda_R] * u[R_idx]
    affect14!(integrator) = (integrator.u[R_idx] -= 1)
    jump14 = ConstantRateJump(rate14, affect14!)
    
    # Reaction 15: IRF1 mRNA degradation (f -> 0)
    rate15(u, p, t) = p[:lambda_f] * u[f_idx]
    affect15!(integrator) = (integrator.u[f_idx] -= 1)
    jump15 = ConstantRateJump(rate15, affect15!)
    
    # Reaction 16: IRF1 protein degradation (F -> 0)
    rate16(u, p, t) = p[:lambda_F] * u[F_idx]
    affect16!(integrator) = (integrator.u[F_idx] -= 1)
    jump16 = ConstantRateJump(rate16, affect16!)
    
    # Reaction 17: STAT1 mRNA degradation (a -> 0)
    rate17(u, p, t) = p[:lambda_a] * u[a_idx]
    affect17!(integrator) = (integrator.u[a_idx] -= 1)
    jump17 = ConstantRateJump(rate17, affect17!)
    
    # Reaction 18: A degradation (A -> 0)
    rate18(u, p, t) = p[:lambda_stat] * u[A_idx]
    affect18!(integrator) = (integrator.u[A_idx] -= 1)
    jump18 = ConstantRateJump(rate18, affect18!)
    
    # Reaction 19: Ap_c degradation (Ap_c -> 0)
    rate19(u, p, t) = p[:lambda_stat] * u[Ap_c_idx]
    affect19!(integrator) = (integrator.u[Ap_c_idx] -= 1)
    jump19 = ConstantRateJump(rate19, affect19!)
    
    # Reaction 20: Ap_n degradation (Ap_n -> 0)
    rate20(u, p, t) = p[:lambda_stat] * u[Ap_n_idx]
    affect20!(integrator) = (integrator.u[Ap_n_idx] -= 1)
    jump20 = ConstantRateJump(rate20, affect20!)
    
    # Collect all jumps
    jumps = [jump1, jump2, jump3, jump4, jump5, jump6, jump7, jump8, jump9, jump10, 
             jump11, jump12, jump13, jump14, jump15, jump16, jump17, jump18, jump19, jump20]
    
    # Create jump problem with all jumps
    jump_prob = JumpProblem(dprob, Direct(), jumps...)
    
    return jump_prob
end

# Function to run multiple stochastic simulations
function run_stochastic_simulations(u0, tspan, p_values, num_simulations)
    # Create the jump problem 
    jump_prob = create_jump_problem(u0, tspan, p_values)
    
    # Run multiple simulations with different seeds
    solutions = []
    
    # Create a progress bar
    p = Progress(num_simulations, desc="Running simulations: ")
    
    for i in 1:num_simulations
        # Set random seed for reproducibility but different for each run
        Random.seed!(42 + i)
        
        # Solve the jump problem
        sol = solve(jump_prob, SSAStepper())
        
        # Store the solution
        push!(solutions, sol)
        
        # Update progress bar
        next!(p)
    end
    
    return solutions
end

# Function to interpolate solutions to common timepoints
function interpolate_to_common_timepoints(solutions, tspan; num_points=500)
    # Create common timepoints for interpolation
    common_times = range(tspan[1], tspan[2], length=num_points)
    
    # Initialize arrays to store interpolated data
    interp_data = zeros(length(solutions), length(species_names), num_points)
    
    # Loop through each solution
    for (sim_idx, sol) in enumerate(solutions)
        for species_idx in 1:length(species_names)
            # Create linear interpolation function
            interp_func = LinearInterpolation(sol.t, sol[species_idx, :], 
                                              extrapolation_bc=Flat())
            
            # Interpolate to common timepoints
            for (t_idx, t) in enumerate(common_times)
                interp_data[sim_idx, species_idx, t_idx] = interp_func(t)
            end
        end
    end
    
    return common_times, interp_data
end

# Function to calculate statistics from interpolated data
function calculate_statistics(interp_data)
    num_sims, num_species, num_timepoints = size(interp_data)
    
    # Initialize arrays for statistics
    means = zeros(num_species, num_timepoints)
    medians = zeros(num_species, num_timepoints)
    lower_quartiles = zeros(num_species, num_timepoints)
    upper_quartiles = zeros(num_species, num_timepoints)
    mins = zeros(num_species, num_timepoints)
    maxs = zeros(num_species, num_timepoints)
    stds = zeros(num_species, num_timepoints)
    
    # Calculate statistics for each species at each timepoint
    for species_idx in 1:num_species
        for t_idx in 1:num_timepoints
            # Extract all simulation values for this species at this timepoint
            values = interp_data[:, species_idx, t_idx]
            
            # Calculate statistics
            means[species_idx, t_idx] = mean(values)
            medians[species_idx, t_idx] = median(values)
            lower_quartiles[species_idx, t_idx] = quantile(values, 0.25)
            upper_quartiles[species_idx, t_idx] = quantile(values, 0.75)
            mins[species_idx, t_idx] = minimum(values)
            maxs[species_idx, t_idx] = maximum(values)
            stds[species_idx, t_idx] = std(values)
        end
    end
    
    return means, medians, lower_quartiles, upper_quartiles, mins, maxs, stds
end

# Function to save IRF1 trajectories to a dedicated file
function save_irf1_trajectories(common_times, interp_data)
    # Create DataFrame to store IRF1 data
    irf1_df = DataFrame(time_min = common_times, time_hours = common_times / 60.0)
    
    # Number of simulations
    num_sims = size(interp_data, 1)
    
    # Add IRF1 mRNA trajectories for each simulation
    for sim_idx in 1:num_sims
        irf1_df[!, "IRF1_mRNA_sim$(sim_idx)"] = interp_data[sim_idx, f_idx, :]
    end
    
    # Add IRF1 protein trajectories for each simulation
    for sim_idx in 1:num_sims
        irf1_df[!, "IRF1_protein_sim$(sim_idx)"] = interp_data[sim_idx, F_idx, :]
    end
    
    # Save DataFrame to CSV
    CSV.write(joinpath(@__DIR__, "irf1_trajectories.csv"), irf1_df)
end

# Function to plot the results with improved visualization
function plot_ensemble_results(common_times, interp_data, means, medians, lower_quartiles, upper_quartiles, stds)
    # Convert time from minutes to hours for plotting
    time_hours = common_times / 60.0
    
    # Plot individual time series for pSTAT1 cytoplasmic
    p1 = plot(title="pSTAT1 cytoplasmic - $num_simulations stochastic simulations",
              xlabel="Time (hours)", ylabel="Molecules", size=(800, 600), legend=:topright,
              dpi=300)
    
    # Plot a subset of individual trajectories with low opacity
    for i in 1:min(20, size(interp_data, 1))
        plot!(p1, time_hours, interp_data[i, Ap_c_idx, :], 
              label=i==1 ? "Individual sims" : "", linewidth=1, color=:lightblue, alpha=0.2)
    end
    
    # Plot statistics
    plot!(p1, time_hours, means[Ap_c_idx, :], linewidth=2, color=:blue, label="Mean")
    plot!(p1, time_hours, medians[Ap_c_idx, :], linewidth=2, color=:black, label="Median")
    
    # Plot confidence ribbon (25-75 percentile)
    plot!(p1, time_hours, lower_quartiles[Ap_c_idx, :], fillrange=upper_quartiles[Ap_c_idx, :],
          fillalpha=0.3, color=:blue, label="25-75th percentiles")
    
    # Create a multi-panel plot for other species
    species_to_plot = [
        (S_idx, "Receptor (S)"),
        ((A_idx, Ap_c_idx, Ap_n_idx), "Total STAT1 (A + Ap_c + Ap_n)"),
        (r_idx, "SOCS1 mRNA"),
        (R_idx, "SOCS1 protein (R)"),
        (f_idx, "IRF1 mRNA"),
        (F_idx, "IRF1 protein (F)"),
        (a_idx, "STAT1 mRNA (a)")  # Add STAT1 mRNA to the plot panels
    ]
    
    p_combined = plot(layout=(length(species_to_plot), 1), size=(800, 200*length(species_to_plot)), dpi=300)
    
    # For each species or combination, plot statistics
    for (i, (idx_or_idxs, title)) in enumerate(species_to_plot)
        if isa(idx_or_idxs, Int)
            # Single species
            plot!(p_combined[i], time_hours, means[idx_or_idxs, :], linewidth=2, color=:blue, label="Mean")
            plot!(p_combined[i], time_hours, medians[idx_or_idxs, :], linewidth=2, color=:black, label="Median")
            plot!(p_combined[i], time_hours, lower_quartiles[idx_or_idxs, :], fillrange=upper_quartiles[idx_or_idxs, :],
                  fillalpha=0.3, color=:blue, label="25-75th percentiles", title=title)
            
            # For STAT1 mRNA, additionally plot individual trajectories to show stochasticity
            if idx_or_idxs == a_idx
                for j in 1:min(5, size(interp_data, 1))  # Plot 5 individual mRNA trajectories
                    plot!(p_combined[i], time_hours, interp_data[j, a_idx, :], 
                          linewidth=0.5, color=:lightblue, alpha=0.3, label=j==1 ? "Individual sims" : "")
                end
            end
        else
            # Combination of species (sum)
            combined_means = sum(means[idx, :] for idx in idx_or_idxs)
            combined_medians = sum(medians[idx, :] for idx in idx_or_idxs)
            combined_lower = sum(lower_quartiles[idx, :] for idx in idx_or_idxs)
            combined_upper = sum(upper_quartiles[idx, :] for idx in idx_or_idxs)
            
            # Plot statistics
            plot!(p_combined[i], time_hours, combined_means, linewidth=2, color=:blue, label="Mean")
            plot!(p_combined[i], time_hours, combined_medians, linewidth=2, color=:black, label="Median")
            plot!(p_combined[i], time_hours, combined_lower, fillrange=combined_upper,
                  fillalpha=0.3, color=:blue, label="25-75th percentiles", title=title)
        end
        
        # For the last subplot, add x-label
        if i == length(species_to_plot)
            plot!(p_combined[i], xlabel="Time (hours)")
        end
        
        # Add y-label to first subplot only
        if i == 1
            plot!(p_combined[i], ylabel="Molecules")
        end
        
        # Add legend to first subplot only, others can use same colors
        if i == 1
            plot!(p_combined[i], legend=:topright)
        else
            plot!(p_combined[i], legend=false)
        end
    end
    
    # Create an additional panel specifically for mRNA species
    p_mrna = plot(layout=(3,1), size=(800, 600), dpi=300)
    
    # Plot all mRNAs: STAT1 mRNA, SOCS1 mRNA, IRF1 mRNA
    mrna_idx = [a_idx, r_idx, f_idx]
    mrna_names = ["STAT1 mRNA (a)", "SOCS1 mRNA (r)", "IRF1 mRNA (f)"]
    
    for i in 1:3
        plot!(p_mrna[i], title=mrna_names[i])
        
        # Plot mean
        plot!(p_mrna[i], time_hours, means[mrna_idx[i], :], linewidth=2, color=:blue, label="Mean")
        
        # Plot median
        plot!(p_mrna[i], time_hours, medians[mrna_idx[i], :], linewidth=2, color=:black, label="Median")
        
        # Plot quartile range
        plot!(p_mrna[i], time_hours, lower_quartiles[mrna_idx[i], :], 
              fillrange=upper_quartiles[mrna_idx[i], :],
              fillalpha=0.3, color=:blue, label="25-75th percentiles")
        
        # Plot a few individual trajectories
        for j in 1:min(5, size(interp_data, 1))
            plot!(p_mrna[i], time_hours, interp_data[j, mrna_idx[i], :], 
                  linewidth=0.5, color=:lightblue, alpha=0.2, label=j==1 ? "Individual sims" : "")
        end
        
        if i == 3
            plot!(p_mrna[i], xlabel="Time (hours)")
        end
        
        if i == 1
            plot!(p_mrna[i], legend=:topright)
        else
            plot!(p_mrna[i], legend=false)
        end
        
        plot!(p_mrna[i], ylabel="Molecules")
    end
    
    return p1, p_combined, p_mrna
end

# Run the stochastic simulations
println("Running $num_simulations stochastic JAK-STAT model simulations...")
solutions = run_stochastic_simulations(u0, tspan, p_values, num_simulations)
println("Stochastic simulations completed successfully!")

# Interpolate solutions to common timepoints
println("Interpolating solutions to common timepoints...")
common_times, interp_data = interpolate_to_common_timepoints(solutions, tspan)
println("Interpolation completed!")

# Calculate statistics
println("Calculating statistics across $num_simulations simulations...")
means, medians, lower_quartiles, upper_quartiles, mins, maxs, stds = calculate_statistics(interp_data)
println("Statistics calculations completed!")

# Save IRF1 trajectories to a dedicated file
println("Saving IRF1 trajectories to file...")
save_irf1_trajectories(common_times, interp_data)
println("IRF1 trajectories saved successfully!")

# Plot the results
println("Generating plots...")
oscillation_plot, combined_plot, mrna_plot = plot_ensemble_results(common_times, interp_data, means, 
                                                    medians, lower_quartiles, upper_quartiles, stds)

# Display plots
display(oscillation_plot)
display(combined_plot)
display(mrna_plot)

# Save plots
savefig(oscillation_plot, joinpath(@__DIR__, "pstat1_oscillation_ensemble100.png"))
savefig(combined_plot, joinpath(@__DIR__, "ifn_stat_dynamics_ensemble100.png"))
savefig(mrna_plot, joinpath(@__DIR__, "mrna_dynamics_ensemble100.png"))


# Create a pathway diagram using ASCII art
pathway_diagram = """
JAK-STAT Pathway Stochastic Model (Pertsovskaya et al.)

External Signal
     |
     ↓
    [S] → Receptor activation
     |
     ↓
[A] → [Ap_c] → [Ap_n] → Gene activation
 ↑      ↑        |
 |      |        ↓
 |      ↓        |
[a] ← [F] ← [f] ←+ 
 |      ↑  
 |      |
 +→ [R] ← [r] ←+
     |         |
     +→ Inhibition

Note: This is an ensemble of $num_simulations stochastic simulations using the Gillespie algorithm
"""

println(pathway_diagram)