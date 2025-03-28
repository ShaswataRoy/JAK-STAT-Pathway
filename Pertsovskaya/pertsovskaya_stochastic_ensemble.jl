using DifferentialEquations
using JumpProcesses
using Plots
using DataFrames
using CSV
using Random
using ProgressMeter
using Statistics
using Interpolations

# Define parameter values as a dictionary for easier access in jump system
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
    :k_F => 1.3e+03,       # Scaling factor for STAT1 transcription
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
    100.0,        # S: Active receptors
    1e5,          # A: Non-phosphorylated STAT1
    10.0,         # Ap_c: Phosphorylated cytoplasmic STAT1
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
num_simulations = 10

# Create a discrete problem and jump system for stochastic simulations
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

# Run the ensemble of stochastic simulations
println("Running ensemble of $num_simulations stochastic JAK-STAT simulations...")

# Create jump problem
jump_prob = create_jump_problem(u0, tspan, p_values)

# Run the ensemble simulations
ensemble_solutions = []
p = Progress(num_simulations, 1, "Running simulations: ")
for i in 1:num_simulations
    # Set different random seed for each simulation
    Random.seed!(42 + i)
    
    # Solve the stochastic model with regular sampling for easier visualization
    sol = solve(jump_prob, SSAStepper(), saveat=10.0)
    
    # Store the solution
    push!(ensemble_solutions, sol)
    
    next!(p)
end

println("Ensemble of $num_simulations stochastic simulations completed successfully!")

# Function to interpolate solutions to common timepoints for consistent visualization
function interpolate_solutions(ensemble_solutions, tspan; num_points=500)
    # Generate common time points
    common_times = range(tspan[1], tspan[2], length=num_points)
    
    # Initialize array to store interpolated data
    interp_data = zeros(length(ensemble_solutions), length(species_names), num_points)
    
    # Interpolate each solution to common time points
    for (sim_idx, sol) in enumerate(ensemble_solutions)
        for species_idx in 1:length(species_names)
            # Create linear interpolation function
            if length(sol.t) > 1  # Ensure we have at least 2 points for interpolation
                interp_func = LinearInterpolation(sol.t, sol[species_idx, :], extrapolation_bc=Flat())
                
                # Interpolate to common timepoints
                for (t_idx, t) in enumerate(common_times)
                    interp_data[sim_idx, species_idx, t_idx] = interp_func(t)
                end
            else
                # Handle the case of only one point (unlikely but possible)
                for (t_idx, _) in enumerate(common_times)
                    interp_data[sim_idx, species_idx, t_idx] = sol[species_idx, 1]
                end
            end
        end
    end
    
    return common_times, interp_data
end

# Interpolate solutions to common timepoints
println("Interpolating solutions to common timepoints...")
common_times, interp_data = interpolate_solutions(ensemble_solutions, tspan)
println("Interpolation completed!")

# Save all individual trajectories rather than just statistics
function save_all_trajectories(common_times, interp_data)
    num_sims, num_species, num_timepoints = size(interp_data)
    
    # Create one DataFrame per simulation
    for sim_idx in 1:num_sims
        # Initialize DataFrame with timepoints
        df = DataFrame(time_min = common_times, time_hours = common_times / 60.0)
        
        # Add all species for this simulation
        for (species_idx, species_name) in enumerate(species_names)
            df[!, string(species_name)] = interp_data[sim_idx, species_idx, :]
        end
        
        # Save to CSV
        CSV.write(joinpath(@__DIR__, "trajectory_$(sim_idx).csv"), df)
    end
    
    println("All $num_sims individual trajectories saved to CSV files")
    
    # Also create a summary DataFrame with means across runs
    summary_df = DataFrame(time_min = common_times, time_hours = common_times / 60.0)
    
    # Calculate means for each species
    for (species_idx, species_name) in enumerate(species_names)
        summary_df[!, "$(species_name)_mean"] = mean(interp_data[:, species_idx, :], dims=1)[1, :]
    end
    
    # Save summary
    CSV.write(joinpath(@__DIR__, "trajectories_summary.csv"), summary_df)
    
    return summary_df
end

# Save all trajectories to CSV files
println("Saving all individual trajectories...")
summary_df = save_all_trajectories(common_times, interp_data)
println("All trajectories saved!")

# Function to plot all individual trajectories and overlay them
function plot_all_trajectories(common_times, interp_data)
    # Convert time from minutes to hours
    time_hours = common_times / 60.0
    
    # Create plots for key species
    key_species = [(Ap_c_idx, "pSTAT1 cytoplasmic"),
                  (Ap_n_idx, "pSTAT1 nuclear"),
                  (a_idx, "STAT1 mRNA"),
                  (r_idx, "SOCS1 mRNA"),
                  (R_idx, "SOCS1 protein"),
                  (f_idx, "IRF1 mRNA"),
                  (F_idx, "IRF1 protein")]
    
    plots = []
    
    # Generate a plot for each key species with all trajectories
    for (idx, title) in key_species
        p = plot(title="$title - All $num_simulations trajectories",
                xlabel="Time (hours)", ylabel="Molecules",
                legend=false, size=(800, 500), dpi=300)
        
        # Plot each trajectory with a different color
        for sim_idx in 1:size(interp_data, 1)
            plot!(p, time_hours, interp_data[sim_idx, idx, :], 
                  linewidth=1.5, alpha=0.7)
        end
        
        push!(plots, p)
    end
    
    # Create a combined figure for mRNAs
    p_mrna = plot(layout=(3,1), size=(800, 900), dpi=300)
    
    mrna_idx = [a_idx, r_idx, f_idx]
    mrna_titles = ["STAT1 mRNA", "SOCS1 mRNA", "IRF1 mRNA"]
    
    for i in 1:3
        plot!(p_mrna[i], title=mrna_titles[i], ylabel="Molecules")
        
        # Plot each trajectory with a different color
        for sim_idx in 1:size(interp_data, 1)
            plot!(p_mrna[i], time_hours, interp_data[sim_idx, mrna_idx[i], :], 
                 linewidth=1.5, alpha=0.7)
        end
        
        if i == 3
            plot!(p_mrna[i], xlabel="Time (hours)")
        end
    end
    
    # Create a combined figure for proteins
    p_protein = plot(layout=(3,1), size=(800, 900), dpi=300)
    
    protein_idx = [Ap_c_idx, R_idx, F_idx]
    protein_titles = ["pSTAT1 cytoplasmic", "SOCS1 protein", "IRF1 protein"]
    
    for i in 1:3
        plot!(p_protein[i], title=protein_titles[i], ylabel="Molecules")
        
        # Plot each trajectory with a different color
        for sim_idx in 1:size(interp_data, 1)
            plot!(p_protein[i], time_hours, interp_data[sim_idx, protein_idx[i], :], 
                 linewidth=1.5, alpha=0.7)
        end
        
        if i == 3
            plot!(p_protein[i], xlabel="Time (hours)")
        end
    end
    
    return plots, p_mrna, p_protein
end

# Plot all trajectories
println("Generating plots of all individual trajectories...")
individual_plots, mrna_plot, protein_plot = plot_all_trajectories(common_times, interp_data)

# Display the mRNA and protein plots
display(mrna_plot)
display(protein_plot)

# Save plots
savefig(mrna_plot, joinpath(@__DIR__, "mrna_trajectories_all.png"))
savefig(protein_plot, joinpath(@__DIR__, "protein_trajectories_all.png"))

# Save individual species plots
for (i, (idx, title)) in enumerate([(Ap_c_idx, "pSTAT1_cytoplasmic"),
                                    (Ap_n_idx, "pSTAT1_nuclear"),
                                    (a_idx, "STAT1_mRNA"),
                                    (r_idx, "SOCS1_mRNA"),
                                    (R_idx, "SOCS1_protein"),
                                    (f_idx, "IRF1_mRNA"),
                                    (F_idx, "IRF1_protein")])
    savefig(individual_plots[i], joinpath(@__DIR__, "$(title)_trajectories.png"))
end

# Create a pathway diagram using ASCII art
pathway_diagram = """
JAK-STAT Pathway Model (Pertsovskaya et al.)

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
"""

println(pathway_diagram)

# Print model summary
println("Stochastic ensemble simulation completed successfully!")
println("Number of simulations: $num_simulations")
println("Simulation time span: 0 to $(tspan[2]/60.0) hours ($(tspan[2]) minutes)")
println("All trajectories have been saved individually as CSV files")
println("All trajectories have been plotted and saved as PNG files")