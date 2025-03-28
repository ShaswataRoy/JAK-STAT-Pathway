using Catalyst
using OrdinaryDiffEq
using Plots
using DataFrames
using CSV

# Define our custom hill function to avoid numerical issues
function safe_hill(x, vmax, n, K)
    # Ensure x is non-negative 
    x_safe = max(x, 0.0)
    return vmax * (x_safe^n) / (K^n + x_safe^n)
end

# Define our custom phosphorylation rate function
function phospho_rate(S, A, R, b_ph, k_A, k_I, q)
    S_safe = max(S, 0.0)
    A_safe = max(A, 0.0)
    R_safe = max(R, 0.0)
    
    # Calculate inhibition term safely
    inhibition = 1.0 + (A_safe/k_A) + (R_safe/k_I)^q
    
    return (b_ph * S_safe * (A_safe/k_A)) / inhibition
end

# Define the JAK-STAT model directly as an ODE system
function pertsovskaya_ode!(du, u, p, t)
    # Extract state variables
    S = u[1]      # Active receptors
    A = u[2]      # Non-phosphorylated STAT1
    Ap_c = u[3]   # Phosphorylated STAT1 (cytoplasm)
    Ap_n = u[4]   # Phosphorylated STAT1 (nucleus)
    r = u[5]      # SOCS1 mRNA
    R = u[6]      # SOCS1 protein
    f = u[7]      # IRF1 mRNA
    F = u[8]      # IRF1 protein
    a = u[9]      # STAT1 mRNA
    
    # Extract parameters
    b_S = p[1]
    b_exp = p[2]
    b_imp = p[3]
    b_ph = p[4]
    b_deph = p[5]
    k_A = p[6]
    k_I = p[7]
    k_r = p[8]
    k_f = p[9]
    k_F = p[10]
    b_A = p[11]
    b_r = p[12]
    b_R = p[13]
    b_f = p[14]
    b_F = p[15]
    b_a = p[16]
    lambda_S = p[17]
    lambda_r = p[18]
    lambda_R = p[19]
    lambda_f = p[20]
    lambda_F = p[21]
    lambda_a = p[22]
    q = p[23]
    n = p[24]
    m = p[25]
    u_param = p[26]  # renamed to avoid conflict with u
    lambda_stat = p[27]
    
    # Calculate rates using our safe custom functions
    phosphorylation = phospho_rate(S, A, R, b_ph, k_A, k_I, q)
    SOCS1_transcription = safe_hill(Ap_n/k_r, b_r, n, 1)
    IRF1_transcription = safe_hill(Ap_n/k_f, b_f, m, 1)
    STAT1_transcription = safe_hill(F/k_F, b_a, u_param, 1)
    
    # S: Active receptor dynamics
    du[1] = b_S - lambda_S * S
    
    # A: STAT1 non-phosphorylated form
    du[2] = b_exp * Ap_n + b_deph * Ap_c + b_A * a - phosphorylation - lambda_stat * A
    
    # Ap_c: Phosphorylated STAT1 (cytoplasm)
    du[3] = phosphorylation - b_imp * Ap_c - b_deph * Ap_c - lambda_stat * Ap_c
    
    # Ap_n: Phosphorylated STAT1 (nucleus)
    du[4] = b_imp * Ap_c - b_exp * Ap_n - lambda_stat * Ap_n
    
    # r: SOCS1 mRNA
    du[5] = SOCS1_transcription - lambda_r * r
    
    # R: SOCS1 protein
    du[6] = b_R * r - lambda_R * R
    
    # f: IRF1 mRNA
    du[7] = IRF1_transcription - lambda_f * f
    
    # F: IRF1 protein
    du[8] = b_F * f - lambda_F * F
    
    # a: STAT1 mRNA
    du[9] = STAT1_transcription  - lambda_a * a
    
    # Ensure non-negativity for all species
    for i in 1:length(du)
        if u[i] <= 0 && du[i] < 0
            du[i] = 0.0
        end
    end
end

# Define parameter values from the MATLAB file
p_values = [
    0.0,                # b_S
    0.08, 0.013,        # b_exp, b_imp 
    1300, 0.036,        # b_ph, b_deph 
    4680, 82680,        # k_A, k_I 
    23400, 7.3e+03,     # k_r, k_f
    # 130e+03,            # k_F  
    1.3e+03,            # k_F 
    65,                 # b_A
    12.8, 1.0e+02,      # b_r, b_R
    2.7, 10,            # b_f, b_Fl
    1e-01,              # b_a   
    0.02,               # lambda_S
    0.03, 0.02,         # lambda_r, lambda_R  
    0.017, 0.01,        # lambda_f, lambda_F  
    5.8e-4,              # lambda_a  
    4, 3, 2, 1,         # q, n, m, u
    6.9e-04  # lambda_stat
]

# Initial conditions from the MATLAB file
u0 = [
    500.0,        # S: Active receptors
    1e+05,         # A: Non-phosphorylated STAT1
    10.0,          # Ap_c: Phosphorylated cytoplasmic STAT1
    1.0,           # Ap_n: Phosphorylated nuclear STAT1
    1.0,           # r: SOCS1 mRNA
    1.0,           # R: SOCS1 protein
    1.0,           # f: IRF1 mRNA
    1.0,           # F: IRF1 protein
    0.0            # a: STAT1 mRNA
]

# Define the time span (24 hours = 1440 minutes)
tspan = (0.0, 24*60.0)

# Create an ODE problem
ode_problem = ODEProblem(pertsovskaya_ode!, u0, tspan, p_values)

# Solve using a stiff solver
sol = solve(ode_problem, Rodas5(), abstol=1e-8, reltol=1e-8, saveat=1.0)

# Function to create separate plots for each variable
function plot_results(sol)
    # Convert time from minutes to hours for plotting
    time_hours = sol.t / 60.0
    
    # Variable names for titles
    var_names = ["Active Receptor", "STAT1 (Non-phospho)", "pSTAT1 (Cytoplasmic)", 
                 "pSTAT1 (Nuclear)", "SOCS1 mRNA", "SOCS1 Protein", 
                 "IRF1 mRNA", "IRF1 Protein", "STAT1 mRNA"]
    
    # Define colors for each plot
    colors = [:black, :blue, :red, :purple, :green, :orange, :darkblue, :brown, :magenta]
    
    # Create individual plots for each variable
    plots = []
    
    # Font settings for consistency
    font_settings = (fontfamily="sans-serif", titlefontsize=10, guidefontsize=8)
    
    for i in 1:9
        p = plot(time_hours, sol[i, :], 
                 label="", 
                 title=var_names[i],
                 xlabel="Time (hours)",
                 ylabel="Concentration",
                 linewidth=2,
                 color=colors[i];
                 font_settings...)
        push!(plots, p)
    end
    
    # Add plots for derived variables
    
    # Total pSTAT1
    p_total_pstat = plot(time_hours, sol[3, :] + sol[4, :], 
                         label="", 
                         title="Total Phosphorylated STAT1",
                         xlabel="Time (hours)",
                         ylabel="Concentration",
                         linewidth=2,
                         color=:cyan;
                         font_settings...)
    push!(plots, p_total_pstat)
    
    # Total STAT1
    p_total_stat = plot(time_hours, sol[2, :] + sol[3, :] + sol[4, :], 
                        label="", 
                        title="Total STAT1",
                        xlabel="Time (hours)",
                        ylabel="Concentration",
                        linewidth=2,
                        color=:darkgreen;
                        font_settings...)
    push!(plots, p_total_stat)
    
    # Create a grid layout of all plots
    grid_plot = plot(plots..., layout=(4, 3), size=(900, 800), 
                    plot_title="JAK-STAT Pathway Components")
    
    return grid_plot
end

# Generate and display the plots
grid_plot = plot_results(sol)
display(grid_plot)

# Save the plot
savefig(grid_plot, joinpath(@__DIR__, "all_variables_plot.png"))

# Create and save a CSV file with the simulation data
function save_simulation_data(sol)
    df = DataFrame(
        time_min = sol.t,
        time_hours = sol.t / 60.0,  # Add time in hours
        S = sol[1, :],       # Active receptors
        A = sol[2, :],       # Non-phosphorylated STAT1
        Ap_c = sol[3, :],    # Phosphorylated cytoplasmic STAT1
        Ap_n = sol[4, :],    # Phosphorylated nuclear STAT1
        r = sol[5, :],       # SOCS1 mRNA
        R = sol[6, :],       # SOCS1 protein
        f = sol[7, :],       # IRF1 mRNA
        F = sol[8, :],       # IRF1 protein
        a = sol[9, :],       # STAT1 mRNA
        total_STAT1 = sol[2, :] + sol[3, :] + sol[4, :],  # Total STAT1
        total_pSTAT1 = sol[3, :] + sol[4, :]              # Total phosphorylated STAT1
    )
    
    CSV.write(joinpath(@__DIR__, "ifn_stat_simulation_data.csv"), df)
    println("Simulation data saved to ifn_stat_simulation_data.csv")
    
    return df
end

# Save the simulation data
# simulation_data = save_simulation_data(sol)

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
println("Model simulation completed successfully!")
println("Simulation time span: 0 to $(tspan[2]/60.0) hours ($(tspan[2]) minutes)")
println("Number of species: 9")
println("Maximum phosphorylated STAT1 in cytoplasm: ", maximum(sol[3, :]))
println("Maximum phosphorylated STAT1 in nucleus: ", maximum(sol[4, :]))
println("Final SOCS1 protein level: ", sol[6, end])
println("Final IRF1 protein level: ", sol[8, end])