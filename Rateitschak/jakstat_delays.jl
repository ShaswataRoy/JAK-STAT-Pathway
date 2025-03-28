using Catalyst
using DifferentialEquations
using BenchmarkTools
using Plots
using DelayDiffEq
using Distributions
using QuadGK
using ProgressMeter
theme(:mute);

# Parameters - keeping the same constants
const k0 = 0.0039
const k1 = 0.03
const k2 = 0.05
const k3 = 2.32
const k4 = 31.14
const k5 = 0.32
const k6 = 2.82
const k7 = 5.49
const k8 = 4.37
const k9 = 18.75
const k10 = 0.0083
const k11 = 0.24
const k12 = 0.3
const k13 = 1.33
const k14 = 0.0048
const k15 = 0.019
const k16 = 0.025
const k17 = 0.022

# Average delay times (in minutes)
const tau_avg_1 = 39  # STAT1Dn → SOCS1
const tau_avg_2 = 22  # STAT1Dn → IRF1 mRNA
const tau_avg_3 = 164 # SOCS1 → Inhibition of STAT1c phosphorylation
const tau_avg_4 = 185 # IFNRJ → Total STAT1

# Gamma distribution parameters
# With shape = 4, scale = tau_avg/4
const shape = 4
const scale_1 = tau_avg_1/4
const scale_2 = tau_avg_2/4
const scale_3 = tau_avg_3/4
const scale_4 = tau_avg_4/4

# Distribution objects
const gamma_dist_1 = Gamma(shape, scale_1)
const gamma_dist_2 = Gamma(shape, scale_2)
const gamma_dist_3 = Gamma(shape, scale_3)
const gamma_dist_4 = Gamma(shape, scale_4)

# Total receptor concentration
const I_tot = 0.0038

# Simulation time settings
const T = 10*60.0  # 24 hours in minutes
const dt = 0.1
tspan = (0.0, T)

# Store history values in a global dictionary for efficient lookup
const history_buffer = Dict{Float64, Vector{Float64}}()

# Define the gamma kernel functions
gamma_kernel_1(τ) = pdf(gamma_dist_1, τ)
gamma_kernel_2(τ) = pdf(gamma_dist_2, τ)
gamma_kernel_3(τ) = pdf(gamma_dist_3, τ)
gamma_kernel_4(τ) = pdf(gamma_dist_4, τ)

# Define the history function for the DDEs
function history_function(p, t)
    if t < 0
        # Return initial conditions for t ≤ 0
        return u0
    else
        return history_buffer[t]
    end
end

# Helper function to save state at each time point
function store_history!(integrator)
    t = integrator.t
    history_buffer[t] = copy(integrator.u)
    return nothing
end

# Helper function to compute distributed delay response
function compute_distributed_delay(kernel, species_idx, h, p, t)
    if t <0
        return 0.0
    end
    
    # Compute the convolution integral numerically
    integral, _ = quadgk(τ -> kernel(τ) * h(p, t - τ)[species_idx], 0, t)
    
    return integral
end

# Define the delay differential equation system with distributed delays
function dde_model!(du, u, h, p, t)
    # u: current state
    # h: history function to access past values h(p, t-tau)
    
    # Current state variables
    IFNg = u[1]         # IFN gamma
    IFNRJ = u[2]        # IFNRJ - complex IFNγ and receptor
    STAT1c = u[3]       # Cytoplasmic STAT1 Dimer
    STAT1n = u[4]       # Free Nuclear STAT1 Dimer
    STAT1Dn = u[5]      # Nuclear STAT1 Dimer bound to DNA
    STAT1n_unp = u[6]   # Nuclear unphosphorylated STAT1
    SOCS1 = u[7]        # SOCS1
    IRF1_mRNA = u[8]    # IRF1 mRNA
    STAT1_total = u[9]  # Total STAT1
    IFNg_R = u[10]      # Free IFNγ receptor
    STAT1c_unp = u[11]  # Cytoplasmic unphosphorylated STAT1
    
    # Direct reactions (no delay)
    v1 = -k0 * IFNg                           # Degradation of IFNg
    v2 = k1 * IFNg * IFNg_R - k2 * IFNRJ      # Formation/dissociation of IFNRJ complex
    v4 = k6 * STAT1c                          # Transport of STAT1c to nucleus
    v5 = k7 * STAT1n - k8 * STAT1Dn           # Binding/unbinding of STAT1n to DNA
    v6 = k5 * STAT1n                          # Dephosphorylation of STAT1n
    v7 = 3 * k11 * STAT1c_unp - k12 * STAT1n_unp  # Transport and dephosphorylation
    
    # Distributed delayed reactions using gamma kernels
    # 1. STAT1Dn affecting SOCS1 production (kernel 1)
    STAT1Dn_response = compute_distributed_delay(gamma_kernel_1, 5, h, p, t)
    v8 = k13 + k9 * STAT1Dn_response - k10 * SOCS1
    
    # 2. STAT1Dn affecting IRF1 mRNA production (kernel 2)
    STAT1Dn_response2 = compute_distributed_delay(gamma_kernel_2, 5, h, p, t)
    v9 = k15 + k16 * STAT1Dn_response2 - k17 * IRF1_mRNA
    
    # 3. SOCS1 inhibiting STAT1c phosphorylation (kernel 3)
    SOCS1_response = compute_distributed_delay(gamma_kernel_3, 7, h, p, t)
    v3 = k4 * IFNRJ * STAT1c_unp / (1 + k14 * SOCS1_response)
    
    # 4. IFNRJ affecting Total STAT1 (kernel 4)
    IFNRJ_response = compute_distributed_delay(gamma_kernel_4, 2, h, p, t)
    v10 = k3 * IFNRJ_response
    
    # Differential equations
    du[1] = v1 - v2      # IFN gamma
    du[2] = v2           # IFNRJ
    du[3] = v3 - v4      # Cytoplasmic STAT1 Dimer
    du[4] = 3*v4 - v5 - v6  # Free Nuclear STAT1 Dimer
    du[5] = v5           # Nuclear STAT1 Dimer bound to DNA
    du[6] = v6 + v7      # Nuclear unphosphorylated STAT1
    du[7] = v8           # SOCS1
    du[8] = v9           # IRF1 mRNA
    du[9] = v10          # Total STAT1
    du[10] = -v2         # Free IFNγ receptor
    du[11] = (v10 - (v6+3*v3+v4+v5))/3  # Cytoplasmic unphosphorylated STAT1
end

# Initial conditions (same as in original model)
u0 = [
    100.0,  # IFN gamma
    0.0,    # IFNRJ - complex IFNγ and receptor
    0.005,  # Cytoplasmic STAT1 Dimer
    0.0,    # Free Nuclear STAT1 Dimer
    0.0,    # Nuclear STAT1 Dimer bound to DNA
    0.5,    # Nuclear unphosphorylated STAT1
    73.6,   # SOCS1
    0.0,    # IRF1 mRNA
    1.0,    # Total STAT1
    I_tot,  # Free IFNγ receptor
    (0.3-(0.2+3*0.005+0.0+0.0))/3  # Cytoplasmic unphosphorylated STAT1
]

# Set up the DDE problem with distributed delays
println("Setting up DDE problem with distributed delays...")
prob = DDEProblem(dde_model!, u0, history_function, tspan, nothing)

# Create a progress meter for the solver
println("Solving DDE problem...")
p = Progress(100, desc="Solving JAK-STAT model: ", dt=1.0)

# Create callbacks for progress tracking and history storage
prog_times = range(0, T, length=100)
prog_callback = PresetTimeCallback(prog_times, integrator -> begin
    next!(p, showvalues = [(:time, "$(round(integrator.t/60, digits=2)) hours")])
    return nothing
end)

# Create callback to store history at each saved time point
save_times = 0:dt:T
save_callback = PresetTimeCallback(save_times, store_history!)

# Combine callbacks
cb = CallbackSet(prog_callback, save_callback)

sol = solve(prob, MethodOfSteps(Tsit5()), dt=dt, 
            callback=cb,
            saveat=save_times, reltol=1e-4, abstol=1e-6)
finish!(p)

# Extract solution for plotting
x_sol = Array(sol)

# Plot results - reusing your plot function
function plot_all_species(x_hist, sol)
    # Get time points from solution object
    t_points = sol.t
    
    # Convert time to hours for better readability
    t_hours = t_points ./ 60
    
    # Species names for better labeling
    species_names = [
        "IFN gamma", 
        "IFNRJ complex", 
        "Cytoplasmic STAT1 Dimer",
        "Free Nuclear STAT1 Dimer",
        "Nuclear STAT1 Dimer bound to DNA",
        "Nuclear unphosphorylated STAT1",
        "SOCS1", 
        "IRF1 mRNA", 
        "Total STAT1",
        "Free IFNγ receptor",
        "Cytoplasmic unphosphorylated STAT1"
    ]
    
    # Calculate derived quantities
    num_points = length(t_points)
    RSNC = zeros(num_points)
    STAT1p = zeros(num_points)
    
    for i in 1:num_points
        # RSNC = (x4+x5+x6)/(3*(x3+x11)) - Ratio of nuclear to cytoplasmic STAT1
        RSNC[i] = (x_hist[4,i] + x_hist[5,i] + x_hist[6,i]) / (3 * (x_hist[3,i] + x_hist[11,i]))
        
        # STAT1p = (3*x3+x4+x5) - Phosphorylated STAT1
        STAT1p[i] = 3 * x_hist[3,i] + x_hist[4,i] + x_hist[5,i]
    end
    
    # Create a 4×3 grid of plots
    font_settings = (fontfamily="sans-serif", titlefontsize=10, guidefontsize=8, legendfontsize=6)
    plots = []
    
    for i in 1:11
        p = plot(t_hours, x_hist[i, :], 
                 label="", 
                 title=species_names[i],
                 xlabel=(i > 8) ? "Time (h)" : "",
                 ylabel="Concentration",
                 linewidth=1.5;
                 font_settings...)
        push!(plots, p)
    end
    
    # Add the derived quantities plots
    p_rsnc = plot(t_hours, RSNC, 
                 label="", 
                 title="RSNC - N/C STAT1 Ratio",
                 xlabel="Time (h)",
                 ylabel="Ratio",
                 linewidth=1.5,
                 color=:darkblue;
                 font_settings...)
    push!(plots, p_rsnc)
    
    p_stat1p = plot(t_hours, STAT1p, 
                  label="", 
                  title="STAT1p - Phosphorylated STAT1",
                  xlabel="Time (h)",
                  ylabel="Concentration",
                  linewidth=1.5,
                  color=:darkred;
                  font_settings...)
    push!(plots, p_stat1p)
    
    # Create the grid plot with 5×3 layout to accommodate the 13 plots
    all_plot = plot(plots..., layout=(5, 3), size=(900, 1000), 
                    plot_title="JAK-STAT Pathway with Gamma-Distributed Delays")
    
    display(all_plot)
    savefig(all_plot, joinpath(@__DIR__, "distributed_delay_model_plot.png"))
    
    # Also create a separate plot just for RSNC and STAT1p for better visibility
    derived_plot = plot(layout=(2, 1), size=(800, 600),
                        plot_title="Derived Quantities in JAK-STAT Pathway")
    
    plot!(derived_plot[1], t_hours, RSNC, 
          label="RSNC", 
          title="Nuclear to Cytoplasmic STAT1 Ratio",
          xlabel="",
          ylabel="Ratio",
          linewidth=2,
          color=:darkblue)
          
    plot!(derived_plot[2], t_hours, STAT1p, 
          label="STAT1p", 
          title="Phosphorylated STAT1",
          xlabel="Time (h)",
          ylabel="Concentration",
          linewidth=2,
          color=:darkred)
    
    display(derived_plot)
    savefig(derived_plot, joinpath(@__DIR__, "derived_quantities_plot.png"))
    
    return all_plot, derived_plot
end

println("Plotting results...")
# Extract solution for plotting using the actual time points from solution
x_sol = Array(sol)
all_species_plot, derived_plot = plot_all_species(x_sol, sol)

println("Completed!")