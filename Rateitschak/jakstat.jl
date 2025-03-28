using Catalyst
using DifferentialEquations
using BenchmarkTools
using Plots
using DelaySSAToolkit
using Distributions
using ProgressMeter
theme(:mute);

# Parameters
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

const tau_avg_1 = 39
const tau_avg_2 = 22
const tau_avg_3 = 164
const tau_avg_4 = 185

const I_tot = 0.0038

const q1 = 4.0/tau_avg_1
const q2 = 4.0/tau_avg_2
const q3 = 4.0/tau_avg_3
const q4 = 4.0/tau_avg_4

const T = 24.0*60.0
const dt = 0.1
tspan = (0.0, T)
t_points = 0:dt:T

# Kernel of Gamma function
function Gam(tau, p, q)
    gamma_dist = Gamma(p, 1/q)
    pdf_value = pdf(gamma_dist, tau)
    return pdf_value
end

function response(t, q, n, x_hist)
    integral_response = 0.0
    t_i = round(Int, t/dt) + 1  # Julia arrays are 1-indexed

    for i in 1:t_i-1
        tau = (i-1)*dt
        integral_response += Gam(tau, 4, q) * x_hist[n, t_i-i+1] * dt
    end
    
    return integral_response
end

# Reaction rates
function reaction_rates(x, t, x_hist)
    i = min(round(Int, t/dt) + 1, size(x_hist, 2))
    
    v = zeros(10)
    v[1] = -k0*x[1]
    v[2] = k1*x[1]*x[10] - k2*x[2]
    v[3] = k4*x[2]*x[11]/(1+k14*response(t, q3, 7, x_hist))
    v[4] = k6*x[3]
    v[5] = k7*x[4] - k8*x[5]
    v[6] = k5*x[4]
    v[7] = 3*k11*x[11] - k12*x[6]
    v[8] = k13 + k9*response(t, q1, 5, x_hist) - k10*x[7]
    v[9] = k15 + k16*response(t, q2, 5, x_hist) - k17*x[8]
    v[10] = k3*response(t, q4, 2, x_hist)

    return v
end

# Define the system of differential equations
function model!(dx, x, p, t)
    x_hist = p  # Pass historical x values as a parameter
    
    v = reaction_rates(x, t, x_hist)
    dx[1] = v[1] - v[2]    # IFN gamma
    dx[2] = v[2]           # IFNRJ - complex IFNγ and receptor
    dx[3] = v[3] - v[4]    # Cytoplasmic STAT1 Dimer
    dx[4] = 3*v[4] - v[5] - v[6]  # Free Nuclear STAT1 Dimer
    dx[5] = v[5]           # Nuclear STAT1 Dimer bound to DNA
    dx[6] = v[6] + v[7]    # Nuclear unphosphorylated STAT1
    dx[7] = v[8]           # SOCS1
    dx[8] = v[9]           # IRF1 mRNA
    dx[9] = v[10]          # Total STAT1
    dx[10] = -v[2]         # Free IFNγ receptor
    dx[11] = (v[9] - (v[6]+3*v[3]+v[4]+v[5]))/3  # Cytoplasmic unphosphorylated STAT1
end

# Initial conditions
x0 = [
    100.0,  # IFN gamma
    0.0,    # IFNRJ - complex IFNγ and receptor
    0.005,  # Cytoplasmic STAT1 Dimer
    0.0,    # Free Nuclear STAT1 Dimer
    0.0,    # Nuclear STAT1 Dimer bound to DNA
    0.5,    # Nuclear unphosphorylated STAT1
    73.6,    # SOCS1
    0.0,    # IRF1 mRNA
    1.0,    # Total STAT1
    I_tot,  # Free IFNγ receptor (initialized to I_tot since x[1,0] is 0)
    (0.3-(0.2+3*0.005+0.0+0.0))/3  # Cytoplasmic unphosphorylated STAT1
]

# Create a matrix to store the solution
x_hist = zeros(11, length(t_points))
x_hist[:, 1] = x0

# Custom Euler solver to handle the memory-dependent system
@showprogress "Simulating..." for i in 1:length(t_points)-1
    t = t_points[i]
    dx = zeros(11)
    model!(dx, x_hist[:, i], x_hist, t)
    x_hist[:, i+1] = x_hist[:, i] + dx * dt
end

# Plot results
font_settings = (fontfamily="sans-serif", titlefontsize=12, guidefontsize=10, legendfontsize=8)

# After simulation is complete, create plots for all species
function plot_all_species(x_hist, t_points)
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
    
    # Create a 4×3 grid of plots (with one empty space)
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
    
    # Create the grid plot
    all_plot = plot(plots..., layout=(4, 3), size=(900, 800), plot_title="JAK-STAT Pathway Species")
    
    display(all_plot)
    savefig(all_plot, joinpath(@__DIR__, "all_species_plot.png"))
    
    return all_plot
end

# Call the plotting function after simulation is complete
all_species_plot = plot_all_species(x_hist, t_points)

