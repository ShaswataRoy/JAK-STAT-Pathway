using CSV
using DataFrames
using Catalyst
using DifferentialEquations
using Plots
using Interpolations

# Load the IRF data
irf_data = CSV.read(joinpath(@__DIR__, "irf_data.csv"), DataFrame)

# Create interpolation function for IRF values
irf_interpolation = LinearInterpolation(irf_data.time_hours, irf_data.irf)

# Define the reaction system
@parameters k_deg k_prod
@variables t
@species mRNA(t)

# Create reaction network with time-dependent IRF
rx = @reaction_network begin
    ($irf_interpolation)(t), 0 --> mRNA    # Production driven by IRF (note the $ prefix)
    k_deg, mRNA --> 0                   # Degradation
end

# Create parameter function for time-dependent production rate
function paramfun(p,t,u)
    p[1] = irf_interpolation(t)  # k_prod varies with time
    return p
end

# Set parameters and initial conditions
p = [
    k_prod => 1.0,  # Initial value will be updated by paramfun
    k_deg => 0.005    # Degradation rate constant
]
u0 = [
    mRNA => 0.0   # Initial mRNA concentration
]

# Define timespan (matching the IRF data timespan)
tspan = (0.0, maximum(irf_data.time_hours))

# Create and solve the problem
prob = ODEProblem(rx, u0, tspan, p, paramfun=paramfun)
sol = solve(prob, Tsit5())

# Plot the mRNA time course
plt = plot(sol, 
            vars=[mRNA], 
            lw=2, 
            xlabel="Time (hours)", 
            ylabel="mRNA Concentration", 
            title="mRNA Expression Dynamics",
            label="mRNA")

# Save the plot
savefig(plt,joinpath(@__DIR__,"mRNA_dynamics.png"))