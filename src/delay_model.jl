using Catalyst
using DifferentialEquations
using BenchmarkTools
using Plots
using DelaySSAToolkit
using Distributions

# Import constants
include("constants.jl")

# Set plot theme
theme(:mute)
# Define the reactions with simplified rate constants
rxs = @reaction_network begin
    # Degradation of IFNg
    k0, IFNg --> 0
    # Binding of IFNg to Ir to form IIr
    k1, IFNg + Ir --> IIr
    # Dissociation of IIr to form IFNg and Ir
    k2, IIr --> IFNg + Ir
    # Production of STAT1_int from IIr
    k3, IIr --> IIr + STAT1_int

    # Translocation of STAT1Dn to STAT1Un
    k5, STAT1Dn --> STAT1Un
    # Translocation of STAT1Dc to STAT1Dn
    k6, STAT1Dc --> STAT1Dn
    # Binding of STAT1Dn to form STAT1Db
    k7, STAT1Dn --> STAT1Db
    # Degradation of STAT1Db
    k8, STAT1Db --> 0
    
    # STAT1Db delay triggers SOCS1 production
    k9, STAT1Db --> STAT1Db + SOCS1_int
    # Creation of SOCS1
    k13, 0 --> SOCS1
    # Degradation of SOCS1
    k10, SOCS1 --> 0

    # Translocation of STAT1Uc to STAT1Un
    k11, STAT1Uc --> STAT1Un
    # Degradation of STAT1Un
    k12, STAT1Un --> 0

    # STAT1Db delay triggers IRF1 production
    k15, STAT1Db --> STAT1Db + IRF1_int
    # Creation of IRF1
    k16, 0 --> IRF1
    # Degradation of IRF1
    k17, IRF1 --> 0

    # Michaelis-Menten reaction with inhibition: Binding of IIr and STAT1Uc to form IIr_STAT1Uc
    ke1, IIr + STAT1Uc --> IIr_STAT1Uc
    ke2, IIr_STAT1Uc --> IIr + STAT1Uc
    # Dissociation of IIr_STAT1Uc to form IIr and STAT1Dc
    kp, IIr_STAT1Uc --> IIr + STAT1Dc

    # Inhibition reaction: Binding of SOCS1 and IIr to form SOCS1_IIr
    ki1, SOCS1 + IIr --> SOCS1_IIr_int
    # Dissociation of SOCS1_IIr to form SOCS1 and IIr
    ki2, SOCS1_IIr --> SOCS1 + IIr

    # Degradation of STAT1
    k18, STAT1 --> 0
end

# Convert to JumpSystem
jumpsys_jakstat = convert(JumpSystem, rxs; combinatoric_ratelaws=false)

# Get the state variables from the system
state_vars = states(jumpsys_jakstat)

# Print state variables to verify names
println("State variables in system:")
foreach(println, state_vars)

# Delay in STAT1 production due to IIr
delay_trigger_affect1! = function (integrator, rng)
    τ = rand(Gamma(p, q4))
    append!(integrator.de_chan[1], τ)
end

# Delay in SOCS1 production due to STAT1Db
delay_trigger_affect2! = function (integrator, rng)
    τ = rand(Gamma(p, q1))
    append!(integrator.de_chan[2], τ)
end

# Delay in IRF1 production due to STAT1Db
delay_trigger_affect3! = function (integrator, rng)
    τ = rand(Gamma(p, q2))
    append!(integrator.de_chan[3], τ)
end

# Delay in SOCS1 production due to IIr
delay_trigger_affect4! = function (integrator, rng)
    τ = rand(Gamma(p, q3))
    append!(integrator.de_chan[4], τ)
end

delay_trigger = Dict(4 => delay_trigger_affect1!,
                     9 => delay_trigger_affect2!,
                     14 => delay_trigger_affect3!,
                     20 => delay_trigger_affect4!)
                     
delay_complete = Dict(1 => [4 => -1,17 => 1],
                      2 => [9 => -1,10 => 1],
                      3 => [14 => -1,13 => 1],
                      4 => [15 => -1,16 => 1])

delay_interrupt = Dict()

delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)

# Create initial condition mapping
u0_map = [
    Symbol(state_vars[1]) => x1_0,
    Symbol(state_vars[2]) => x10_0,
    Symbol(state_vars[3]) => x2_0,
    Symbol(state_vars[4]) => 0.0,
    Symbol(state_vars[5]) => x4_0,
    Symbol(state_vars[6]) => x6_0,
    Symbol(state_vars[7]) => x3_0,
    Symbol(state_vars[8]) => x5_0,
    Symbol(state_vars[9]) => 0.0,
    Symbol(state_vars[10]) => x7_0,
    Symbol(state_vars[11]) => x8_0,
    Symbol(state_vars[12]) => 0.0,
    Symbol(state_vars[13]) => x9_0,
    Symbol(state_vars[14]) => 0.0,
    Symbol(state_vars[15]) => 0.0,
    Symbol(state_vars[16]) => 0.0,
    Symbol(state_vars[17]) => x11_0
]

ps = [k0, k1, k2, k3, k5, k6, k7, k8, k9, k13, k10, k11, k12, k15, k16, k17, ke1, ke2, kp, ki1, ki2, k18]
tspan = (0.0, 1000.)
dprob = DiscreteProblem(jumpsys_jakstat, u0_map, tspan, ps)

de_chan0 = [[],[],[],[]]
algo = DelayRejection()
jprob = DelayJumpProblem(
    jumpsys_jakstat, dprob, algo, delayjumpset, de_chan0; save_positions=(true, true)
)
sol = solve(jprob, SSAStepper(); seed=1234);

# Define initial conditions and time span
# Add your simulation code here

# Create plots
function plot_results(sol)
    p1 = plot(sol, idxs=[12], label="STAT1", 
        xlabel="Time (min)", ylabel="Molecules", 
        title="STAT1 Dynamics")

    p2 = plot(sol, idxs=[10], label="SOCS1", 
        xlabel="Time (min)", ylabel="Molecules", 
        title="SOCS1 Dynamics")

    plot(p1, p2, layout=(2,1), size=(800,600))
end

plot_results(sol)
# Run simulation and plot results
# Add your simulation execution code here 