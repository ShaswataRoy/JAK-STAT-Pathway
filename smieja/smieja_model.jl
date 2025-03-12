# smieja_model.jl
using Catalyst
using DifferentialEquations
using Plots
using ModelingToolkit
using Sundials  # Add this line

# Register any special functions needed
@register_symbolic hill(s, ka, kb, n)

# Define the JAK-STAT pathway model from Smieja
@variables t
rn = @reaction_network begin
    # Define compartment parameters
    @parameters compartment_cytoplasm=1.0 compartment_nucleus=1.0
    
    # Define model parameters
    @parameters ks1phos ks1dephc ks2phos ks2dephc ks1phos_sat ks2phos_sat
    @parameters k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14
    
    # IFN receptor binding and activation
    (k1, k2), IFN + R <--> IIr        # Binding and dissociation
    k3, IIr --> 0                     # Degradation
    
    # STAT reactions
    k4, IIr + STAT1Uc --> IIr + STAT1Dc  # STAT1 phosphorylation
    k5, STAT1Dc --> STAT1Dn              # STAT1 nuclear transport
    k6, STAT1Dn --> STAT1Uc              # STAT1 dephosphorylation and export
    
    # IRF-1 related reactions
    k7, STAT1Dn --> STAT1Dn + IRF        # IRF transcription
    k8, IRF --> 0                        # IRF degradation
    
    # SOCS related reactions
    k9, STAT1Dn --> STAT1Dn + SOCS       # SOCS transcription
    k10, SOCS --> 0                      # SOCS degradation
    
    # Inhibition mechanism
    (k11, k12), SOCS + IIr <--> SOCS_IIr # SOCS binding to IIr
    k13, SOCS_IIr --> SOCS               # IIr degradation when bound to SOCS
    
    # STAT2 reactions (if included in the model)
    k14, IIr + STAT2 --> IIr + STAT2P    # STAT2 phosphorylation
end

# Convert network to ODESystem
osys = convert(ODESystem, rn)

# Function to run a simulation with JAK-STAT model
function run_smieja_simulation(tspan=(0.0, 3.0*3600.0); ifn_dose=10.0)
    println("Running simulation with IFN dose = $ifn_dose")
    
    # Calculate time-dependent parameters
    tphos1 = 15  # in minutes
    ts1pc = 15   # in minutes
    tphos2 = 15  # for STAT2
    
    ks1phos_val = log(2)/tphos1/60  # phosphorylation rate for STAT1
    ks1dephc_val = log(2)/ts1pc/60  # dephosphorylation of STAT1 in cytoplasm
    ks2phos_val = log(2)/tphos2/60  # phosphorylation rate for STAT2
    ks2dephc_val = log(2)/ts1pc/60  # dephosphorylation of STAT2 in cytoplasm
    
    # Parameter values
    p = [
        1.0,        # compartment_cytoplasm
        1.0,        # compartment_nucleus
        ks1phos_val,# ks1phos
        ks1dephc_val,# ks1dephc
        ks2phos_val,# ks2phos
        ks2dephc_val,# ks2dephc
        1e4,        # ks1phos_sat
        1e4,        # ks2phos_sat
        0.1,        # k1 - IFN + R binding rate
        0.05,       # k2 - IIr dissociation rate
        0.003,      # k3 - IIr degradation
        0.008,      # k4 - STAT1 phosphorylation
        0.005,      # k5 - Nuclear import
        0.05,       # k6 - Nuclear export
        0.01,       # k7 - IRF transcription
        0.0005,     # k8 - IRF degradation
        0.01,       # k9 - SOCS transcription
        0.0005,     # k10 - SOCS degradation
        0.02,       # k11 - SOCS-IIr binding
        0.1,        # k12 - SOCS-IIr dissociation
        0.003,      # k13 - IIr degradation when bound to SOCS
        0.008       # k14 - STAT2 phosphorylation
    ]
    
    # Print parameters for debugging
    param_names = ["compartment_cytoplasm", "compartment_nucleus", "ks1phos", "ks1dephc", "ks2phos", "ks2dephc", 
                 "ks1phos_sat", "ks2phos_sat", "k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", 
                 "k10", "k11", "k12", "k13", "k14"]
    
    println("Parameters:")
    for (i, name) in enumerate(param_names)
        println("  $name = $(p[i])")
    end
    
    # Get species names
    species_names = ["IFN", "R", "IIr", "STAT1Uc", "STAT1Dc", "STAT1Dn", "IRF", "SOCS", "SOCS_IIr", "STAT2", "STAT2P"]
    
    # Initial conditions
    u0 = zeros(length(species_names))
    u0[1] = ifn_dose  # IFN
    u0[2] = 10.0       # R
    u0[4] = 1000.0    # STAT1Uc
    u0[10] = 0.0    # STAT2
    
    println("Initial conditions:")
    for (i, name) in enumerate(species_names)
        println("  $name = $(u0[i])")
    end
    
    # Create the ODE problem with stricter tolerances
    prob = ODEProblem(osys, u0, tspan, p)

    # Solve using CVODE_BDF (stiff solver) with adjusted tolerances
    println("Solving ODE system...")
    @time sol = solve(prob, CVODE_BDF(), 
                     abstol=1e-8, 
                     reltol=1e-6,
                     maxiters=Int(1e6));
    println("Simulation complete.")
    
    return sol, species_names
end

# Function to plot results
function plot_smieja_results(sol, species_names)
    println("Creating plots...")
    
    # Create a layout for plotting all species
    n_species = length(species_names)
    n_rows = ceil(Int, n_species / 4)
    plt = plot(layout=(n_rows, 4), size=(1200, 300*n_rows), legend=:topright)
    
    # Plot each species
    for (i, name) in enumerate(species_names)
        plot!(plt, sol.t/3600.0, sol[i, :], subplot=i, label=name, 
              xlabel="Time (hour)", ylabel="Concentration")
        title!(plt, name, subplot=i)
    end
    
    # Save plot
    savefig(plt, joinpath(@__DIR__, "smieja_jak_stat_all_species.png"))
    println("Saved plot: smieja/smieja_jak_stat_all_species.png")
    
    return plt
end

# Main function to execute the simulation with different scenarios
function main()
    println("=============================================")
    println("JAK-STAT Pathway Simulation (Smieja Model)")
    println("=============================================")
    
    println("\nRunning simulation with constant IFNγ level of 100...")
    sol, species_names = run_smieja_simulation(ifn_dose=100.0)
    print(sol.t[end])
    plot_smieja_results(sol, species_names)
    
    println("\nSimulation completed successfully!")
end

# Show that the script is ready to run
println("Smieja model script loaded. Run main() to execute simulations.")