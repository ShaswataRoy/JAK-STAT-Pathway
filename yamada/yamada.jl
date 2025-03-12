using Catalyst
using DifferentialEquations
using Plots
using ModelingToolkit
using CSV
using DataFrames

# Register the hill function for symbolic use
@register_symbolic hill(s, ka, kb, n)

# Define the reaction system with all species, parameters and reactions
jak_stat = @reaction_network begin
    # Define parameters that will be used in the model
    @parameters compartment_cytoplasm compartment_nucleus
    
    # v1: Receptor-JAK Binding
    (0.1, 0.05), R + JAK <--> RJ
    
    # v2: Interferon-Receptor Binding
    (0.02, 0.02), RJ + IFN <--> IFNRJ
    
    # v3: IFN-Receptor complex dimerization
    (0.04, 0.2), 2IFNRJ <--> IFNRJ2
    
    # v4: INF-Receptor complex activation
    0.005, IFNRJ2 --> IFNRJ2_star
    
    # v5: Activated INFRJ2-STAT1c binding
    (0.008, 0.8), STAT1c + IFNRJ2_star <--> IFNRJ2_star_STAT1c
    
    # v6: STAT1c activation
    0.4, IFNRJ2_star_STAT1c --> IFNRJ2_star + STAT1c_star
    
    # v7: Activated IFNRJ2-STAT1c binding
    (0.005, 0.5), IFNRJ2_star + STAT1c_star <--> IFNRJ2_star_STAT1c_star
    
    # v8: Activated STAT1c dimerization
    (0.02, 0.1), 2STAT1c_star <--> STAT1c_star_STAT1c_star
    
    # v9: SHP2 binding
    (0.001, 0.2), IFNRJ2_star + SHP2 <--> IFNRJ2_star_SHP2
    
    # v10: IFNJR2 dephosphorylation
    0.003, IFNRJ2_star_SHP2 --> IFNRJ2 + SHP2
    
    # v11: Phosphorylated STAT1c-PPX binding
    (0.001, 0.2), PPX + STAT1c_star <--> STAT1c_star_PPX
    
    # v12: STAT1c dephosphorylation
    0.003, STAT1c_star_PPX --> STAT1c + PPX
    
    # v13: PPX binding
    (0.001, 0.2), PPX + STAT1c_star_STAT1c_star <--> STAT1c_star_STAT1c_star_PPX
    
    # v14: STAT1c dimer dephosphorylation
    0.003, STAT1c_star_STAT1c_star_PPX --> STAT1c_STAT1c_star + PPX
    
    # v15: STAT1c-phosphorylated STAT1c binding
    (2.0e-7, 0.2), STAT1c + STAT1c_star <--> STAT1c_STAT1c_star
    
    # v16: STAT1c-nuclear transport
    0.005, STAT1c_star_STAT1c_star --> STAT1n_star_STAT1n_star
    
    # v17: Phosphorylated STAT1n dimerization
    (0.02, 0.1), 2STAT1n_star <--> STAT1n_star_STAT1n_star
    
    # v18: PPN binding
    (0.001, 0.2), PPN + STAT1n_star <--> STAT1n_star_PPN
    
    # v19: STAT1n dephosphorylation
    0.005, STAT1n_star_PPN --> STAT1n + PPN
    
    # v20: PPN binding
    (0.001, 0.2), PPN + STAT1n_star_STAT1n_star <--> STAT1n_star_STAT1n_star_PPN
    
    # v21: STAT1n dephosphorylation
    0.005, STAT1n_star_STAT1n_star_PPN --> STAT1n_STAT1n_star + PPN
    
    # v22: STAT1n-phosphorylated STAT1n dimerization
    (2.0e-7, 0.2), STAT1n + STAT1n_star <--> STAT1n_STAT1n_star
    
    # v23: STAT1n transport to cytoplasm
    0.05, STAT1n --> STAT1c
    
    # v24: Transcription (Hill function)
    hill(STAT1n_star_STAT1n_star, 0.01, 400.0, 1), 0 --> mRNAn
    
    # v25: mRNA transport to cytoplasm
    0.001, mRNAn --> mRNAc
    
    # v26: SOCS1 synthesis
    0.01, mRNAc --> mRNAc + SOCS1
    
    # v27: mRNAc degradation
    5.0e-4, mRNAc --> 0
    
    # v28: SOCS1 degradation
    5.0e-4, SOCS1 --> 0
    
    # v29: phosphorylated IFNRJ2-SOCS1 binding
    (0.02, 0.1), SOCS1 + IFNRJ2_star <--> IFNRJ2_star_SOCS1
    
    # v30: STAT1c binding
    (0.008, 0.8), STAT1c + IFNRJ2_star_SOCS1 <--> IFNRJ2_star_SOCS1_STAT1c
    
    # v31: SHP2 binding
    (0.001, 0.2), SHP2 + IFNRJ2_star_SOCS1_STAT1c <--> IFNRJ2_star_SHP2_SOCS1_STAT1c
    
    # v32: IFNRJ2 dephosphorylation
    0.003, IFNRJ2_star_SHP2_SOCS1_STAT1c --> IFNRJ2 + SHP2 + SOCS1 + STAT1c
    
    # v33: SOCS1 unbinding
    5.0e-4, IFNRJ2_star_SHP2_SOCS1_STAT1c --> IFNRJ2_star_SHP2_STAT1c
    
    # v34: SHP2 binding
    (0.001, 0.2), SHP2 + IFNRJ2_star_SOCS1 <--> IFNRJ2_star_SHP2_SOCS1
    
    # v35: STAT1c binding
    (0.008, 0.8), STAT1c + IFNRJ2_star_SHP2_SOCS1 <--> IFNRJ2_star_SHP2_SOCS1_STAT1c
    
    # v36: SHP2 binding
    (0.001, 0.2), SHP2 + IFNRJ2_star_STAT1c <--> IFNRJ2_star_SHP2_STAT1c
    
    # v37: IFNRJ2 dephosphorylation
    0.003, IFNRJ2_star_SHP2_STAT1c --> IFNRJ2 + SHP2 + STAT1c
    
    # v38: SOCS1 unbinding
    5.0e-4, IFNRJ2_star_SOCS1_STAT1c --> IFNRJ2_star_STAT1c + SOCS1
    
    # v39: SOCS1 unbinding
    5.0e-4, IFNRJ2_star_SHP2_SOCS1 --> IFNRJ2_star_SHP2
    
    # v40: IFNRJ2 dephosphorylation
    0.003, IFNRJ2_star_SHP2_SOCS1 --> IFNRJ2 + SHP2 + SOCS1
    
    # v41: SOCS1 unbinding
    5.0e-4, IFNRJ2_star_SOCS1 --> IFNRJ2_star
    
    # v42: SOCS1 binding
    (0.02, 0.1), SOCS1 + IFNRJ2_star_STAT1c <--> IFNRJ2_star_SOCS1_STAT1c
    
    # v43: SOCS1 binding
    (0.02, 0.1), SOCS1 + IFNRJ2_star_SHP2 <--> IFNRJ2_star_SHP2_SOCS1
    
    # v44: SOCS1 binding
    (0.02, 0.1), SOCS1 + IFNRJ2_star_SHP2_STAT1c <--> IFNRJ2_star_SHP2_SOCS1_STAT1c
    
    # v45: Interferon-receptor binding
    (0.02, 0.02), IFN + R <--> IFNR
    
    # v46: IFNR-JAK binding
    (0.1, 0.05), IFNR + JAK <--> IFNRJ

    # IRF1 dynamics
    (0.01),  STAT1n_STAT1n_star --> STAT1n_STAT1n_star + IRF
    5e-4, IRF --> 0
end

# Define initial conditions
u0 = [
    :R => 10.0,            # Receptor
    :JAK => 10.0,          # JAK
    :RJ => 0.0,            # Receptor JAK complex
    :IFNRJ => 0.0,         # Interferon-Receptor-JAK complex
    :IFNRJ2 => 0.0,        # IFNRJ dimer
    :IFNRJ2_star => 0.0,   # Activated IFNRJ complex
    :STAT1c => 1000.0,     # STAT1c
    :IFNRJ2_star_STAT1c => 0.0,
    :STAT1c_star => 0.0,
    :IFNRJ2_star_STAT1c_star => 0.0,
    :STAT1c_star_STAT1c_star => 0.0,
    :SHP2 => 100.0,
    :IFNRJ2_star_SHP2 => 0.0,
    :PPX => 50.0,
    :STAT1c_star_PPX => 0.0,
    :STAT1c_STAT1c_star => 0.0,
    :STAT1n_star_STAT1n_star => 0.0,
    :STAT1n_star => 0.0,
    :PPN => 60.0,
    :STAT1n_star_PPN => 0.0,
    :STAT1n => 0.0,
    :STAT1n_STAT1n_star => 0.0,
    :mRNAn => 0.0,
    :mRNAc => 0.0,
    :SOCS1 => 0.0,
    :IFNRJ2_star_SOCS1 => 0.0,
    :IFNRJ2_star_SHP2_SOCS1_STAT1c => 0.0,
    :STAT1c_star_STAT1c_star_PPX => 0.0,
    :STAT1n_star_STAT1n_star_PPN => 0.0,
    :IFNRJ2_star_SOCS1_STAT1c => 0.0,
    :IFNRJ2_star_SHP2_STAT1c => 0.0,
    :IFNRJ2_star_SHP2_SOCS1 => 0.0,
    :IFNR => 0.0,
    :IFN => 100.0,
    :IRF => 0.0
]

# Define parameters
p = [
    :compartment_cytoplasm => 1.0,
    :compartment_nucleus => 1.0
]

# Define problem
tspan = (0.0, 10*3600.0)
prob = ODEProblem(jak_stat, u0, tspan, p)


# This function can be used to replicate the original MATLAB code's behavior
function simulate_jak_stat()
    @time sol = solve(prob, Tsit5(), abstol=1e-3)
    
    # Convert time to hours
    time_in_hours = sol.t ./ 3600
    
    # Create DataFrame with time and IRF data
    df = DataFrame(
        time_hours = time_in_hours,
        irf = sol[Symbol("IRF")]
    )
    
    # Save to CSV
    output_path = joinpath(@__DIR__, "irf_data.csv")
    CSV.write(output_path, df)
    
    # Continue with existing plotting code
    plt = plot(layout=(6,6), size=(1200,1000), legend=false)
    species_list = [
        :R, :JAK, :RJ, :IFNRJ, :IFNRJ2, :IFNRJ2_star, :STAT1c, :IFNRJ2_star_STAT1c, :STAT1c_star,
        :IFNRJ2_star_STAT1c_star, :STAT1c_star_STAT1c_star, :SHP2, :IFNRJ2_star_SHP2,
        :PPX, :STAT1c_star_PPX, :STAT1c_STAT1c_star, :STAT1n_star_STAT1n_star, :STAT1n_star,
        :PPN, :STAT1n_star_PPN, :STAT1n, :STAT1n_STAT1n_star, :mRNAn, :mRNAc, :SOCS1,
        :IFNRJ2_star_SOCS1, :IFNRJ2_star_SHP2_SOCS1_STAT1c, :STAT1c_star_STAT1c_star_PPX,
        :STAT1n_star_STAT1n_star_PPN, :IFNRJ2_star_SOCS1_STAT1c, :IFNRJ2_star_SHP2_STAT1c,
        :IFNRJ2_star_SHP2_SOCS1, :IFNR, :IFN, :IRF
    ]
    
    for (i, sp) in enumerate(species_list)
        if i <= length(species_list)
            plot!(plt[i], time_in_hours, sol[sp], title=string(sp), titlefontsize=8)
        end
    end
    
    display(plt)
    savefig(plt, joinpath(@__DIR__,"jak_stat_all_species.png"))
end

# Call the simulation function
simulate_jak_stat()