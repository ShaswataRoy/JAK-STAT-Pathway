using Catalyst
using JumpProcesses
using Plots
using ModelingToolkit
using CSV
using DataFrames
using ProgressMeter
using Base.GC
using Statistics
using Profile
using StochasticDiffEq
using StochasticDiffEq

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
    (0.01),  STAT1n_STAT1n_star --> STAT1n_STAT1n_star + IRF1
    5e-4, IRF1 --> 0
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
    :IRF1 => 0.0
]

# Define parameters
p = [
    :compartment_cytoplasm => 1.0,
    :compartment_nucleus => 1.0
]

# Define simulation timespan (24 hours in seconds)
tspan = (0.0, 24*3600.0)

# Define problem as discrete jump process
disc_prob = DiscreteProblem(jak_stat, u0, tspan, p)

# This function can be used to replicate the original MATLAB code's behavior
function simulate_jak_stat()
    # Measure memory before simulation
    GC.gc()  # Force garbage collection
    mem_before = Sys.total_memory() - Sys.free_memory()

    # Sample every 6 minutes (360 seconds) over 24 hours
    time_points = 0.0:360.0:tspan[2]  # This already samples every 6 minutes (360 seconds)
    
    # Convert time to hours for plotting
    time_in_hours = time_points ./ 3600

    # Define number of ensemble simulations
    ensembles = 1000

    # Get the species by name from the states
    species_names = string.(states(jak_stat))
    IRF1_index = findfirst(x -> x == "IRF1(t)", species_names)

    # Initialize arrays for IRF1 data
    IRF1_values = zeros(length(time_points), ensembles)
    averaged_IRF1 = zeros(length(time_points))
    std_dev_IRF1 = zeros(length(time_points))
    
    # Progress bar
    p = Progress(ensembles, desc="Running simulations: ")
    
    for i in 1:ensembles
        # Method 1: Use SSAStepper with a SortingDirect jump problem (reliable)
        local_jump_prob = JumpProblem(jak_stat, disc_prob, SortingDirect())
        sol = solve(local_jump_prob, SSAStepper())
        
        # Method 2: Use SimpleTauLeaping (requires special setup)
        # local_jump_prob = JumpProblem(jak_stat, disc_prob, Direct(), save_positions=(false,false))
        # sol = solve(local_jump_prob, SimpleTauLeaping(), dt=300)
        
        # Sample only at the specific time points
        for (j, t) in enumerate(time_points)
            sol_at_t = sol(t)
            IRF1_values[j, i] = sol_at_t[IRF1_index] # Use the index directly
        end
        
        # Force cleanup to reduce memory usage
        sol = nothing
        if i % 10 == 0
            GC.gc()  # Periodic garbage collection
        end
        
        next!(p)
    end
    
    # Calculate mean and standard deviation
    for j in 1:length(time_points)
        averaged_IRF1[j] = mean(IRF1_values[j, :])
        std_dev_IRF1[j] = std(IRF1_values[j, :])
    end
    
    # Measure memory after simulation
    GC.gc()
    mem_after = Sys.total_memory() - Sys.free_memory()
    mem_used = (mem_after - mem_before) / (1024^2)  # Convert to MB
    
    println("Memory used for simulation: approximately $(round(mem_used, digits=2)) MB")
    
    # Plot only IRF1 with error bands - cleaner version
    plt = plot(time_in_hours, averaged_IRF1, 
    ribbon=std_dev_IRF1, 
    fillalpha=0.2,
    lw=3, 
    xlabel="Time (hours)", 
    ylabel="IRF1 Concentration", 
    title="IRF1 Expression Dynamics",
    label="Mean IRF1",
    legend=:topright,
    linecolor=:blue,
    fillcolor=:blue,
    framestyle=:box,
    grid=false,
    fontfamily="Arial",
    dpi=300)
    
    display(plt)
    savefig(plt, joinpath(@__DIR__,"IRF1_dynamics.png"))
end

function simulate_and_plot_all_species()
    # Measure memory before simulation
    GC.gc()
    mem_before = Sys.total_memory() - Sys.free_memory()

    # Sample every 6 minutes (360 seconds) over 24 hours
    time_points = 0.0:60.0:tspan[2]  
    time_in_hours = time_points ./ 3600

    # Define number of ensemble simulations (use a small number for testing first)
    ensembles = 1000

    # Get all species from the model
    species_names = string.(states(jak_stat))
    num_species = length(species_names)
    
    # Initialize arrays for all species data
    species_values = zeros(length(time_points), num_species, ensembles)
    averaged_species = zeros(length(time_points), num_species)
    std_dev_species = zeros(length(time_points), num_species)
    
    # Find the IRF1 index for easier access
    IRF1_index = findfirst(x -> x == "IRF1(t)", species_names)
    
    # Progress bar
    p = Progress(ensembles, desc="Running simulations: ")
    
    for i in 1:ensembles
        # Create jump problem
        local_jump_prob = JumpProblem(jak_stat, disc_prob, DirectCR())
        sol = solve(local_jump_prob, SSAStepper())
        
        # Sample all species at each time point
        for (j, t) in enumerate(time_points)
            sol_at_t = sol(t)
            species_values[j, :, i] = sol_at_t
        end
        
        # Force cleanup
        sol = nothing
        if i % 5 == 0
            GC.gc()
        end
        
        next!(p)
    end
    
    # Calculate statistics for all species
    for j in 1:length(time_points)
        for k in 1:num_species
            averaged_species[j, k] = mean(species_values[j, k, :])
            std_dev_species[j, k] = std(species_values[j, k, :])
        end
    end
    
    # Create output directory
    output_dir = joinpath(@__DIR__, "output")
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    
    # Export all IRF1 trajectories to CSV
    IRF1_df = DataFrame()
    IRF1_df.Time_hours = time_in_hours
    
    # Add each simulation as a column
    for i in 1:ensembles
        IRF1_df[!, "Simulation_$(i)"] = species_values[:, IRF1_index, i]
    end
    
    # Add mean and standard deviation columns
    IRF1_df.Mean = averaged_species[:, IRF1_index]
    IRF1_df.StdDev = std_dev_species[:, IRF1_index]
    
    # Save to CSV
    CSV.write(joinpath(output_dir, "IRF1_trajectories.csv"), IRF1_df)
    println("IRF1 trajectories saved to: $(joinpath(output_dir, "IRF1_trajectories.csv"))")
    
    # Continue with the rest of your existing function...
    # Define human-readable names for species groups
    species_groups = Dict(
        "Receptor Complexes" => [
            "R", "JAK", "RJ", "IFN", "IFNR", "IFNRJ", 
            "IFNRJ2", "IFNRJ2_star"
        ],
        "Cytoplasmic STAT1" => [
            "STAT1c", "STAT1c_star", "STAT1c_star_STAT1c_star", 
            "STAT1c_STAT1c_star", "STAT1c_star_PPX"
        ],
        "Nuclear STAT1" => [
            "STAT1n", "STAT1n_star", "STAT1n_star_STAT1n_star", 
            "STAT1n_STAT1n_star", "STAT1n_star_PPN"
        ],
        "Regulatory Proteins" => [
            "SHP2", "PPX", "PPN", "SOCS1"
        ],
        "Gene Expression Products" => [
            "mRNAn", "mRNAc", "IRF1"
        ],
        "IFNRJ2 Receptor Complexes" => [
            "IFNRJ2_star_STAT1c", "IFNRJ2_star_SHP2", 
            "IFNRJ2_star_SOCS1", "IFNRJ2_star_SOCS1_STAT1c"
        ],
        "Multi-Protein Complexes" => [
            "STAT1c_star_STAT1c_star_PPX", "STAT1n_star_STAT1n_star_PPN",
            "IFNRJ2_star_SHP2_STAT1c", "IFNRJ2_star_SHP2_SOCS1"
        ]
    )
    
    # Define display names for better plot labels
    display_names = Dict(
        "R" => "Receptor",
        "JAK" => "JAK",
        "RJ" => "Receptor-JAK",
        "IFN" => "Interferon",
        "IFNR" => "IFN-Receptor",
        "IFNRJ" => "IFN-R-JAK",
        "IFNRJ2" => "IFN-R-JAK₂",
        "IFNRJ2_star" => "IFN-R-JAK₂*",
        "STAT1c" => "STAT1 (cyto)",
        "STAT1c_star" => "STAT1* (cyto)",
        "STAT1c_star_STAT1c_star" => "STAT1*-STAT1* (cyto)",
        "STAT1c_STAT1c_star" => "STAT1-STAT1* (cyto)",
        "STAT1c_star_PPX" => "STAT1*-PPX (cyto)",
        "STAT1n" => "STAT1 (nuc)",
        "STAT1n_star" => "STAT1* (nuc)",
        "STAT1n_star_STAT1n_star" => "STAT1*-STAT1* (nuc)",
        "STAT1n_STAT1n_star" => "STAT1-STAT1* (nuc)",
        "STAT1n_star_PPN" => "STAT1*-PPN (nuc)",
        "SHP2" => "SHP2",
        "PPX" => "PPX",
        "PPN" => "PPN",
        "SOCS1" => "SOCS1",
        "mRNAn" => "mRNA (nuc)",
        "mRNAc" => "mRNA (cyto)",
        "IRF1" => "IRF1",
        "IFNRJ2_star_STAT1c" => "IFN-R-JAK₂*-STAT1",
        "IFNRJ2_star_SHP2" => "IFN-R-JAK₂*-SHP2",
        "IFNRJ2_star_SOCS1" => "IFN-R-JAK₂*-SOCS1",
        "IFNRJ2_star_SOCS1_STAT1c" => "IFN-R-JAK₂*-SOCS1-STAT1",
        "STAT1c_star_STAT1c_star_PPX" => "STAT1*-STAT1*-PPX",
        "STAT1n_star_STAT1n_star_PPN" => "STAT1*-STAT1*-PPN",
        "IFNRJ2_star_SHP2_STAT1c" => "IFN-R-JAK₂*-SHP2-STAT1",
        "IFNRJ2_star_SHP2_SOCS1" => "IFN-R-JAK₂*-SHP2-SOCS1"
    )
    
    # Plot each group in a separate figure with subplots
    for (group_name, group_species) in species_groups
        # Find indices of species in this group
        species_indices = Int[]
        species_display_names = String[]
        
        for sp in group_species
            for (idx, name) in enumerate(species_names)
                if occursin(sp, name)
                    push!(species_indices, idx)
                    
                    # Find the matching display name
                    display_name = get(display_names, sp, replace(sp, "_" => "-"))
                    push!(species_display_names, display_name)
                    break
                end
            end
        end
        
        if isempty(species_indices)
            continue
        end
        
        # Determine layout
        n_species = length(species_indices)
        n_cols = min(3, n_species)
        n_rows = ceil(Int, n_species / n_cols)
        
        plt = plot(layout=(n_rows, n_cols), size=(300*n_cols, 250*n_rows),
                   margin=5Plots.mm)
        
        for (i, (idx, display_name)) in enumerate(zip(species_indices, species_display_names))
            plot!(plt[i], time_in_hours, averaged_species[:, idx],
                 ribbon=std_dev_species[:, idx],
                 fillalpha=0.2,
                 lw=2,
                 xlabel="Time (hours)",
                 ylabel="Molecules",
                 title=display_name,
                 label=nothing,
                 linecolor=:blue,
                 fillcolor=:blue)
        end
        
        # Add overall title
        # title!(plt, group_name, titlefontsize=16)
        
        # Save figure with clean filename
        clean_filename = replace(group_name, " " => "_")
        display(plt)
        savefig(plt, joinpath(output_dir, "$(clean_filename).png"))
    end
    
    # Create an overview plot of key species in the pathway
    key_species = ["IFN", "STAT1n_star_STAT1n_star", "mRNAc", "IRF1"] 
    key_indices = Int[]
    key_display_names = String[]
    
    for sp in key_species
        for (idx, name) in enumerate(species_names)
            if occursin(sp, name)
                push!(key_indices, idx)
                # Find the matching display name
                display_name = get(display_names, sp, replace(sp, "_" => "-"))
                push!(key_display_names, display_name)
                break
            end
        end
    end

    print(key_display_names)
    
    plt_key = plot(layout=(2,2), size=(800, 600), margin=5Plots.mm)
    
    for (i, (idx, display_name)) in enumerate(zip(key_indices, key_display_names))
        plot!(plt_key[i], time_in_hours, averaged_species[:, idx],
             ribbon=std_dev_species[:, idx],
             fillalpha=0.2,
             lw=3,
             xlabel="Time (hours)",
             ylabel="Molecules",
             title=display_name,
             label=nothing,
             linecolor=:blue,
             fillcolor=:blue)
    end
    
    # title!(plt_key, "Key JAK-STAT Pathway Components", titlefontsize=16)
    display(plt_key)
    savefig(plt_key, joinpath(output_dir, "Key_Pathway_Components.png"))
    
    # Measure memory after simulation
    GC.gc()
    mem_after = Sys.total_memory() - Sys.free_memory()
    mem_used = (mem_after - mem_before) / (1024^2)
    println("Memory used: $(round(mem_used, digits=2)) MB")
    
    # # Return the aggregated data
    # return Dict(
    #     "species_names" => species_names,
    #     "time_points" => time_points,
    #     "time_in_hours" => time_in_hours,
    #     "averaged_species" => averaged_species,
    #     "std_dev_species" => std_dev_species,
    #     "display_names" => display_names
    # );
end

# Run the function
all_species_data = @time simulate_and_plot_all_species()
