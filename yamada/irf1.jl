using CSV
using DataFrames
using Catalyst
using DifferentialEquations
using JumpProcesses
using Plots
using Interpolations
using Statistics
using ProgressMeter
using Random

"""
Load IRF1 trajectories data from CSV file
"""
function load_irf1_data(file_path="yamada/output/IRF1_trajectories.csv")
    df = CSV.read(file_path, DataFrame)
    
    # Extract time and mean IRF1 values
    times = df.Time_hours
    irf1_values = df.Mean
    
    println("Loaded IRF1 data with $(length(times)) time points")
    plot(times, irf1_values, label="IRF1 data")
    return times, irf1_values
end

