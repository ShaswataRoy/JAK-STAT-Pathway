using CSV
using DataFrames
using Statistics
using Plots
using ProgressMeter
using StatsBase

"""
Read all FISH quantification data from a folder and store in a DataFrame.
This program:
1. Recursively finds all FISH quant files in the specified folder
2. Reads each file and extracts relevant data
3. Combines all data into a single DataFrame with metadata
4. Performs basic statistical analysis
5. Generates summary visualizations
"""

# Configuration
const FISH_QUANT_DIR = joinpath(@__DIR__, "fish_quant")  # Change this to your actual directory
const OUTPUT_DIR = joinpath(@__DIR__, "fish_quant_analysis")
const FILE_EXTENSIONS = [".csv", ".txt"]  # Only use CSV and text files

# Create output directory if it doesn't exist
mkpath(OUTPUT_DIR)

function find_fish_quant_files(dir_path)
    all_files = String[]
    
    println("Searching for FISH quant files in: $dir_path")
    
    # Walk through directory recursively
    for (root, dirs, files) in walkdir(dir_path)
        for file in files
            # Check if file has one of our target extensions
            if any(endswith(lowercase(file), ext) for ext in FILE_EXTENSIONS)
                push!(all_files, joinpath(root, file))
            end
        end
    end
    
    println("Found $(length(all_files)) potential FISH quant files.")
    return all_files
end

function extract_metadata_from_path(file_path)
    # Extract useful metadata from the file path
    # This function should be customized based on your specific file naming convention
    parts = splitpath(file_path)
    
    # Default metadata
    metadata = Dict(
        "file_name" => basename(file_path),
        "full_path" => file_path
    )
    
    # Extract folder names as potential metadata
    for (i, part) in enumerate(parts)
        if i < length(parts)  # Not the filename
            # Try to identify if this folder name represents metadata
            if occursin("sample", lowercase(part))
                metadata["sample"] = part
            elseif occursin("time", lowercase(part)) || occursin("hour", lowercase(part))
                metadata["time_point"] = part
            elseif occursin("treatment", lowercase(part)) || occursin("condition", lowercase(part))
                metadata["condition"] = part
            elseif occursin("replicate", lowercase(part)) || occursin("rep", lowercase(part))
                metadata["replicate"] = part
            end
        end
    end

    return metadata
end

function read_fish_quant_file(file_path)
    metadata = extract_metadata_from_path(file_path)
    
    try
        # Try to read as CSV with various possible delimiters
        for delimiter in [',', '\t', ';', ' ']
            try
                data = CSV.read(file_path, DataFrame, delim=delimiter)
                
                # If we get here, reading worked
                metadata["format"] = "delimited"
                metadata["delimiter"] = string(delimiter)
                
                # Check if we got actual data (more than 0 rows and at least 2 columns)
                if nrow(data) > 0 && ncol(data) >= 2
                    return data, metadata
                end
            catch
                # Try the next delimiter
                continue
            end
        end
        
        # If all delimiter attempts failed, try more flexible reading
        data = CSV.read(file_path, DataFrame, delim=nothing)
        metadata["format"] = "auto-detected"
        return data, metadata
        
    catch e
        @warn "Error reading file $file_path: $e"
        return nothing, metadata
    end
end

function process_fish_quant_data(data, metadata)
    # If data reading failed, return empty DataFrame with metadata
    if data === nothing
        processed = DataFrame()
        # Add metadata as a single row
        for (key, value) in metadata
            processed[!, "meta_$key"] = [value]
        end
        return processed
    end
    
    # Initialize processed DataFrame by copying the original data
    processed = copy(data)
    
    # Rename columns based on their content
    for col in names(processed)
        col_lower = lowercase(col)
        if contains(col_lower, r"cell|id")
            rename!(processed, col => "cell_id")
        elseif contains(col_lower, r"count|spot|intens|signal")
            rename!(processed, col => "value_$col")
        elseif contains(col_lower, r"pos|coord|x|y|z")
            rename!(processed, col => "pos_$col")
        else
            rename!(processed, col => "data_$col")
        end
    end
    
    # Add metadata columns - ensure we create arrays of the correct length
    n_rows = nrow(processed)
    for (key, value) in metadata
        processed[!, "meta_$key"] = fill(value, n_rows)
    end
    
    return processed
end


function (
    # Main function to run the entire workflow
    
    # 1. Find all FISH quant files
    quant_files = find_fish_quant_files(FISH_QUANT_DIR)
    
    if isempty(quant_files)
        println("No FISH quant files found. Check the directory path and file extensions.")
        return
    end
    
    # 2. Read and process each file
    all_data = DataFrame()
    
    println("Processing FISH quant files...")
    p = Progress(length(quant_files), desc="Reading files: ", dt=0.5)
    
    for file_path in quant_files
        # Read the file
        data, metadata = read_fish_quant_file(file_path)
        
        # Process the data
        processed_data = process_fish_quant_data(data, metadata)
        
        # Add to the combined DataFrame
        if nrow(processed_data) > 0
            if nrow(all_data) == 0
                all_data = processed_data
            else
                # Ensure all columns exist in both DataFrames
                for col in setdiff(names(processed_data), names(all_data))
                    all_data[!, col] = missing
                end
                for col in setdiff(names(all_data), names(processed_data))
                    processed_data[!, col] = missing
                end
                
                append!(all_data, processed_data)
            end
        end
        
        next!(p)
    end
    
    # 3. Save the combined data
    println("Saving combined data...")
    CSV.write(joinpath(OUTPUT_DIR, "all_fish_quant_data.csv"), all_data)
    
    println("\nSummary:")
    println("- Processed $(length(quant_files)) FISH quant files")
    println("- Combined data has $(nrow(all_data)) rows and $(ncol(all_data)) columns")
    println("- Results saved to $OUTPUT_DIR")
end
