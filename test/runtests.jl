
import NMRHamiltonian
const HAM = NMRHamiltonian

using DataDeps
using Tar, CodecZlib

using LinearAlgebra

using Test

import JSON3
import BSON

include("../examples/helpers/data.jl") # use DataDeps.jl and Tar.jl to download and extract the sample coupling values.
include("../examples/helpers/utils.jl")
include("../examples/helpers/roundtrip.jl")

# for the simulation tests, only tests whether it runs.
@testset "simulate, round-trip, MSPs" begin
    # Write your tests here.
    
    root_data_path = getdatapath() # coupling values data repository root path

    H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

    molecule_mapping_root_path = joinpath(
        root_data_path,
        "molecule_name_mapping",
    )


    mapping_path_set = Vector{String}(undef, 2)
    mapping_path_set[1] = joinpath(
        molecule_mapping_root_path,
        "select_molecules.json",
    )
    mapping_path_set[2] = joinpath(
        molecule_mapping_root_path,
        "select_molecules_multiple_ME.json",
    )

    entries_set = Vector{Vector{String}}(undef, 2)
    entries_set[1] = [
        "L-Valine";
        # #"L-Leucine"; # long. # all trues.
        # #"L-Isoleucine"; # might be long.
        "Epinephrine";
        "Dopamine";
        "Gamma-Aminobutyric acid";
        "L-Cysteine";
        "L-Aspartic acid";
        "L-Asparagine";
        "L-Proline";
        "Purine";
        "Serotonin";
        "Putrescine";
        "beta-Alanine";
        "L-Alanine";
        "L-(+) Lactic acid";
        "L-Glutamine";
        "Ethanol";
        "L-Methionine";     
        "L-Phenylalanine";
        "L-Glutathione reduced";
        "L-Glutathione oxidized";  
        "L-Histidine";    
        "alpha-D-Glucose";
        "beta-D-Glucose";
        "DSS";
        "Singlet - 4.7 ppm";
    ]
    entries_set[2] = [
        "N1-(1,3-thiazolan-2-yliden)-2-ethoxyaniline Maybridge_Ro3_Fragment_12_C09 simulation_1";
        "2-Piperazin-1-ylaniline Maybridge_Ro3_Fragment_01_E09 simulation_1";
    ]

    spectrometer_tags = ["400"; "500"; "600"; "700"; "900"]

    Ts = [Float32; Float64]
    
    for path_ind in eachindex(mapping_path_set)

        molecule_mapping_file_path = mapping_path_set[path_ind]
        molecule_entries = entries_set[path_ind]

        for spectrometer_tag in spectrometer_tags
            for T in Ts
                @show molecule_mapping_file_path, molecule_entries, T, spectrometer_tag
                
                coherence_sum_zero_tol::T = convert(T, 1e-6)
                roundtrip_zero_tol::T = convert(T, 1e-6)
                if T == Float64
                    roundtrip_zero_tol = convert(T, 1e-14)
                    coherence_sum_zero_tol = convert(T, 1e-14)
                end

                # the amount from -1 can sum(a Δc) and each entry of a Δc vector be, whilst not being flagged as an invalid Δc vector.
                Δc_valid_atol::T = convert(T, 0.1)


                config = HAM.SHConfig{T}(
                    coherence_tol = convert(T, 0.01),
                    relative_α_threshold = convert(T, 0.005),
                    tol_radius_1D = convert(T, 0.1), # strictly between 0 and 1. The lower, the better the approximation, but would a larger partition (i.e. more resonance groups).
                    nuc_factor = convert(T, 1.5),
                )
                unique_cs_atol = convert(T, 1e-6)
                
                Phys, As, MSPs = HAM.loadandsimulate(
                    T,
                    spectrometer_tag,
                    molecule_entries,
                    H_params_path,
                    molecule_mapping_file_path;
                    config = config,
                    unique_cs_atol = unique_cs_atol
                )

                # tests
                roundtripJSON(As)
                roundtripBSON(As)

                roundtripJSON(Phys, molecule_entries)
                roundtripBSON(Phys, molecule_entries)

                # round-trip test from cs_sys to ζs.
                for n in eachindex(Phys)

                    test = deepcopy(Phys[n])

                    for i in eachindex(Phys[n].cs_sys)
                        ζs = NMRHamiltonian.readbasechemshifts(Phys[n])
                        cs = similar(Phys[n].cs_sys)
                
                        fill!(test.cs_sys[i], zero(T))
                        NMRHamiltonian.writebasechemshifts!(test, ζs)

                        @test norm(test.cs_sys[i]- Phys[n].cs_sys[i]) < roundtrip_zero_tol
                        # end round-trip test.
                    end
                end
            end
        end
    end

end

