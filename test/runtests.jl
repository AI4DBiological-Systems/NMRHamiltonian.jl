
import NMRHamiltonian

using DataDeps
using Tar, CodecZlib

using LinearAlgebra

using Test

import JSON3
import BSON

include("../examples/helpers/data.jl") # use DataDeps.jl and Tar.jl to download and extract the sample coupling values.
include("../examples/helpers/utils.jl")

include("./helpers/setup.jl")
include("./helpers/roundtrip.jl")


@testset "serialization round-trip" begin
    # Write your tests here.
    
    T = Float64
    coherence_sum_zero_tol = 1e-14
    roundtrip_zero_tol = 1e-14

    # the amount from -1 can sum(a Δc) and each entry of a Δc vector be, whilst not being flagged as an invalid Δc vector.
    Δc_valid_atol = 0.1
    
    root_data_path = getdatapath() # coupling values data repository root path
    H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.
    molecule_mapping_root_path = joinpath(root_data_path, "molecule_name_mapping")
    #molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "GISSMO_names.json")

    ### common molecules.
    molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "select_molecules.json")

    molecule_entries = [
    #"alpha-D-Glucose";
    # "HEPES";
    "L-Valine";
    # #"L-Leucine"; # long. # all trues.
    # #"L-Isoleucine"; # might be long.
    # "ATP";
    # "Epinephrine";
    # "Dopamine";
    # "Gamma-Aminobutyric acid";
    # "L-Cysteine";
    # "L-Aspartic acid";
    # "L-Asparagine";
    # "L-Proline";
    # "Purine";
    # "Serotonin";
    # "Putrescine";
    # "beta-Alanine";
    # "L-Alanine";
    # "L-(+) Lactic acid";
    # "L-Glutamine";
    # "Ethanol";
    # "L-Methionine";     
    # "L-Phenylalanine";
    # "L-Glutathione reduced";
    # "L-Glutathione oxidized";  
    # "L-Histidine";     
    "DSS";
    "Singlet - 4.7 ppm";
]

    Phys, As, Rs, MSPs = setupSH(molecule_entries, H_params_path, molecule_mapping_file_path)
    
    # about the number of resonance groups
    status, n, i = NMRHamiltonian.checkNgroups(As; verbose = true)

    println("Does all spin systems each have a larger number of resonance groups than its nuclei (with magnetic equivalence taken into account)?")
    println(status)
    if !status
        println("An example spin system where this does not hold is (compound index, system index) = ($n,$i).")
    end

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

    ###########################################

    ### multiple spin systems that have magnetic equivalence.
    molecule_mapping_file_path = joinpath(
        molecule_mapping_root_path, "select_molecules_multiple_ME.json")
    molecule_entries = [
        "N1-(1,3-thiazolan-2-yliden)-2-ethoxyaniline Maybridge_Ro3_Fragment_12_C09 simulation_1";
        "2-Piperazin-1-ylaniline Maybridge_Ro3_Fragment_01_E09 simulation_1";
    ]
    
    Phys, As, Rs, MSPs = setupSH(molecule_entries, H_params_path, molecule_mapping_file_path)
    
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
    
    #coherence_diagnostics = NMRHamiltonian.getΔcdiagnostics(As)

    AS_Δc_valid, coherence_diagnostics = NMRHamiltonian.checkcoherences(As; atol = Δc_valid_atol)
    @test AS_Δc_valid

    # test the partial contributions and full quantum numbers agree.
    MSPs_valid = NMRHamiltonian.checkcoherences(MSPs; coherence_sum_zero_tol = coherence_sum_zero_tol)
    @test MSPs_valid

    # about the number of resonance groups
    status, n, i = NMRHamiltonian.checkNgroups(As; verbose = true)
    
    println("Does all spin systems each have a larger number of resonance groups than its nuclei (with magnetic equivalence taken into account)?")
    println(status)
    if !status
        println("An example spin system where this does not hold is (compound index, system index) = ($n,$i).")
    end

end
