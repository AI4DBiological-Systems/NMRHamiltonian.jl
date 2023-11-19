# run a.jl first.

include("./helpers/utils.jl")

PLT.close("all")
fig_num = 1

#PLT.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

#T = Float64
T = Float32
#T = Float16 # will result in numerical precision-related errors. Please use Float32 or above.

### user inputs.

molecule_entries = [
    "L-Serine";
    "alpha-D-Glucose";
    "beta-D-Glucose";
    "Ethanol";
    "L-Methionine";     
    "L-Phenylalanine";
    "L-Glutathione reduced";
    "L-Glutathione oxidized";       
    "Uridine";
    "L-Glutamine";
    "L-Valine";
    "DSS";
]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs, SW, ν_0ppm = HAM.getpresetspectrometer(T, "700")

## pull the sample coupling values into dictionary structures.

root_data_path = DS.getdatapath(DS.NMR2023()) # coupling values data repository root path

H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

molecule_mapping_root_path = joinpath(
    root_data_path,
    "molecule_name_mapping",
)
molecule_mapping_file_path = joinpath(
    molecule_mapping_root_path,
    "select_molecules.json",
)

config = HAM.SHConfig{T}(
    coherence_tol = convert(T, 0.01),
    relative_α_threshold = convert(T, 0.005),
    tol_radius_1D = convert(T, 0.1), # strictly between 0 and 1. The lower, the better the approximation, but would a larger partition (i.e. more resonance groups).
    nuc_factor = convert(T, 1.5),
)

### end inputs.



println("Timing: getphysicalparameters")
@time Phys = NMRHamiltonian.getphysicalparameters(
    T,
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path;
    unique_cs_atol = convert(T, 1e-6),
)

println("Timing: simulate()")
@time As, MSPs = HAM.simulate(
    Phys,
    molecule_entries,
    fs,
    SW,
    ν_0ppm,
    config,
)
println()

cs_values = collect( HAM.readbasechemshifts(Phys[n]) for n in eachindex(Phys) )
println("The reference chemical shift values (in ppm) for the Δc features: [compoune][spin system][ME nuclei]")
display(cs_values)
println()

del_c_bars = collect( As[n].Δc_bar for n in eachindex(As) )
N_res_groups = collect( length.(As[n].Δc_bar) for n in eachindex(As) )
N_nucs = collect( 
    collect(
        length(As[n].Δc[i][begin])
        for i in eachindex(As[n].Δc)
    ) 
    for n in eachindex(As)
)
println("[molecule_entries N_nucs N_res_groups]:")
display([molecule_entries N_nucs N_res_groups])
println()

# typically, the number of nuclei in a spin system should less than or equal to the number of resonance groups for that spin system.
@show N_nucs .<= N_res_groups # DSS is an exception it seems.

### visualize a target molecule and spin group, given a T-2 decay parameter.
λ0 = convert(T, 3.4) # T-2 decay parameter.
molecule_select = 1
spin_system_select = 1
ppm_offset = convert(T, 0.2) # for display padding.
N_viz = 50000



a = As[molecule_select].αs[spin_system_select]
F = As[molecule_select].Ωs[spin_system_select]

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

ΩS_ppm = hz2ppmfunc.( F ./ HAM.twopi(T) )
ΩS_ppm_sorted = sort(ΩS_ppm)

u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - ppm_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + ppm_offset)

P_min = hz2ppmfunc(u_min)
P_max = hz2ppmfunc(u_max)

P = LinRange{T}(P_min, P_max, N_viz)
U = ppm2hzfunc.(P)
U_rad = U .* HAM.twopi(T)

# absorption Lorentzian.

function evalzerophasecl1Darray(u_rad, αs::Vector{T}, Ωs::Vector{T}, λ::T)::Complex{T} where T <: AbstractFloat

    out = zero(Complex{T})
    for l = 1:length(αs)
        out += evalzerophaseclpartitionelement(u_rad, αs[l], Ωs[l], λ)
    end

    return out
end

function evalzerophaseclpartitionelement(r,
    α::T, Ω::T, λ::T)::Complex{T} where T <: AbstractFloat

    out = α/(λ+im*(r-Ω))

    return out
end

q = uu->evalzerophasecl1Darray(uu, a, F, λ0)
q_U = q.(U_rad)

# plot.
PLT.figure(fig_num)
fig_num += 1

PLT.plot(P, real.(q_U))

PLT.gca().invert_xaxis()
PLT.ylabel("real part")
PLT.title("spectrum of $(molecule_entries[molecule_select]), spin system $(spin_system_select)")

###### round-trip tests.

include("./helpers/roundtrip.jl")
roundtripJSON(As)
roundtripBSON(As)

roundtripJSON(Phys, molecule_entries)
roundtripBSON(Phys, molecule_entries)

# save to disk.
Phys_filename = "Phys.json"
As_filename = "As.json"

S_Phys = HAM.serializephysicalparams(Phys, molecule_entries)
HAM.saveasJSON(
    Phys_filename,
    S_Phys,
)

S_As = HAM.serializemixture(As)
HAM.saveasJSON(
    As_filename,
    S_As,
)

# load from disk.
dict_As = HAM.readJSON(As_filename)
As_rec = HAM.deserializemixture(dict_As)

dict_Phys = HAM.readJSON(Phys_filename)
Phys_rec, molecule_entries_rec = HAM.deserializephysicalparams(Dict(dict_Phys))

@assert typeof(As_rec) == typeof(As)
@assert typeof(Phys_rec) == typeof(Phys)

println("The following should be close to zero.")
@show sum( norm(As_rec[n].Δc - As[n].Δc) for n in eachindex(As) ) 

# clean-update
rm(As_filename)
rm(Phys_filename)

nothing