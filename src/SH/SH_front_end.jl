# """
# fetchSHparameters(molecule_name::String, record_name::String, base_path::String)
#
# Get the chemical shift and J-coupling parameters from the locally-stored GISSMO record.
# """
# function fetchSHparameters(molecule_name::String, record_name::String,
#     base_path::String)
#
#     file_name = "$(molecule_name)_$(record_name)"
#
#     css_sys, J_IDs_sys, J_vals_sys, cs_singlets,
#     J_lb_sys, J_ub_sys, cs_lb_sys, cs_ub_sys,
#     cs_LUT, reference_concentration, solute_concentration,
#     _ = GISSMOReader.loadGISSMOmolecule(base_path, file_name)
#
#     #N_groups = length(cs_group)
#     p_cs_sys = cstopcs(css_sys, cs_LUT)
#
#     N_spins_singlet = Vector{Int}(undef, length(cs_singlets))
#     cs_singlets_compact = Vector{Float64}(undef, length(cs_singlets))
#     for m = 1:length(cs_singlets)
#         cs_singlets_compact[m] = cs_singlets[m][1]
#         N_spins_singlet[m] = length(cs_singlets[m])
#     end
#
#     N_sys = length(css_sys)
#     @assert length(p_cs_sys) == N_sys
#
#     N_spins_sys = collect( length(css_sys[m]) for m = 1:length(css_sys) )
#     intermediates_sys = prepcouplingalgorithm(N_spins_sys)
#
#     # prepare initial guess and bounds.
#     #tmp = collect( [cs_group[i]; (J_lb_group[i] + J_ub_group[i]) ./ 2] for i = 1:N_groups )
#     css = combinevectors(p_cs_sys)
#     push!(css, cs_singlets_compact...)
#
#     N_p_cs_sys = collect( length(p_cs_sys[i]) for i = 1:N_sys )
#     p0 = copy(css)
#
#     # could be used for dim-reduced p.
#     p_compound = Vector{Vector{Float64}}(undef, N_sys + length(cs_singlets_compact))
#     for i = 1:N_sys
#         #p_compound[i] = p_cs_sys[i]
#         p_compound[i] = copy(p_cs_sys[i])
#     end
#     for i = 1:length(cs_singlets_compact)
#         p_compound[N_sys+i] = [ cs_singlets_compact[i] ]
#     end
#
#     # the full cs version of p_compound.
#     cs_compound = Vector{Vector{Float64}}(undef, N_sys + length(cs_singlets))
#     for i = 1:N_sys
#         cs_compound[i] = copy(p_cs_sys[i])
#     end
#     for i = 1:length(cs_singlets)
#         cs_compound[N_sys+i] = copy(cs_singlets[i])
#     end
#
#     return p_cs_sys, css_sys, cs_singlets, J_vals_sys, J_IDs_sys, intermediates_sys,
#     cs_LUT, cs_singlets_compact, N_spins_singlet, css, p_compound, cs_compound
# end
