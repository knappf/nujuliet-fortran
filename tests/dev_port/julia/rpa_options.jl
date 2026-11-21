# dev-port test case 'rpa_options': TDA/RPA options not covered by the other cases, bundled interaction (hw = 16,
# Nmax = 3, N2max = 6, N3max = 9).
#
#   io/O16/        16O HF -> TDA Ortho = false (Name TDA-noOrtho), RPA Ortho = false (Name RPA-noOrtho),
#                  RPA TBCF = true, RPA_Orb = true, Transitions = 'Exact' (Name RPA-TBCF: RPA_sigma incl. fix F2, density
#                  energy, natural-orbital export, exact transitions with the RPA Sigma)
#   io/O16_MBPT/   16O HF + MBPT(Export, TBCF, ScH) -> RPA Reference = 'HF-MBPT', Reference_TBCF = true, TBCF = true,
#                  Transitions = 'Exact' (Name RPA-dress: Sigma_AB dressing energy), and the same with DEV = true
#                  (Name RPA-DEV-TBCF: RPA_allocate_SE path)
#   stage 3 part 2 (external routines & the DEV path):
#   io/O16/        externals on TDA-noOrtho / RPA-noOrtho: analysis of a non-EJ multipole (2-, Cut 0.05), 2+ transition
#                  densities (File_Name D2), 0+ / 3- operator densities, 2+ currents with M = 1 on a 20 x 30 map,
#                  TDA/RPA_transitions_recalculate; RPA Transitions = 'Exact' (Name RPA-Exact) + 1- densities, analysis,
#                  currents & recalculation (exact Linear-free kernel)
#   io/O16_MBPT/   RPA Reference_TBCF, TBCF, DEV, 'Exact' (Name RPA-DEV-dress: RPA_allocate_SE in A & in the dressing energy)
#                  + all RPA externals with the RPA sigma.bin; TDA Reference_TBCF (Name TDA-MBPT) + operator densities,
#                  analysis & recalculation with the reference Sigma
#                  RPA Reference_TBCF, DEV, Tol_par = 0.9, Tol_hol = 0.02 (Name RPA-DEV-frac: the natural orbitals 0p1/2,
#                  1s1/2 & 0d are particle & hole orbitals, backward coupling B_12 of RPA_allocate_SE)
#   io/C12/        12C HF -> RPA TBCF = true (Name RPA-TBCF) and HF_ScRPA (Name HF-ScRPA, final Sigma -> E_corr density; with TBCF = true
#                  12C diverges in dev, the initial RPA is unstable): the particle
#                  0p1/2 precedes the hole 0p3/2 in the orbital list, so the h-p-p-h terms of RPA_sigma fall into the computed
#                  Ket <= Bra triangle in the pp, nn & pn branches (fix F2: Y in every hpph term, Sigma symmetric)
#
# Run from the work directory: julia -t N,1 rpa_options.jl  (see ../run_julia.sh)

include(joinpath(@__DIR__, "common.jl"))

const IntP = BUNDLED_INT
calc(; kw...) = Calculation_Parameters(; A = 16, Z = 8, hw = 16.0, Nmax = 3, N2max = 6, N3max = 9, kw...)

let P = "O16"
    io = "io/" * P
    Params = Parameters(IntP, calc(Path = P,
        TDA = TDA_Parameters(Name = "TDA-noOrtho", Ortho = false),
        RPA = RPA_Parameters(Name = "RPA-noOrtho", Ortho = false)))
    step(() -> HF(Params), "O16 HF"; expect = [io * "/HF/bin/orb.bin"])
    step(() -> TDA(Params), "O16 TDA Ortho=false"; expect = [io * "/TDA-noOrtho/bin"])
    step(() -> RPA(Params), "O16 RPA Ortho=false"; expect = [io * "/RPA-noOrtho/bin"])

    ParamsT = Parameters(IntP, calc(Path = P,
        RPA = RPA_Parameters(Name = "RPA-TBCF", TBCF = true, RPA_Orb = true, Transitions = "Exact")))
    step(() -> RPA(ParamsT), "O16 RPA TBCF RPA_Orb Exact"; expect = [io * "/RPA-TBCF/bin/sigma.bin", io * "/RPA-TBCF/bin/H1b.bin"])

    step("O16 TDA/RPA-noOrtho externals") do
        TDA_transition_analysis(Params; JP = "2-", nu_list = [1, 2], Cut = 0.05)
        RPA_transition_analysis(Params; JP = "2-", nu_list = [1, 2], Cut = 0.05)
        TDA_transition_densities(Params; JP = "2+", nu_list = [1, 3], File_Name = "D2")
        RPA_transition_densities(Params; JP = "2+", nu_list = [1, 3], File_Name = "D2")
        TDA_transition_operator_densities(Params; JP = "0+", nu_list = [1, 2])
        RPA_transition_operator_densities(Params; JP = "3-", nu_list = [1, 2])
        TDA_transition_currents(Params; JP = "2+", M = 1, nu_list = [1, 2], xN_grid = 20, zN_grid = 30)
        RPA_transition_currents(Params; JP = "2+", M = 1, nu_list = [1, 2], xN_grid = 20, zN_grid = 30)
        TDA_transitions_recalculate(Params)
        RPA_transitions_recalculate(Params)
    end

    ParamsX = Parameters(IntP, calc(Path = P, RPA = RPA_Parameters(Name = "RPA-Exact", Transitions = "Exact")))
    step(() -> RPA(ParamsX), "O16 RPA Transitions=Exact"; expect = [io * "/RPA-Exact/bin"])
    step("O16 RPA-Exact externals") do
        RPA_transition_densities(ParamsX; JP = "1-", nu_list = [5, 6, 7])
        RPA_transition_analysis(ParamsX; JP = "1-", nu_list = [5, 6, 7])
        RPA_transition_currents(ParamsX; JP = "1-", M = 0, nu_list = [5, 6, 7])
        RPA_transitions_recalculate(ParamsX)
    end
end

let P = "O16_MBPT"
    io = "io/" * P
    HFP = HF_Parameters(MBPT = HF_MBPT_Parameters(Export = true, TBCF = true, ScH = true))
    Params = Parameters(IntP, calc(Path = P, HF = HFP,
        RPA = RPA_Parameters(Name = "RPA-dress", Reference = "HF-MBPT", Reference_TBCF = true, TBCF = true,
                             Transitions = "Exact")))
    step(() -> HF(Params), "O16_MBPT HF+MBPT export"; expect = [io * "/HF-MBPT/bin/sigma.bin"])
    step(() -> RPA(Params), "O16_MBPT RPA dress"; expect = [io * "/RPA-dress/bin/sigma.bin"])

    ParamsD = Parameters(IntP, calc(Path = P, HF = HFP,
        RPA = RPA_Parameters(Name = "RPA-DEV-TBCF", Reference = "HF-MBPT", Reference_TBCF = true, DEV = true)))
    step(() -> RPA(ParamsD), "O16_MBPT RPA DEV Reference_TBCF"; expect = [io * "/RPA-DEV-TBCF/bin"])

    ParamsE = Parameters(IntP, calc(Path = P, HF = HFP,
        TDA = TDA_Parameters(Name = "TDA-MBPT", Reference = "HF-MBPT", Reference_TBCF = true),
        RPA = RPA_Parameters(Name = "RPA-DEV-dress", Reference = "HF-MBPT", Reference_TBCF = true, TBCF = true, DEV = true,
                             Transitions = "Exact")))
    step(() -> RPA(ParamsE), "O16_MBPT RPA DEV dress"; expect = [io * "/RPA-DEV-dress/bin/sigma.bin"])
    step(() -> TDA(ParamsE), "O16_MBPT TDA Reference_TBCF"; expect = [io * "/TDA-MBPT/bin"])
    ParamsF = Parameters(IntP, calc(Path = P, HF = HFP,
        RPA = RPA_Parameters(Name = "RPA-DEV-frac", Reference = "HF-MBPT", Reference_TBCF = true, DEV = true,
                             Tol_par = 0.9, Tol_hol = 0.02)))
    step(() -> RPA(ParamsF), "O16_MBPT RPA DEV overlapping particle/hole sets"; expect = [io * "/RPA-DEV-frac/bin"])

    step("O16_MBPT RPA-DEV-dress & TDA-MBPT externals") do
        RPA_transition_densities(ParamsE; JP = "1-", nu_list = [5, 6, 7])
        RPA_transition_operator_densities(ParamsE; JP = "1-", nu_list = [5, 6, 7])
        RPA_transition_analysis(ParamsE; JP = "1-", nu_list = [5, 6, 7])
        RPA_transition_currents(ParamsE; JP = "1-", M = 0, nu_list = [5, 6, 7])
        RPA_transitions_recalculate(ParamsE)
        TDA_transition_operator_densities(ParamsE; JP = "1-", nu_list = [5, 6, 7])
        TDA_transition_analysis(ParamsE; JP = "1-", nu_list = [5, 6, 7])
        TDA_transitions_recalculate(ParamsE)
    end
end

let P = "C12"
    io = "io/" * P
    calc12(; kw...) = Calculation_Parameters(; A = 12, Z = 6, hw = 16.0, Nmax = 3, N2max = 6, N3max = 9, Path = P, kw...)
    Params = Parameters(IntP, calc12(RPA = RPA_Parameters(Name = "RPA-TBCF", TBCF = true)))
    step(() -> HF(Params), "C12 HF"; expect = [io * "/HF/bin/orb.bin"])
    step(() -> RPA(Params), "C12 RPA TBCF"; expect = [io * "/RPA-TBCF/bin/sigma.bin"])
    ParamsS = Parameters(IntP, calc12(RPA = RPA_Parameters(Reference = "HF", Name = "HF-ScRPA", Reference_TBCF = false,
        Ortho = true, Transitions = "Linear", TBCF = false, Tol_hol = 0.3, Tol_par = 0.85, Tol = 1e-5, Imax = 100,
        OBDM_diagonal = false, ScH = true)))
    step(() -> HF_ScRPA(ParamsS), "C12 HF_ScRPA"; expect = [io * "/HF-ScRPA/RPA_summary.dat"])
end

finish("rpa_options")
