# dev-port test case 'scrpa_options': HF-ScRPA options of no other case (stage 4), bundled interaction (hw = 16, Nmax = 3,
# N2max = 6, N3max = 9). The He4 steps use the sample-script parameters of 'bundled' (Tol_hol = 0.3, Tol_par = 0.85,
# Imax = 100) plus:
#
#   io/He4/ScRPA-noScH         ScH = false (reference H1b/H2b rotated to the NAT basis; Fortran skips the bare reads):
#                              the outer iteration oscillates and stops after Imax = 100 (non-convergence summary)
#   io/He4/ScRPA-noScH-TBCF    ScH = false, TBCF = true, TBCF_alpha = 0.5, TBCF_beta = 3.0 (Sigma only in A, B & energy)
#   io/He4/ScRPA-ramp          TBCF = true, TBCF_alpha = 1.0 (no Sigma mixing), TBCF_beta = 2.0, OBDM_alpha = 0.8,
#                              RPA_Orb = true, Transitions = 'Exact'
#   io/He4/ScRPA-DEV           TBCF = true, DEV = true (RPA_allocate_SE in A and in the dressing energy)
#   io/He4/ScRPA-noconv        TBCF = true, Imax = 3: non-convergence (termination summary RPA_summary.dat only)
#   io/He4/HF-ScRPA            default sample run + RPA externals on the ScRPA output directory
#   io/O16_MBPT/ScRPA-MBPT     HF + MBPT(Export, TBCF, ScH) -> HF_ScRPA Reference = 'HF-MBPT' (fractional reference
#                              occupations), default RPA_Parameters tolerances
#   io/O20/ScRPA-frac          20O HF ESH = false (0d5/2 neutrons 4/6) -> HF_ScRPA Tol_par = 0.9, Tol_hol = 0.5, Imax = 30
#
# Run from the work directory: julia -t N,1 scrpa_options.jl  (see ../run_julia.sh)

include(joinpath(@__DIR__, "common.jl"))

const IntP = BUNDLED_INT
calc(; kw...) = Calculation_Parameters(; A = 4, Z = 2, hw = 16.0, Nmax = 3, N2max = 6, N3max = 9, Path = "He4", kw...)
scrpa(; kw...) = RPA_Parameters(; Reference = "HF", Name = "HF-ScRPA", Reference_TBCF = false, Ortho = true,
                                Transitions = "Linear", TBCF = false, Tol_hol = 0.3, Tol_par = 0.85, Tol = 1e-5,
                                Imax = 100, OBDM_diagonal = false, ScH = true, kw...)
he4(; kw...) = Parameters(IntP, calc(RPA = scrpa(; kw...)))

let io = "io/He4"
    step(() -> HF(he4()), "He4 HF"; expect = [io * "/HF/bin/orb.bin"])
    step(() -> HF_ScRPA(he4(Name = "ScRPA-noScH", ScH = false)), "He4 HF_ScRPA ScH=false";
         expect = [io * "/ScRPA-noScH/RPA_summary.dat"])
    step(() -> HF_ScRPA(he4(Name = "ScRPA-noScH-TBCF", ScH = false, TBCF = true, TBCF_alpha = 0.5, TBCF_beta = 3.0)),
         "He4 HF_ScRPA ScH=false TBCF"; expect = [io * "/ScRPA-noScH-TBCF/bin/sigma.bin"])
    step(() -> HF_ScRPA(he4(Name = "ScRPA-ramp", TBCF = true, TBCF_alpha = 1.0, TBCF_beta = 2.0, OBDM_alpha = 0.8,
                            RPA_Orb = true, Transitions = "Exact")),
         "He4 HF_ScRPA ramp OBDM_alpha RPA_Orb Exact"; expect = [io * "/ScRPA-ramp/bin/H1b.bin"])
    step(() -> HF_ScRPA(he4(Name = "ScRPA-DEV", TBCF = true, DEV = true)), "He4 HF_ScRPA TBCF DEV";
         expect = [io * "/ScRPA-DEV/bin/sigma.bin"])
    step(() -> HF_ScRPA(he4(Name = "ScRPA-noconv", TBCF = true, Imax = 3)), "He4 HF_ScRPA non-convergence";
         expect = [io * "/ScRPA-noconv/RPA_summary.dat"])

    P = he4()
    step(() -> HF_ScRPA(P), "He4 HF_ScRPA"; expect = [io * "/HF-ScRPA/bin"])
    step("He4 HF-ScRPA RPA externals") do
        RPA_transition_densities(P; JP = "1-", nu_list = [1, 2])
        RPA_transition_operator_densities(P; JP = "2+", nu_list = [1, 2])
        RPA_transition_analysis(P; JP = "1-", nu_list = [1, 2], Cut = 0.05)
        RPA_transition_currents(P; JP = "1-", M = 0, nu_list = [1, 2])
        RPA_transitions_recalculate(P)
    end
end

let io = "io/O16_MBPT"
    HFP = HF_Parameters(MBPT = HF_MBPT_Parameters(Export = true, TBCF = true, ScH = true))
    P = Parameters(IntP, calc(A = 16, Z = 8, Path = "O16_MBPT", HF = HFP,
                              RPA = RPA_Parameters(Name = "ScRPA-MBPT", Reference = "HF-MBPT")))
    step(() -> HF(P), "O16_MBPT HF+MBPT export"; expect = [io * "/HF-MBPT/bin/H2b.bin"])
    step(() -> HF_ScRPA(P), "O16_MBPT HF_ScRPA Reference=HF-MBPT"; expect = [io * "/ScRPA-MBPT/RPA_summary.dat"])
end

let io = "io/O20"
    P = Parameters(IntP, calc(A = 20, Z = 8, Path = "O20", HF = HF_Parameters(ESH = false),
                              RPA = RPA_Parameters(Name = "ScRPA-frac", Tol_par = 0.9, Tol_hol = 0.5, Imax = 30)))
    step(() -> HF(P), "O20 HF (ESH=false)"; expect = [io * "/HF/bin/orb.bin"])
    step(() -> HF_ScRPA(P), "O20 HF_ScRPA fractional"; expect = [io * "/ScRPA-frac/RPA_summary.dat"])
end

finish("scrpa_options")
