# dev-port test case 'ca40_5sh': real interaction DN2LOGO394 (NN + NNN, hw = 12 MeV, file truncation emax = 16,
# e2max = 32, e3max = 28), 40Ca in 5 oscillator shells (Nmax = 4, N2max = 8, N3max = 8).
#
#   io/Ca40/   HF + MBPT(Export) -> TDA -> RPA (+ transition densities of the lowest 1- states)
#              -> HF_ScRPA (Name HF-ScRPA, sample-script parameters)
#
# Run from the work directory: julia -t N,1 ca40_5sh.jl  (see ../run_julia.sh)

include(joinpath(@__DIR__, "common.jl"))

const IntP = DN2LOGO394_INT(12)
calc(; kw...) = Calculation_Parameters(; A = 40, Z = 20, hw = 12.0, Nmax = 4, N2max = 8, N3max = 8, Path = "Ca40", kw...)

let io = "io/Ca40"
    Params = Parameters(IntP, calc(HF = HF_Parameters(MBPT = HF_MBPT_Parameters(Export = true))))
    step(() -> HF(Params), "Ca40 HF+MBPT export"; expect = [io * "/HF/HF_summary.dat", io * "/HF/bin/orb.bin",
                                                           io * "/HF-MBPT/bin"])
    step(() -> TDA(Params), "Ca40 TDA"; expect = [io * "/TDA/bin"])
    step(() -> RPA(Params), "Ca40 RPA"; expect = [io * "/RPA/bin"])
    step("Ca40 TDA/RPA transition densities") do
        TDA_transition_densities(Params; JP = "1-", nu_list = [1, 2, 3])
        RPA_transition_densities(Params; JP = "1-", nu_list = [1, 2, 3])
    end

    ParamsS = Parameters(IntP, calc(RPA = RPA_Parameters(Reference = "HF", Name = "HF-ScRPA", Ortho = true,
                         Transitions = "Linear", TBCF = false, Tol_hol = 0.3, Tol_par = 0.85, Tol = 1e-5, Imax = 100,
                         OBDM_diagonal = false, ScH = true)))
    step(() -> HF_ScRPA(ParamsS), "Ca40 HF_ScRPA"; expect = [io * "/HF-ScRPA/bin"])
end

finish("ca40_5sh")
