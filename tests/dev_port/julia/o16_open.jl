# dev-port test case 'o16_open': open-shell / fractional-occupation case on the bundled interaction (hw = 16,
# Nmax = 3, N2max = 6, N3max = 9). 20O = 16O core + 4 neutrons in 0d5/2 (occupation 4/6), HF with ESH = false.
#
#   io/O20/         HF -> TDA -> RPA with the default Tol_par = 0.1 / Tol_hol = 0.7: the 2/3-filled 0d5/2 lies
#                   between the two thresholds and is dropped from the 1p1h space (reported upstream behaviour)
#   io/O20/TDA-frac, RPA-frac   TDA / RPA with Tol_par = 0.9, Tol_hol = 0.5: the 0d5/2 is both a particle and a hole
#                   orbital, so the fractional occupations enter the 1p1h matrices and transitions
#
# Run from the work directory: julia -t N,1 o16_open.jl  (see ../run_julia.sh)

include(joinpath(@__DIR__, "common.jl"))

const IntP = BUNDLED_INT
calc(; kw...) = Calculation_Parameters(; A = 20, Z = 8, hw = 16.0, Nmax = 3, N2max = 6, N3max = 9, Path = "O20",
                                       HF = HF_Parameters(ESH = false), kw...)

let io = "io/O20"
    Params = Parameters(IntP, calc())
    step(() -> HF(Params), "O20 HF (ESH=false)"; expect = [io * "/HF/HF_summary.dat", io * "/HF/bin/orb.bin"])
    step(() -> TDA(Params), "O20 TDA"; expect = [io * "/TDA/bin"])
    step(() -> RPA(Params), "O20 RPA"; expect = [io * "/RPA/bin"])

    ParamsF = Parameters(IntP, calc(TDA = TDA_Parameters(Name = "TDA-frac", Tol_par = 0.9, Tol_hol = 0.5),
                                    RPA = RPA_Parameters(Name = "RPA-frac", Tol_par = 0.9, Tol_hol = 0.5)))
    step(() -> TDA(ParamsF), "O20 TDA fractional 1p1h"; expect = [io * "/TDA-frac/bin"])
    step(() -> RPA(ParamsF), "O20 RPA fractional 1p1h"; expect = [io * "/RPA-frac/bin"])
    step("O20 TDA/RPA-frac transition densities") do
        TDA_transition_densities(ParamsF; JP = "2+", nu_list = [1, 2])
        RPA_transition_densities(ParamsF; JP = "2+", nu_list = [1, 2])
    end
end

finish("o16_open")
