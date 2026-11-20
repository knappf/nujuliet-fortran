# dev-port test case 'bundled': the small interaction bundled with dev (io/NN.bin, io/NNN.bin; hw = 16, Nmax = 3,
# N2max = 6, N3max = 9). Mirrors dev scripts/sample_script_run.jl and adds option variants.
#
#   io/O16/          16O  HF -> TDA -> RPA (+ transition densities/currents/operator densities/analysis of both),
#                         RPA Transitions='Exact' (Name RPA-Exact), RPA DEV=true (Name RPA-DEV)
#   io/O16_MBPT/     16O  HF with MBPT(Export, TBCF, ScH) -> TDA, RPA with Reference='HF-MBPT', Reference_TBCF=true
#                         (+ transition densities/currents of both)
#   io/He4/          4He  HF -> HF_ScRPA (sample parameters), HF_ScRPA TBCF=true (Name HF-ScRPA-TBCF),
#                         HF_ScRPA OBDM_diagonal=true (Name HF-ScRPA-OBDMdiag)
#   io/O18_BCS/      18O  BCS -> QTDA -> QRPA (Reference 'BCS')
#   io/O18_HFB/      18O  HFB -> QTDA -> pnQTDA -> QRPA
#   io/parse/        V2b_parse / V3b_no2b_parse to Nmax 2 / N2max 4 / N3max 6
#
# Run from the work directory: julia -t N,1 bundled.jl  (see ../run_julia.sh)

include(joinpath(@__DIR__, "common.jl"))

const IntP = BUNDLED_INT
calc(; kw...) = Calculation_Parameters(; A = 16, Z = 8, hw = 16.0, Nmax = 3, N2max = 6, N3max = 9, kw...)

const JP1 = "1-"
const NU1 = [5, 6, 7]

# ---------------------------------------------------------------------------------------------------------------------
# 16O: HF -> TDA -> RPA (sample script) ...
let P = "O16"
    io = "io/" * P
    Params = Parameters(IntP, calc(Path = P))
    step(() -> HF(Params), "O16 HF"; expect = [io * "/HF/HF_summary.dat", io * "/HF/bin/orb.bin"])
    step(() -> TDA(Params), "O16 TDA"; expect = [io * "/TDA/bin"])
    step(() -> RPA(Params), "O16 RPA"; expect = [io * "/RPA/bin"])

    step("O16 TDA/RPA transition densities") do
        RPA_transition_densities(Params; JP = JP1, nu_list = NU1, File_Name = "Phonon_Densities")
        TDA_transition_densities(Params; JP = JP1, nu_list = NU1, File_Name = "Phonon_Densities")
    end
    step("O16 TDA/RPA transition currents") do
        TDA_transition_currents(Params; JP = JP1, M = 0, nu_list = NU1, File_Name = "Phonon_Currents",
                                xN_grid = 40, zN_grid = 40)
        RPA_transition_currents(Params; JP = JP1, M = 0, nu_list = NU1, File_Name = "Phonon_Currents",
                                xN_grid = 40, zN_grid = 40)
    end
    step("O16 TDA/RPA transition operator densities") do
        TDA_transition_operator_densities(Params; JP = JP1, nu_list = NU1)
        RPA_transition_operator_densities(Params; JP = JP1, nu_list = NU1)
    end
    step("O16 TDA/RPA transition analysis") do
        TDA_transition_analysis(Params; JP = JP1, nu_list = NU1)
        RPA_transition_analysis(Params; JP = JP1, nu_list = NU1)
    end

    # RPA with exact transitions ...
    ParamsX = Parameters(IntP, calc(Path = P, RPA = RPA_Parameters(Name = "RPA-Exact", Transitions = "Exact")))
    step(() -> RPA(ParamsX), "O16 RPA Transitions=Exact"; expect = [io * "/RPA-Exact/bin"])

    # Experimental RPA DEV path (RPA_allocate_SE) ...
    ParamsD = Parameters(IntP, calc(Path = P, RPA = RPA_Parameters(Name = "RPA-DEV", DEV = true)))
    step(() -> RPA(ParamsD), "O16 RPA DEV=true"; expect = [io * "/RPA-DEV/bin"])
end

# ---------------------------------------------------------------------------------------------------------------------
# 16O: HF + MBPT(Export, TBCF, ScH) -> TDA / RPA on the correlated HF-MBPT reference with its TBCF ...
let P = "O16_MBPT"
    io = "io/" * P
    Params = Parameters(IntP, calc(Path = P,
        HF  = HF_Parameters(MBPT = HF_MBPT_Parameters(Export = true, TBCF = true, ScH = true)),
        TDA = TDA_Parameters(Name = "TDA-MBPT", Reference = "HF-MBPT", Reference_TBCF = true),
        RPA = RPA_Parameters(Name = "RPA-MBPT", Reference = "HF-MBPT", Reference_TBCF = true)))
    step(() -> HF(Params), "O16_MBPT HF+MBPT export"; expect = [io * "/HF/bin/orb.bin", io * "/HF-MBPT/bin/sigma.bin"])
    step(() -> TDA(Params), "O16_MBPT TDA"; expect = [io * "/TDA-MBPT/bin"])
    step(() -> RPA(Params), "O16_MBPT RPA"; expect = [io * "/RPA-MBPT/bin"])
    step("O16_MBPT TDA/RPA transition densities+currents") do
        TDA_transition_densities(Params; JP = JP1, nu_list = NU1)
        RPA_transition_densities(Params; JP = JP1, nu_list = NU1)
        TDA_transition_currents(Params; JP = JP1, M = 0, nu_list = NU1, xN_grid = 40, zN_grid = 40)
        RPA_transition_currents(Params; JP = JP1, M = 0, nu_list = NU1, xN_grid = 40, zN_grid = 40)
    end
end

# ---------------------------------------------------------------------------------------------------------------------
# 4He: HF -> HF_ScRPA (sample script) + TBCF and OBDM_diagonal variants ...
let P = "He4"
    io = "io/" * P
    scrpa(; kw...) = RPA_Parameters(; Reference = "HF", Name = "HF-ScRPA", Reference_TBCF = false, Ortho = true,
                                    Transitions = "Linear", TBCF = false, Tol_hol = 0.3, Tol_par = 0.85, Tol = 1e-5,
                                    Imax = 100, OBDM_diagonal = false, ScH = true, kw...)
    Params = Parameters(IntP, calc(A = 4, Z = 2, Path = P, RPA = scrpa()))
    step(() -> HF(Params), "He4 HF"; expect = [io * "/HF/bin/orb.bin"])
    step(() -> HF_ScRPA(Params), "He4 HF_ScRPA"; expect = [io * "/HF-ScRPA/bin"])

    ParamsT = Parameters(IntP, calc(A = 4, Z = 2, Path = P, RPA = scrpa(Name = "HF-ScRPA-TBCF", TBCF = true)))
    step(() -> HF_ScRPA(ParamsT), "He4 HF_ScRPA TBCF=true"; expect = [io * "/HF-ScRPA-TBCF/bin"])

    ParamsO = Parameters(IntP, calc(A = 4, Z = 2, Path = P, RPA = scrpa(Name = "HF-ScRPA-OBDMdiag", OBDM_diagonal = true)))
    step(() -> HF_ScRPA(ParamsO), "He4 HF_ScRPA OBDM_diagonal=true"; expect = [io * "/HF-ScRPA-OBDMdiag/bin"])
end

# ---------------------------------------------------------------------------------------------------------------------
# 18O: BCS -> QTDA -> QRPA ...
let P = "O18_BCS"
    io = "io/" * P
    Params = Parameters(IntP, calc(A = 18, Z = 8, Path = P,
        QTDA = QTDA_Parameters(Ortho = true, Reference = "BCS"),
        QRPA = QRPA_Parameters(Ortho = true, Reference = "BCS")))
    step(() -> BCS(Params), "O18_BCS BCS"; expect = [io * "/BCS"])
    step(() -> QTDA(Params), "O18_BCS QTDA"; expect = [io * "/QTDA"])
    step(() -> QRPA(Params), "O18_BCS QRPA"; expect = [io * "/QRPA"])
end

# ---------------------------------------------------------------------------------------------------------------------
# 18O: HFB -> QTDA -> pnQTDA -> QRPA ...
let P = "O18_HFB"
    io = "io/" * P
    Params = Parameters(IntP, calc(A = 18, Z = 8, Path = P,
        QTDA = QTDA_Parameters(Ortho = true, Reference = "HFB"),
        QRPA = QRPA_Parameters(Ortho = true, Reference = "HFB")))
    step(() -> HFB(Params), "O18_HFB HFB"; expect = [io * "/HFB"])
    step(() -> QTDA(Params), "O18_HFB QTDA"; expect = [io * "/QTDA"])
    step(() -> pnQTDA(Params), "O18_HFB pnQTDA"; expect = [io * "/pnQTDA"])
    step(() -> QRPA(Params), "O18_HFB QRPA"; expect = [io * "/QRPA"])
end

# ---------------------------------------------------------------------------------------------------------------------
# Interaction parsing (sample script) ...
let
    isdir("io/parse") || mkdir("io/parse")
    Params = Parameters(IntP, calc())
    step("parse V2b/V3b_no2b") do
        V2b_parse(Params, 2, 4, "io/parse/NN_parsed.bin")
        V3b_no2b_parse(Params, 2, 4, 6, "io/parse/NNN_parsed.bin")
    end
end

finish("bundled")
