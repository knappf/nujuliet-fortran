# dev-port test case 'qp_options': quasiparticle-chain options not covered by 'bundled' (18O BCS & HFB with default
# options), bundled interaction (hw = 16, Nmax = 3, N2max = 6, N3max = 9).
#
#   io/Ne22_HFB/   22Ne (open proton & neutron shells) HFB (Pairing "Full") -> QTDA (Ortho = false) -> pnQTDA
#                  (QRPA of 22Ne has complex (A-B)(A+B) eigenvectors: dev stops with an InexactError, the port with an
#                  explicit error, so it is not part of the case)
#   io/O20_LN/     20O HFB with Lipkin-Nogami (LNT = true, LNRes = true; Broyden restart at iteration 17) -> QRPA
#   io/O20_MCA/    20O HFB Pairing = "MCA" -> QTDA
#   io/O20_HBCS/   20O HFB Pairing = "BCS" (HFB_diagonalize_BCS with HF_orbitals_make) -> pnQTDA
#   io/O20_ScBCS/  20O BCS with ScBCS = true and user reference A0 = 16, Z0 = 8 -> QTDA, pnQTDA, QRPA (Ortho = false),
#                  all with Reference = "BCS"
#
# Run from the work directory: julia -t N,1 qp_options.jl  (see ../run_julia.sh)

include(joinpath(@__DIR__, "common.jl"))

const IntP = BUNDLED_INT
calc(; kw...) = Calculation_Parameters(; hw = 16.0, Nmax = 3, N2max = 6, N3max = 9, kw...)

let P = "Ne22_HFB"
    io = "io/" * P
    Params = Parameters(IntP, calc(A = 22, Z = 10, Path = P,
        QTDA = QTDA_Parameters(Ortho = false, Reference = "HFB"),
        pnQTDA = pnQTDA_Parameters(Reference = "HFB")))
    step(() -> HFB(Params), "Ne22 HFB"; expect = [io * "/HFB/bin/qpH2b.bin"])
    step(() -> QTDA(Params), "Ne22 QTDA Ortho=false"; expect = [io * "/QTDA/bin/QTDA_X.bin"])
    step(() -> pnQTDA(Params), "Ne22 pnQTDA"; expect = [io * "/pnQTDA/bin/pnQTDA_X.bin"])
end

let P = "O20_LN"
    io = "io/" * P
    Params = Parameters(IntP, calc(A = 20, Z = 8, Path = P, HFB = HFB_Parameters(LNT = true, LNRes = true)))
    step(() -> HFB(Params), "O20 HFB LNT"; expect = [io * "/HFB/bin/qpH2b.bin"])
    step(() -> QRPA(Params), "O20 LNT QRPA"; expect = [io * "/QRPA/bin/QRPA_X.bin"])
end

let P = "O20_MCA"
    io = "io/" * P
    Params = Parameters(IntP, calc(A = 20, Z = 8, Path = P, HFB = HFB_Parameters(Pairing = "MCA")))
    step(() -> HFB(Params), "O20 HFB MCA"; expect = [io * "/HFB/bin/qpH2b.bin"])
    step(() -> QTDA(Params), "O20 MCA QTDA"; expect = [io * "/QTDA/bin/QTDA_X.bin"])
end

let P = "O20_HBCS"
    io = "io/" * P
    Params = Parameters(IntP, calc(A = 20, Z = 8, Path = P, HFB = HFB_Parameters(Pairing = "BCS")))
    step(() -> HFB(Params), "O20 HFB Pairing=BCS"; expect = [io * "/HFB/bin/qpH2b.bin"])
    step(() -> pnQTDA(Params), "O20 HBCS pnQTDA"; expect = [io * "/pnQTDA/bin/pnQTDA_X.bin"])
end

let P = "O20_ScBCS"
    io = "io/" * P
    Params = Parameters(IntP, calc(A = 20, Z = 8, Path = P,
        BCS = BCS_Parameters(ScBCS = true, A0 = 16, Z0 = 8),
        QTDA = QTDA_Parameters(Ortho = true, Reference = "BCS"),
        QRPA = QRPA_Parameters(Ortho = false, Reference = "BCS"),
        pnQTDA = pnQTDA_Parameters(Reference = "BCS")))
    step(() -> BCS(Params), "O20 ScBCS"; expect = [io * "/BCS/bin/qpH2b.bin"])
    step(() -> QTDA(Params), "O20 ScBCS QTDA"; expect = [io * "/QTDA/bin/QTDA_X.bin"])
    step(() -> pnQTDA(Params), "O20 ScBCS pnQTDA"; expect = [io * "/pnQTDA/bin/pnQTDA_X.bin"])
    step(() -> QRPA(Params), "O20 ScBCS QRPA Ortho=false"; expect = [io * "/QRPA/bin/QRPA_X.bin"])
end

finish("qp_options")
