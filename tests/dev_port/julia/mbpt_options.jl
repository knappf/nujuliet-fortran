# dev-port test case 'mbpt_options': HF-MBPT export options not covered by the other cases, bundled interaction
# (hw = 16, Nmax = 3, N2max = 6, N3max = 9), HF + MBPT only.
#
#   io/O28_TBCF/  28O (N = 20 closed sd shell, N != Z) MBPT(Export, TBCF, ScH = false): sigma.bin with different proton &
#                 neutron particle/hole patterns (fix F1: nn channel classified with nO) & the ScH = false Hamiltonian
#   io/O16_ScH/   16O MBPT(Export, ScH, TBCF = false): polarized H1b/H2b without the correlation function
#
# Run from the work directory: julia -t N,1 mbpt_options.jl  (see ../run_julia.sh)

include(joinpath(@__DIR__, "common.jl"))

const IntP = BUNDLED_INT
calc(; kw...) = Calculation_Parameters(; hw = 16.0, Nmax = 3, N2max = 6, N3max = 9, kw...)

let io = "io/O28_TBCF"
    Params = Parameters(IntP, calc(A = 28, Z = 8, Path = "O28_TBCF",
                                   HF = HF_Parameters(MBPT = HF_MBPT_Parameters(Export = true, TBCF = true, ScH = false))))
    step(() -> HF(Params), "O28 HF+MBPT TBCF"; expect = [io * "/HF-MBPT/bin/sigma.bin"])
end

let io = "io/O16_ScH"
    Params = Parameters(IntP, calc(A = 16, Z = 8, Path = "O16_ScH",
                                   HF = HF_Parameters(MBPT = HF_MBPT_Parameters(Export = true, TBCF = false, ScH = true))))
    step(() -> HF(Params), "O16 HF+MBPT ScH"; expect = [io * "/HF-MBPT/bin/H2b.bin"])
end

finish("mbpt_options")
