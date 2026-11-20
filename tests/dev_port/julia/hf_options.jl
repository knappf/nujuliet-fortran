# dev-port test case 'hf_options': HF options not covered by the other cases, bundled interaction (hw = 16, Nmax = 3,
# N2max = 6, N3max = 9), HF only (BMF = false).
#
#   io/O16_VHF/   16O HF with V_HF = true (non-local radial HF potential, CMS1+2B)
#   io/O18_ESH/   18O HF with ESH = true (enforced subshell closure: the 2 valence neutrons skip the 0d5/2 and fill
#                 the 1s1/2) and V_HF = true with CMS = "None" (full kinetic term subtracted; CMS2B crashes in Julia,
#                 porting note 2). (20O with ESH = true stops in both codes with the dev "closed-subshell" error:
#                 the 0f7/2 slot a_max lies below the remaining eligible candidates, see UPSTREAM_BUGS.md.)
#
# Run from the work directory: julia -t N,1 hf_options.jl  (see ../run_julia.sh)

include(joinpath(@__DIR__, "common.jl"))

const IntP = BUNDLED_INT
calc(; kw...) = Calculation_Parameters(; hw = 16.0, Nmax = 3, N2max = 6, N3max = 9, kw...)

let io = "io/O16_VHF"
    Params = Parameters(IntP, calc(A = 16, Z = 8, Path = "O16_VHF", HF = HF_Parameters(V_HF = true, BMF = false)))
    step(() -> HF(Params), "O16 HF V_HF"; expect = [io * "/HF/densities/HF_radial_potential.dat"])
end

let io = "io/O18_ESH"
    Params = Parameters(IntP, calc(A = 18, Z = 8, Path = "O18_ESH", CMS = "None",
                                   HF = HF_Parameters(ESH = true, V_HF = true, BMF = false)))
    step(() -> HF(Params), "O18 HF ESH V_HF CMS None"; expect = [io * "/HF/densities/HF_radial_potential.dat"])
end

finish("hf_options")
