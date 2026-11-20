# Stage-4 unit check of the HF-ScRPA building blocks (no Julia source changes: the patched reference functions are
# called directly). Run from a work dir that contains io/He4/HF (the HF of bundled step 07), see run_unit.sh.
# Dumps (raw Float64/Int64, column-major) to dump/:
#   XY.bin    initial RPA amplitudes per (J,P): Int64 N, Re X, Im X, Re Y, Im Y, Re E, Im E
#   oY4.bin   RPA_OBDM_oY4(Rho_ref,X,Y) (p then n, a_max x a_max); oY4c.bin the same with X, Y times e^(i(0.37 nu + J))
#   nat.bin   C of OBDM_natural_orbitals(oY4) (p, n) ; tr.bin transform_amplitudes(X,Y,C) per (J,P): Re X, Im X, Re Y, Im Y
#   solve.bin HF_ScRPA_OBDM_solve(Orb_ref,...,X,Y,Rho_ref,I): Rho (p,n), C_LHO_NAT (p,n), C_HF_NAT (p,n)
#   H_*       RPA_Hamiltonian for ScH/TBCF variants (O1b_export / O2b_export files)
#   energy.txt E_QBA of the initial RPA and RPA_energy (Density = true, Sigma = RPA_sigma) after one ScRPA step
const NUJULIET_REF = get(ENV, "NUJULIET_REF", "/scratch/local.ssd/knapp/job_25503048.pbs-m1/NuJuliet_dev_ref")
include(joinpath(NUJULIET_REF, "src", "import_NuJuliet.jl"))
using .NuJuliet
using LinearAlgebra, Printf
const NJ = NuJuliet
BLAS.set_num_threads(Threads.nthreads())
NJ.wigner_init_float(75, "Jmax", 9)

IntP = NJ.Interaction_Parameters(NN_File = joinpath(NUJULIET_REF, "io", "NN.bin"), NNN_File = joinpath(NUJULIET_REF, "io", "NNN.bin"),
                                 hw = 16.0, Nmax = 3, N2max = 6, N3max = 9)
mk(; kw...) = NJ.Parameters(IntP, NJ.Calculation_Parameters(; A = 4, Z = 2, hw = 16.0, Nmax = 3, N2max = 6, N3max = 9, Path = "He4",
              RPA = NJ.RPA_Parameters(; Reference = "HF", Name = "HF-ScRPA", Reference_TBCF = false, Ortho = true, Transitions = "Linear",
                                      TBCF = false, Tol_hol = 0.3, Tol_par = 0.85, Tol = 1e-5, Imax = 100, OBDM_diagonal = false,
                                      ScH = true, kw...)))
Params = mk()
mkpath("dump")
IO = "io/He4/HF"
Orb_ref = NJ.orbitals_import(Params, IO * "/bin/orb.bin")
T = NJ.T1b(Params, Orb_ref)
V_NN_bare, Orb_NN = NJ.V2b_read(Params)
V_NNN_bare, Orb_NNN = NJ.V3b_no2b_read(Params)
C_LHO_HF = NJ.O1b_import(Params, Orb_ref, IO * "/bin/C.bin")
Rho_ref = NJ.O1b_import(Params, Orb_ref, IO * "/bin/rho.bin")
h_N = NJ.O1b_import(Params, Orb_ref, IO * "/bin/H1b.bin")
V_NN = NJ.O2b_import(Params, Orb_ref, IO * "/bin/H2b.bin")
Sigma0 = NJ.O2b_initialize_empty()
Orb_1p1h = NJ.orbitals_1p1h_make(Params, Orb_ref; Tol_par = 0.85, Tol_hol = 0.3)
Orb_ph = NJ.RPA_OBDM_orbitals_make(Params, Orb_ref)
TrOp_LHO = NJ.Tr1b_initialize(Params, Orb_ref, 1.0)
TrOp = NJ.Tr1b_transformation(Params, Orb_ref, C_LHO_HF, TrOp_LHO)
A, B, N = NJ.RPA_allocate(Params, Orb_ref, Orb_1p1h, Orb_NN, h_N, V_NN, Rho_ref, Sigma0, false)
E, X, Y, Stab = NJ.RPA_diagonalize(Params, Orb_ref, Orb_1p1h, A, B, N, TrOp)
E_QBA0 = NJ.RPA_energy_QBA(Params, Orb_ref, Orb_1p1h, E, X, Y)

wo(io, O) = (write(io, O.p); write(io, O.n))
open("dump/XY.bin", "w") do io
    for J in 0:size(Orb_1p1h.N,2)-1, P in 1:2
        n = Orb_1p1h.N[P,J+1]
        write(io, Int64(n))
        n == 0 && continue
        write(io, real.(X[P,J+1])); write(io, imag.(X[P,J+1])); write(io, real.(Y[P,J+1])); write(io, imag.(Y[P,J+1]))
        write(io, real.(E[P,J+1])); write(io, imag.(E[P,J+1]))
    end
end
pM = Matrix{Matrix{ComplexF64}}(undef, size(X)); nM = Matrix{Matrix{ComplexF64}}(undef, size(X))
Rho1 = NJ.RPA_OBDM_oY4(Params, Orb_ref, Orb_ph, Orb_1p1h, Rho_ref, X, Y, pM, nM)
open(io -> wo(io, Rho1), "dump/oY4.bin", "w")
# complex amplitudes with generic column phases (dev's M = -1/2 X^T W X^* is not phase invariant: exact formula check) ...
Xc = copy(X); Yc = copy(Y)
for J in 0:size(Orb_1p1h.N,2)-1, P in 1:2
    n = Orb_1p1h.N[P,J+1]; n == 0 && continue
    ph = [exp(im * (0.37 * k + J)) for k in 1:n]
    Xc[P,J+1] = X[P,J+1] .* transpose(ph); Yc[P,J+1] = Y[P,J+1] .* transpose(ph)
end
Rho1c = NJ.RPA_OBDM_oY4(Params, Orb_ref, Orb_ph, Orb_1p1h, Rho_ref, Xc, Yc, pM, nM)
open(io -> wo(io, Rho1c), "dump/oY4c.bin", "w")
_, Cn, _ = NJ.OBDM_natural_orbitals(Params, Orb_ref, Rho1)
open(io -> wo(io, Cn), "dump/nat.bin", "w")
X2, Y2 = NJ.RPA_OBDM_transform_amplitudes(Orb_1p1h, X, Y, Cn)
open("dump/tr.bin", "w") do io
    for J in 0:size(Orb_1p1h.N,2)-1, P in 1:2
        Orb_1p1h.N[P,J+1] == 0 && continue
        write(io, real.(X2[P,J+1])); write(io, imag.(X2[P,J+1])); write(io, real.(Y2[P,J+1])); write(io, imag.(Y2[P,J+1]))
    end
end

function solve_dump(Pr, tag)
    C_HF_NAT = NJ.O1B(diagm(ones(size(Rho_ref.p,1))), diagm(ones(size(Rho_ref.p,1))))
    Rho, C_LHO_NAT, Orb = NJ.HF_ScRPA_OBDM_solve(Pr, Orb_ref, Orb_ph, Orb_1p1h, X, Y, Rho_ref, C_HF_NAT)
    open(io -> (wo(io, Rho); wo(io, C_LHO_NAT); wo(io, C_HF_NAT)), "dump/solve" * tag * ".bin", "w")
    return Rho, C_LHO_NAT, Orb, C_HF_NAT
end
Rho, C_LHO_NAT, Orb, C_HF_NAT = solve_dump(Params, "")
solve_dump(mk(OBDM_diagonal = true), "_diag")
solve_dump(mk(OBDM_alpha = 1.0), "_a1")

Sigma = NJ.RPA_sigma(Params, Orb, Orb_1p1h, Orb_NN, X, Y)       # NAT occupations, reference-basis amplitudes (test input only)
NJ.O2b_export(Params, Orb_NN, Sigma, "dump/Sigma_in.bin")
for (tag, Pr) in (("ScH", Params), ("ScH_TBCF", mk(TBCF = true)), ("HFref", mk(ScH = false)))
    hN, VNN = NJ.RPA_Hamiltonian(Pr, Orb, Orb_NN, Orb_NNN, C_LHO_NAT, Rho, deepcopy(Sigma), T, V_NN_bare, V_NNN_bare)
    NJ.O1b_export(Pr, Orb, hN, "dump/H1b_" * tag * ".bin")
    NJ.O2b_export(Pr, Orb_NN, VNN, "dump/H2b_" * tag * ".bin")
end

# One full ScRPA step (ScH, TBCF = false) & the final energy evaluation ...
hN, VNN = NJ.RPA_Hamiltonian(Params, Orb, Orb_NN, Orb_NNN, C_LHO_NAT, Rho, Sigma0, T, V_NN_bare, V_NNN_bare)
TrOp1 = NJ.Tr1b_transformation(Params, Orb, C_LHO_NAT, TrOp_LHO)
A, B, N = NJ.RPA_allocate(Params, Orb, Orb_1p1h, Orb_NN, hN, VNN, Rho, Sigma0, false)
E1, X1, Y1, Stab1 = NJ.RPA_diagonalize(Params, Orb, Orb_1p1h, A, B, N, TrOp1)
E_QBA1 = NJ.RPA_energy_QBA(Params, Orb, Orb_1p1h, E1, X1, Y1)
Sig1 = NJ.RPA_sigma(Params, Orb, Orb_1p1h, Orb_NN, X1, Y1)
Ec = NJ.RPA_energy(Params, Orb, Orb_1p1h, Orb_NN, E1, X1, Y1, Rho, NJ.O1b_transformation(Rho_ref, C_HF_NAT), hN, Sig1, VNN; Density = true)
open("dump/energy.txt", "w") do io
    println(io, @sprintf("%.17e %.17e %.17e %.17e %.17e %.17e %.17e", E_QBA0, E_QBA1, Ec.E_QBA, Ec.E_sigma, Ec.E_rho, Ec.E_dress, Ec.E_dens))
end
println("UNIT_DONE")
