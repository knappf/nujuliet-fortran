# Dev-port spec, group "common B": interactions, transitions, Hamiltonian

Base: NuJuliet 2595e82 (`Source/MatrixElements/{V2B,V3B_NO2B,Tr1B}.jl`), dev: f3f7a33 (`src/common/{V2b,V3b_no2b,Tr1b,Hamiltonian}.jl`).
Fortran: `src/core/{nj_v2b,nj_v3b,nj_tr1b}.f90` (+ new `src/core/nj_hamiltonian.f90`).

Note on diffs: `Tr1B.jl -> Tr1b.jl` shows as a full rewrite only because of CRLF line endings. Use
`diff --strip-trailing-cr -w` to see the real changes. `T1b.jl` is unchanged. `T2b.jl` only renames `f6j` to `w6j`.
`cg/w3j/w6j/w9j` in dev `functions.jl:16-30` are plain wrappers of `fCG/f3j/f6j/f9j`, so `nj_wigner` stays as it is.
`hat(j) = j + 1` (`functions.jl:40`, NOT sqrt) and `parity(L) = rem(L,2)+1` (`functions.jl:9`).
`Interaction_Parameters` is unchanged (only CRLF differs).

---------------------------------------------------------------------------------------------------------------

## 1. V2b.jl (634 lines) <- base V2B.jl (459) ; Fortran `nj_v2b.f90`

### 1.1 `V2b_read`: signature, orbitals, buffering
* Dev signature `V2b_read(Params) -> (V_NN, Orb_NN)` (was `V2b_read(Params, Orb)`). The orbitals are now made inside:
  `Orb = orbitals_make(Params,Int = true)` (V2b.jl:16). Dev `orbitals_make` defaults to **Calc.Nmax** and uses
  Int.Nmax only with `Int=true` (orb.jl:1-5). Base `orbitals_make` always used Int.Nmax. The `orbitals_make(...,Int)`
  change belongs to group common A, but this group depends on it.
* `T = T1b(Params,Orb)` and `O2b_initialize(Params,Orb,Make_Orb_NN=true)` now get the **Int** orbitals. O2b_initialize
  sizes everything by Calc.Nmax, so `Orb_NN` is the same as before.
* The file is now read in blocks of at most `Buffer_size = 1000000` Float64 (V2b.jl:29, 70-82) instead of `read(io,Float64)`
  per element. The order and decoding (`ME_00, ME_nn, ME_10, ME_pp` per J) are unchanged. `V2b_count` now returns
  `(N_Chunk, N_Chunk_skip)` (V2b.jl:304). Element counts stay the same, and the skip is in bytes (`8*sum N_Chunk[b<a]`).
* The accumulation formulas, CMS1+2B/CMS2B, and the `pG2N/nG2N` contact enhancement are **unchanged**. This includes the base bug
  `T[a,c]` on an `O1B` (V2b.jl:178, 201-203), so port-note deviation 2 is still needed.
* The commented-out `V2b_0p_analysis` export block was deleted. Nothing to port.

Fortran (`V2b_read`):
* New interface `subroutine V2b_read(Params,V_NN,Orb_NN)`, with `Orb_Int = orbitals_make(Params, Int=.true.)` made
  internally (see shared interface S1). Update the callers `nj_hf_driver`, `nj_bcs`, `nj_hfb`, plus the RPA driver
  that replaces `nj_hf_rrpa_driver`/`nj_hf_rpa_driver`. The HF-EFA caller is removed in dev.
* Optionally replace `allocate(Buf(N_Chunk(a)))` with a 1e6-element block buffer, as dev does. The results are identical.
  Today the largest chunk is allocated whole. The fix is ~20 lines.
* Keep the port extensions: `CM_Lawson` (Lawson term, not in Julia), and `T%p(a,c)` instead of `T[a,c]`
  (deviation 2).

### 1.2 `V2b_parse` (NEW, V2b.jl:307-458)
`V2b_parse(Params, N_max_Out, N_2max_Out, Out_File)` cuts `Params.Int.NN_File` (space Int.Nmax/N2max) down to a
smaller space and writes it in the same NuHamil Float64 layout.
* It errors if `N_max_Out > Int.Nmax || N_2max_Out > Int.N2max`, or if `abspath(In)==abspath(Out)`, or if
  `filesize(In) != 8*sum(N_Chunk_In)`. The size check is new; `V2b_read` has none.
* `Orb = orbitals_make(Params, Int=true)`. It counts both spaces with `V2b_count`. `a_max_Out = (N_max_Out+1)(N_max_Out+2)/2`.
* The output file is pre-sized with `truncate(Out, 8*sum(N_Chunk_Out))`. Each chunk `a in 1:a_max_Out` is processed
  independently. The queue is ordered largest-first, which only affects threading.
* Per chunk it loops with the **In** limits (b<=a with N_ab<=N_2max_In, c<=a, d<=d_max, N_cd<=N_2max_In, parity) and gets
  `J_min, J_max` as in the reader. If `N_a,N_b,N_c,N_d <= N_max_Out && N_ab,N_cd <= N_2max_Out`, it copies
  `4*max(J_max-J_min+1,0)` Float64. Otherwise it skips them. The output is written at `N_skip_Out[a]`, in order.
* Check at the end of each chunk: `Chunk_i == N_In && Out_N == N_Out`, else error "Wrong # of matrix elements in the chunk a = ...".

Fortran: new `subroutine V2b_parse(Params,N_max_Out,N_2max_Out,Out_File)` (~110 lines). Use stream I/O with
`pos=` and OpenMP `schedule(dynamic)` over chunks, ordered by `sortperm_rev(N_Chunk_Out)`. Each thread needs its own
units. Writing disjoint byte ranges of one file from several units opened `access='stream', action='readwrite'` works
with ifx. The simplest fallback is to process the chunks serially in order with one output unit; the bytes are identical.
The file must exist and be pre-sized; Fortran cannot truncate, so either write sequentially or write a final zero at
`pos=Size_Out-7` first. An entry point in the app is needed (see open question Q2).

### 1.3 `V2b_residual_no2b` is NON-MUTATING now
Dev V2b.jl:474: `V_NN_res = deepcopy(V_NN)`. All `+=`/`*= cRes` go to `V_NN_res`, which is returned. The loop structure,
formulas and the `cRes` placement (`V_res = (V + Sum)*cRes`) are the same as base.
Callers: HF_solver.jl:45 `V_NN_den = V2b_residual_no2b(...)`, which **keeps the bare `V_NN`** for
`HF_MBPT_export(..., V_NN, V_NN_res, V_NNN)`. BCS_solver.jl:45 and HFB_solver.jl:46 assign back to `V_NN` (same effect as before).

Fortran: change to
`subroutine V2b_residual_no2b(Params,Orb,Orb_NN,Orb_NNN,Rho,V_NN,V_NNN,V_NN_res)` with `V_NN` `intent(in)`,
`V_NN_res` `intent(out)` (`V_NN_res = V_NN`, then the same loops writing into `V_NN_res`). In nj_bcs/nj_hfb use a
temporary and `call move(V_NN_res -> V_NN)`. In the HF driver keep both. The port note "In-place mutation ...
(see HF_EFA)" becomes obsolete for this routine. ~15 changed lines plus callers.

### 1.4 `V2b_correlation_energy` (NEW in common, V2b.jl:572-635; it was only in base `Temp/`)
```
E_corr = sum_{a,b: N_ab<=N2max} sum_{c,d: N_cd<=N2max, P_ab==P_cd} sum_{J=max(|ja-jb|,|jc-jd|)/2}^{min(ja+jb,jc+jd)/2}
         (2J+1) * [ 1/4 O2b_pp(abcdJ;V) O2b_pp(abcdJ;Sigma) + 1/4 O2b_nn(V) O2b_nn(Sigma) + O2b_pn(V) O2b_pn(Sigma) ]
```
a,b,c,d each run over `1:a_max` (Calc), giving all ordered pairs. `P = rem(l_a+l_b,2)+1`. Threads use `:static` with
`threadid()` accumulators, which is race-free. The only dev caller is RPA_energy.jl:227 (`RPA_energy_rho`). It prints
"Computing the correlation ground-state energy ..." and "Performed calculation ...".
Fortran: new `real(dp) function V2b_correlation_energy(Params,Orb,Orb_NN,Sigma,V_NN)` with an OpenMP reduction over a (~50 lines).

---------------------------------------------------------------------------------------------------------------

## 2. V3b_no2b.jl (667 lines) <- base V3B_NO2B.jl (432) ; Fortran `nj_v3b.f90`

### 2.1 Orb3B structure: Dict replaced by a dense map (import_structures.jl:243-250)
```
struct Orb3B
    Map::Vector{Array{Int32,3}}   # Map[ab][TT, J-J_min+1, c],  ab = div(a*(a-1),2)+b (b<=a), J_min = |j_a-j_b|/2
    N::Array{Vector{Int64},4}     # N[P,l+1,J+1,T_ind][j_index]  (unchanged)
    nlj::Array{Int64,3}           # nlj[n+1,l+1,j_index] = orbital index a
end
```
* `TT = T_ab + T_ind`, where `T_ind = (T+1)/2` with T doubled. That gives TT = 1 for (T_ab,T)=(0,1/2), 2 for (1,1/2), 3 for (1,3/2).
  The value 0 means not allowed.
* `Map[ab]` is allocated for **every** `1<=b<=a<=a_max` as `zeros(Int32,3,J_count,a_max)` with
  `J_count = (N_ab<=N_2max) ? (j_a+j_b)/2-|j_a-j_b|/2+1 : 0` (V3b_no2b.jl:11-18).
* `nlj` has size `(div(N_max,2)+1, l_max+1, 2)`. `j_index = 1` if `j == |2l-1|`, else 2 (V3b_no2b.jl:21-25).
* `V3b_no2b_initialize` uses the same loops and counting as base. Instead of `Dic[key] = N[...]` it stores
  `Map_ab[T_ab+T, J-J_min+1, c] = N_Orb_NNN[P,l_c+1,J+1,T][j_index]` (V3b_no2b.jl:61-71). The block sizes `N_Orb_NNN`
  and the V_NNN allocation (`Float32`, length `N(N+1)/2`) are unchanged. `V3b_no2b_key` is no longer used in dev.

### 2.2 Lookup and packed index (CHANGED STORAGE ORDER)
`V3b_no2b_lookup(a,b,c,T_ab,J,T_ind,Orb,Orb_NNN)` (V3b_no2b.jl:298-310) returns 0 if `b>a || b<1 || ab>length(Map) ||
T_ab<T_ind-1 || T_ab>1 || J<J_min || J-J_min+1>size(Map_ab,2) || c<1 || c>size(Map_ab,3)`. Otherwise it returns
`Map_ab[T_ab+T_ind, J-J_min+1, c]`.

`V3b_no2b_index` (V3b_no2b.jl:312-340):
```
Bra = lookup(a,b,c,T_ab,J,T_ind);  Bra==0 -> KeyError
f_c = (l_f==l_c && j_f==j_c) ? f : (n_f+1 <= size(nlj,1) ? nlj[n_f+1,l_c+1,j_index(c)] : 0)   # ket orbital (n_f,l_c,j_c)
Ket = f_c>0 ? lookup(d,e,f_c,T_de,J,T_ind) : 0
P_abc = parity(l_a+l_b+l_c); P_def = parity(l_d+l_e+l_c)
Ket==0 || P!=P_abc || P!=P_def -> KeyError
if Bra < Ket: swap
Ind = Ket + div(Bra*(Bra-1),2)        # <-- packed lower triangle stored ROW BY ROW (dev :338)
```
Base used `Ind = Bra + (Ket-1)*N - div(Ket*(Ket-1),2)` (column by column). The set of stored elements and the values
returned by `V3b_no2b(...)` are identical. Only the position inside each `V_NNN[P,l,J,T][j]` vector changes.
**Hamiltonian.jl indexes V_NNN blocks directly with the new formula** (`H1b_V3b_index`, `H2b_V3b_index`), so the
Fortran storage must switch to the row-wise order. The KeyError semantics are the same as base: Dict misses there
covered both wrong P and missing keys.

`V3b_no2b(...)` (phases for a<b / d<e) is unchanged (dev :342-379).

### 2.3 `V3b_no2b_read(Params)`
* New signature `V3b_no2b_read(Params) -> (V_NNN, Orb_NNN)`. It builds `Orb = orbitals_make(Params, Int=true)` itself (:396).
* `Mmap` was dropped. It uses one Float32 buffer per thread (1e6) with `seek(Int_IO, N_skip + 4*Buffer_0)` + `read!`.
* Only (a,b) chunks with `N_a<=Calc.Nmax && N_b<=Calc.Nmax && N_a+N_b<=Calc.N2max` are queued (:409-418), largest first.
* Inside a chunk (`V3b_no2b_read_chunk!`, :441-552), if any of the Calc-space conditions fails for (c,d,e,f), the
  `5*max(J_max-J_min+1,0)` elements are **skipped without reading** (:537-540). Otherwise it loops J, T_ab, T_de, T
  (5 combinations per J) and stores `V_NNN[P,l_c+1,J+1,T_ind][j_index][Ind] = V`.
* `V3b_no2b_count`: count per (c,d,e,f) is `N += 5*max(J_max-J_min+1,0)` (:645-648), with a running-sum byte skip (:661-665).
  The numbers are identical to base.
* The result is identical to base and to the current Fortran (which already skips non-Calc chunks via `ab_last` and cycles).

Fortran (`nj_v3b.f90`):
1. `orb3b_t` (nj_types): add `type(i3_t), allocatable :: Map(:)` (each `integer(int32), allocatable :: m(:,:,:)`,
   dims `(3, J_count, a_max)`) and `integer, allocatable :: nlj(:,:,:)` with dims `(0:N_max/2, 0:N_max, 2)`. Keep `N`.
   `Dic` is no longer needed for V3b (see Q1).
2. `V3b_no2b_initialize`: allocate `Map` and `nlj` as in 2.1, and fill `Map` instead of `Dic%set` (~40 lines).
3. New `pure integer function V3b_no2b_lookup(a,b,c,T_ab,J,T_ind,Orb,Orb_NNN)` with the bound checks above (~20).
4. `V3b_no2b_index`: dense lookup + `nlj` ket orbital + parity check + **`Ind = Ket + Bra*(Bra-1)/2`** (~25). Keep the
   out-of-line cold error routine `V3b_no2b_missing_key`, which is part of the uncommitted fix. Compute
   `Bra*(Bra-1)/2` and the allocation length `N*(N+1)/2` in `integer(i8)`, or at least guard N <= 46340. The current
   default-integer `N_count*(N_count+1)/2` overflows beyond that.
5. `V3b_no2b_read(Params,V_NNN,Orb_NNN)`: build the Int orbitals internally. Add the "skip non-Calc (c,d,e,f) blocks
   without consuming the buffer" logic: advance `Pos` by `4*5*nJ` and keep the buffer logic consistent. The simplest way is
   to track `Chunk_i` and refill from `pos = 4*(skip+Buffer_0)+1` the way dev does. This only matters for I/O volume.
   Keep the port extensions: empty `NNN_File` gives zero NNN, plus the file-size check.
6. New `V3b_no2b_parse(Params,N_max_Out,N_2max_Out,N_3max_Out,Out_File)` (dev :114-294) works like `V2b_parse`:
   * containment checks for (Nmax,N2max,N3max);
   * `Size_In = N_chunk_skip_In[end] + 4*N_chunk_In[end]` must equal the file size;
   * `ab_In_index[ab_o]` is a monotone search for the same (a,b) in the In chunk list;
   * the inner loops use the In limits (c over the In `a_max`) and copy `5*nJ` Float32 when all 14 Out conditions hold
     (`N_a..N_f<=N_max_Out`, `N_ab,N_ac,N_bc,N_de,N_df,N_ef<=N_2max_Out`, `N_abc,N_def<=N_3max_Out`), else skip;
   * per-chunk check `Chunk_i==N_In && Out_N==N_Out` (~130 lines).
7. `V3b_no2b_key` can be deleted, or kept for the Dic (Q1).

### 2.4 Shared-memory and race notes
* Dev V3b read: chunks (a,b) and (a,e) with the same leading `a` both write the elements {(a b c),(a e f)} (the file stores
  both triangles for d==a). This is benign only if the NuHamil file is exactly symmetric. Otherwise both Julia and Fortran
  keep the value of whichever thread writes last. This is pre-existing and unchanged, so no deviation is needed. A note in
  PORTING_NOTES is optional.

---------------------------------------------------------------------------------------------------------------

## 3. Tr1b.jl (1028 lines) <- base Tr1B.jl (860) ; Fortran `nj_tr1b.f90`

### 3.1 Tr1B / qpTr1B structures gain F and GT (import_structures.jl:309-347)
`Tr1B(E0,E1,E2,E3,E1_VC,E1_VS,E1_TC,E1_TS,E1_sTC,E1_sTS,E1_C,E1_sC,M1,M2,M3,F,GT)`, with the same fields in `qpTr1B`.
Fortran: add `F, GT` to `tr1b_t` and `qptr1b_t` (declared in nj_tr1b, not nj_types). The dev code does not consume
F/GT anywhere outside Tr1b.jl at f3f7a33. pnQTDA (another group) may need them.

### 3.2 Toroidal/vortical E1 operators: BUG FIX in dev ("toroidal fix")
`rME_nabla_dot_rN_YJL` (dev :167-205) is now the convective combination `<a||A.(psi_a^+ grad psi_b - (grad psi_a)^+ psi_b)||b>`:
| term | base | dev |
|---|---|---|
| Bra, l_a-1 (`cg(2l_b,2L,2l_a-2)`, `sqrt(l_a/(2l_a+1))`) | `rME += ME` | `rME -= ME` |
| Bra, l_a+1 | condition `... && l_a >= 1`, factor `sqrt((l_a-1)/(2l_a+1))`, `rME -= ME` | no `l_a>=1`, factor **`sqrt((l_a+1)/(2l_a+1))`**, `rME += ME` |
| Ket, l_b-1 | `+=` | `+=` (unchanged) |
| Ket, l_b+1 | condition `&& l_b >= 1`, `sqrt((l_b-1)/(2l_b+1))`, `-=` | no `l_b>=1`, **`sqrt((l_b+1)/(2l_b+1))`**, `-=` |

`rME_nabla_cross_S_dot_rN_YJL` (dev :215-217) uses the normalization `/(16.0*pi*J*(J+1))` instead of `/(8.0*pi*J*(J+1))`,
so it is normalized to s = sigma/2.

Operator prefactors (dev :293-327). Only the convective parts get `0.5*`:
```
VC  = 0.5 * hc_mc2 * sqrt(3/2)/5 * rG_dot_r2_Y12
TC  = 0.5 * hc_mc2 * (sqrt(3/2)/15 * rG_dot_r2_Y12 + sqrt(3)/6 * rG_dot_r2_Y10)
sTC = -0.5 * hc_mc2 * sqrt(3)/6 * chR2 * rG_dot_Y10
VS, TS, sTS, C, sC : formulas unchanged (VS/TS/sTS change only through the 1/16pi above)
```
Units comment (dev :288-291): combine as `VC + 0.5*g_s*VS` downstream (RPA_rB). Downstream strength formulas in the
RPA/QRPA groups changed accordingly and must be ported **together** with this change.
Fortran: `rME_nabla_dot_rN_YJL` (4 branches), `rME_nabla_cross_S_dot_rN_YJL` (8 -> 16), and `Tr1b_E1_NLO_initialize`
(VC/TC/sTC `0.5_dp*` in front: `0.5_dp * hc_mc2 * ...`, which multiplies exactly, so the results are bitwise 0.5x).
About 15 lines. The CC/sC/VS/TS/sTS lines are otherwise unchanged.

### 3.3 New `Tr1b_BX_initialize` (dev :470-535): Fermi and Gamow-Teller
```
rF(a,b)  = (n_a==n_b && l_a==l_b && j_a==j_b) ? sqrt(j_a+1) : 0
rGT(a,b) = (n_a==n_b && l_a==l_b) ? sqrt(2*(j_a+1)*(j_b+1)) * phase(l_a + (j_a+3)/2) * w6j(1,1,2,j_b,j_a,2*l_a) : 0
F.p = F.n = rF   (only when l,j,n equal);   GT.p = GT.n = rGT  (only when l_a==l_b, n_a==n_b, |j_a-j_b|<=2)
```
`Tr1b_initialize` calls it after `Tr1b_MX_initialize`. Its printed header changed to "Constructing electromagnetic &
weak transition operators E^lambda, M^lambda, M^F & M^GT in the LHO basis ...", plus the corresponding "... succesfully
constructed" line.

### 3.4 `Tr1b_transformation` (dev :536-771): F and GT mix the isospin bases
```
F_new.p[a,b] = sum_{k~a, l~b} F.p[k,l] * C.p[k,a] * C.n[l,b]
F_new.n[a,b] = sum            F.n[k,l] * C.n[k,a] * C.p[l,b]      (same for GT)
```
All other operators transform as before. Fortran: add `transform_pn(O)` (same lj-block loops, with mixed C) for F and GT (~25 lines).

### 3.5 Quasiparticle transforms (qpTr1b_*), CHANGED FORMULAS
`Phase_11 = (-1)^((j_a+j_b)/2 + X)`, `Phase_20 = (-1)^X` (unchanged). `O_kl = TrOp.p[k,l]` and now also `O_lk = TrOp.p[l,k]`.
* **EX, X>0** (dev :875-880):
  `11 = O_kl U_ka U_lb + Phase_11 * O_lk V_ka V_lb` (base: `O_kl (U U - Phase_11 V V)`).
  `20 = O_kl (Phase_T V_ka U_lb + Phase_20 U_ka V_lb)`, with `Phase_T = T_odd ? -1 : +1`.
* **E0 (X=0)** (dev :859-866): `11 = O_kl U U + Phase_11 O_lk V V`. The c-number term `+ j_a_hat*delta_ab*V V` is
  **removed**. `20 = O_kl (V U + U V)` is unchanged, with no Phase_T.
* **MX** (dev :943-948): `11 = O_kl U U + Phase_11 O_lk V V` (changed). `20 = O_kl (V U - Phase_20 U V)` (unchanged).
* `T_odd = true` only for `E1_VC, E1_VS, E1_TC, E1_TS, E1_sTC, E1_sTS` (dev :789-796). `E1_C, E1_sC, E0..E3` stay false.
* New `qpTr1b_BX_initialize(Params,Orb,U,V,TrOp_BX,X)` (dev :970-1028), with F at X=0 and GT at X=1. It requires
  `|j_a-j_b|/2 <= X <= (j_a+j_b)/2` and `rem(l_a+l_b,2)==0`:
```
p11 += BX.p[k,l] pU_ka nU_lb + Phase_11 * BX.n[l,k] pV_ka nV_lb
n11 += BX.n[k,l] nU_ka pU_lb + Phase_11 * BX.p[l,k] nV_ka pV_lb
p20 += BX.p[k,l] pU_ka nV_lb
n20 += BX.n[k,l] nU_ka pV_lb
```
The qpTr1b header strings now say "electromagnetic & weak".
Fortran `qp_transform`: add an `O_lk` read, the new 11 formulas, a `T_odd` argument for the 20 term, removal of the E0
c-number, and a new `qp_transform_bx` (~60 lines total). Callers are nj_qtda and nj_qrpa (unchanged signatures).

---------------------------------------------------------------------------------------------------------------

## 4. Hamiltonian.jl (NEW, 1449 lines) ; Fortran counterpart: none -> new module `src/core/nj_hamiltonian.f90`

Functions, all with `Orb` the Calc orbitals, `Rho` the LHO-basis OBDM (only lj-diagonal blocks used), and `V_NN`/`V_NNN` bare:

| function | dev lines | used by (dev) | port? |
|---|---|---|---|
| `H1b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho)` -> O1B | 443-597 | RPA_Hamiltonian.jl:17, HF_MBPT.jl:2274,2284 | yes |
| `H2b(Params,Orb,Orb_NN,Orb_NNN,V_NN,V_NNN,Rho;Report=true)` -> O2B | 599-760 | RPA_Hamiltonian.jl:18, HF_MBPT.jl:2275, H1b_H2b | yes |
| `H1b_H2b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho,Sigma;Report=true)` -> (O1B,O2B) | 995-1449 | RPA_Hamiltonian.jl:10, HF_MBPT.jl:2269 | yes |
| `H1b(...,Rho,Sigma)` (9 args) | 2-441 | **none** | no (unreachable) |
| `E1b(...)`, `E2b(...)` | 762-993 | **none** | no (unreachable) |

### 4.1 `H1b` (8 args): the mean field h_N[Rho]
This is the same physics and formula as base `HF_allocate` (`Source/HF/HF_Allocate.jl:17-107`, Fortran `nj_hf::HF_allocate`).
* `lj_pairs`: all ordered (a,b) with equal (l,j), a outer and b inner, giving `ad = lj_pairs_low` (b<=a), `be = lj_pairs`, `cf = lj_pairs`.
* The work list is `adbe_tasks = [(ad_i,be_i) for ad_i, be_i]`. There are `N_task = 4*nthreads` tasks, each taking every
  N_task-th element, with task-local accumulators `ph_N_local[ad_i,task]`.
* Per (a,d,b,e): requires `N_ab<=N2max && N_de<=N2max`, then loops `J = |j_d-j_e|/2 : (j_d+j_e)/2`.
  * NN part, if `P_ab==parity(l_d+l_e)`: `J_j_a_hat = (2J+1)/(j_a+1)`, accumulating
    `p += Jh*pRho_be*O2b_pp(a,b,d,e) + Jh*nRho_be*O2b_pn(a,b,d,e)` and `n += Jh*nRho_be*O2b_nn(a,b,d,e) + Jh*pRho_be*O2b_pn(b,a,e,d)`.
  * NNN part: for (c,f) in cf with `N_abc<=N3max && N_def<=N3max && P_abc==parity(l_d+l_e+l_f)`, using the five
    `V3b_no2b` elements ME001, ME101, ME011, ME111, ME113, exactly as in HF_allocate (Hamiltonian.jl:548-554).
* Then add T: `CMS1+2B` gives `T*(1-1/A)`, `CMS != CMS2B` gives `T`, and CMS2B adds nothing. Symmetrize `ph_N[a,d]=ph_N[d,a]`.
* Prints "Evaluating 1-body mean-field Hamiltonian h_N ..." and "\tCompleted the calculation of 1-body mean-field Hamiltonian h_N ...".
Fortran: port `H1b` into nj_hamiltonian as an OpenMP loop over (ad,be) pairs with a reduction into `ph(ad_count)`
(~110 lines). Do NOT call `nj_hf::HF_allocate`: core must not depend on hf, and the dev HF_allocate got extra
l_max/j_max pruning (HF group).

### 4.2 `H2b`: density-dependent residual NN (NO2B), the fast equivalent of `V2b_residual_no2b`
* `H_NN = deepcopy(V_NN)`. Loops `J in 0:J_max (=N2max+1)`, `P in 1:2`, and threads over `Bra in 1:N_pp = Orb_NN.N[2,P,J+1]`
  (pp pairs a>=b). Each Bra writes only its own elements, so it is race-free (checked: pn element {xy,zw} is reached from exactly one
  (Bra>=Ket) pp pair and one of the 4 orderings).
* Per Bra=(a,b) and Ket=(d,e) in `1:Bra`:
  `Phase_ab = phase((j_a+j_b)/2 - J)`, `Phase_de` likewise, `Map_ab = Orb_NNN.Map[ab]`, `J_ab = J - |j_a-j_b|/2 + 1`, `J_de` likewise.
  For `(c,f) in cf` (all ordered lj-equal pairs, using `Orb_c.a` order):
  ```
  abc_0,def_0 = Map_ab[1,J_ab,c], Map_de[1,J_de,f];  skip if either == 0
  abc_1,def_1 = Map[2,...]; abc_3,def_3 = Map[3,...]
  V_T1 = V_NNN[P3B,l_c+1,J+1,1][j_index[c]], V_T3 = V_NNN[P3B,l_c+1,J+1,2][j_index[c]],  P3B = parity(l_a+l_b+l_c)
  ME001 = V_T1[idx(abc_0,def_0)], ME101 = V_T1[idx(abc_1,def_0)], ME011 = V_T1[idx(abc_0,def_1)],
  ME111 = V_T1[idx(abc_1,def_1)], ME113 = V_T3[idx(abc_3,def_3)]     # idx(B,K) = B>=K ? K + B(B-1)/2 : B + K(K-1)/2
  ppSum += ME113*pRho_fc + (2ME111+ME113)/3*nRho_fc
  nnSum += ME113*nRho_fc + (2ME111+ME113)/3*pRho_fc
  pnSum_abde += (0.5ME001 + ME111/6 + ME113/3)(p+n) + ( ME101+ME011)/sqrt(12)(n-p)
  pnSum_bade += (0.5ME001 - ME111/6 - ME113/3)(p+n) + (-ME101+ME011)/sqrt(12)(n-p)
  pnSum_abed += (0.5ME001 - ME111/6 - ME113/3)(p+n) + ( ME101-ME011)/sqrt(12)(n-p)
  pnSum_baed += (0.5ME001 + ME111/6 + ME113/3)(p+n) - ( ME101+ME011)/sqrt(12)(n-p)      (p,n = pRho_fc, nRho_fc)
  ```
  The writes use `is_J_hat = 1/(2J+1)` (Hamiltonian.jl:728-749):
  ```
  pp/nn[O2b_index(a,b,d,e,J,P,1)] = cRes*(H + is_J_hat*pp/nnSum)
  pn[(a,b,d,e)] = cRes*(H + is_J_hat*pnSum_abde)
  if a!=b:                 pn[(b,a,d,e)] = cRes*(H + is_J_hat*Phase_ab*pnSum_bade)
  if d!=e && Ket!=Bra:     pn[(a,b,e,d)] = cRes*(H + is_J_hat*Phase_de*pnSum_abed)
  if a!=b && d!=e:         pn[(b,a,e,d)] = cRes*(H + is_J_hat*Phase_ab*Phase_de*pnSum_baed)
  ```
* If `Report`: prints "\nCalculating density dependent residual 2-body interaction Hamiltonian H_NN ... (NO2B approximation
  of V_NNN)", then a per-(J,P) line `"\tCalculating   ...   J = J/J_max\tP = +|-"`, then "\tCompleted ...".
* The result equals `V2b_residual_no2b` up to summation order (I checked the pn formula algebraically against V2B.jl:516-519).
Fortran: `subroutine H2b(Params,Orb,Orb_NN,Orb_NNN,V_NN,V_NNN,Rho,H_NN,Report)` (~120 lines). It needs `Orb_NNN%Map`
and the row-wise V_NNN packing (section 2). It uses Fortran `Orb_NN%Dic(a,b,J,T)`, which is the same as dev
`Orb_NN.Map[a,b,J+1,T+1]`.

### 4.3 `H1b_H2b`: the Sigma-dressed mean field (used when `RPA.ScH && RPA.TBCF`, or `HF.MBPT.ScH && TBCF`)
1. `H_NN = H2b(..., Rho, Report)`.
2. `h_N` = the same Rho part as `H1b` (4.1), with **bare `V_NN`** and T/CMS.
3. Plus the NNN-Sigma contraction (Hamiltonian.jl:1226-1336):
   * `NNN_triples[P,l_c+1,J+1][j_index]`: a list of `(a,b,c,ab,i0,i1,i3)` for every `a>=b` (a outer, b inner), every
     `J_i in 1:size(Map_ab,2)` (`J = J_min+J_i-1`), and every `c in 1:a_max` with `Map_ab[1,J_i,c]!=0`.
     `ab = Orb_NN.Map[a,b,J+1,2]` (pp index), `i0,i1,i3 = Map_ab[1..3,J_i,c]`, `P = parity(l_a+l_b+l_c)`.
   * For `J in 0:J_max` and `P2B in 1:2`: skip if `N_ab = Orb_NN.N[2,P2B,J+1] == 0`. Unpack Sigma (pp, nn: N_ab^2; pn:
     N_ab_pn^2) into dense symmetric matrices (`H1b_unpack`, column-packed `Bra + (Ket-1)N - Ket(Ket-1)/2`).
     Threads are tasks over the triples. For `l_c in 0:N_max`, `P3B = parity(P2B-1+l_c)`, and each `j_index`: with
     `V_T1 = V_NNN[P3B,l_c+1,J+1,1][j]`, `V_T3 = ...[2][j]`, loop `abc` in task:N_task:end and `def in 1:abc`:
   ```
   ab_pn, ba_pn = Orb_NN.Map[a,b,J+1,1], Orb_NN.Map[b,a,J+1,1];  de_pn, ed_pn likewise
   Phase_ab = phase((j_a+j_b)/2 - J); Phase_de likewise
   ME001=V_T1[idx(i0,i0')], ME101=V_T1[idx(i1,i0')], ME011=V_T1[idx(i0,i1')], ME111=V_T1[idx(i1,i1')], ME113=V_T3[idx(i3,i3')]
   ppS = ppSigma[de,ab] * (a!=b ? 2 : 1) * (d!=e ? 2 : 1);  nnS likewise
   s_abde = pnSigma[de_pn,ab_pn]
   s_bade = a!=b ? Phase_ab*pnSigma[de_pn,ba_pn] : 0
   s_abed = d!=e ? Phase_de*pnSigma[ed_pn,ab_pn] : 0
   s_baed = (a!=b && d!=e) ? Phase_ab*Phase_de*pnSigma[ed_pn,ba_pn] : 0
   S00 = s_abde+s_bade+s_abed+s_baed;  S10 = s_abde-s_bade+s_abed-s_baed
   S01 = s_abde+s_bade-s_abed-s_baed;  S11 = s_abde-s_bade-s_abed+s_baed
   pSum = ME113*ppS + (ME113+2ME111)/3*nnS + 2*(ME001*S00 - (ME101*S10+ME011*S01)/sqrt(3) + (ME111+2ME113)/3*S11)
   nSum = ME113*nnS + (ME113+2ME111)/3*ppS + 2*(ME001*S00 + (ME101*S10+ME011*S01)/sqrt(3) + (ME111+2ME113)/3*S11)
   ph[c,f] += pSum; nh[c,f] += nSum;  if def != abc: ph[f,c] += pSum; nh[f,c] += nSum
   ```
   * Final: for (c,f) in cf: `h_N[c,f] += 0.25/(j_c+1) * ph_cf[c,f]`.
4. The V_NN-Sigma part (`H1b_dense`) is **commented out** in `H1b_H2b` (Hamiltonian.jl:1338-1441). Do not port it.
5. Returns `(O1B(ph_N,nh_N), H_NN)`. If `Report`, it prints "\nEvaluating 2-body Hamiltonian H_NN ...", "\nEvaluating 1-body
   mean-field Hamiltonian h_N ...", and "\tCompleted the calculation of 1-body mean-field Hamiltonian h_N ...".
Fortran: `subroutine H1b_H2b(Params,Orb,Orb_NN,Orb_NNN,T,V_NN,V_NNN,Rho,Sigma,h_N,H_NN,Report)` (~220 lines). Store the
triples per (P,l,J,j) as `integer, allocatable :: t(:,:)` of shape (7,n), built in two passes (count, fill). Use an OpenMP
`parallel do schedule(dynamic)` over `abc` with per-thread `ph(a_max,a_max)` reductions. The dense Sigma blocks are
N_pn^2 doubles per (J,P), roughly 70 MB at N_pn ~ 3000, built one block at a time.

### 4.4 Not ported (unreachable in dev)
`H1b(...,Rho,Sigma)` (9 args, including the V_NN-Sigma `H1b_dense` term), `E1b`, `E2b`. Add them to the "Not ported" list.

---------------------------------------------------------------------------------------------------------------

## 5. Shared-module interface changes (src/core, nj_types, nj_params, nj_input)

S1. `nj_orbitals::orbitals_make(Params, Int)`: add `logical, intent(in), optional :: Int`. The default uses **Calc.Nmax**
    (dev orb.jl:1-5); `Int=.true.` uses Int.Nmax. Owned by group common A. This group needs it in V2b_read,
    V3b_no2b_read and the parse routines.
S2. `nj_o1b::T1b`: dev still sizes by `Params.Int.Nmax` (T1b.jl:3) while dev solvers now pass Calc-size orbitals. That is a
    Julia out-of-bounds read when Int.Nmax > Calc.Nmax. Fortran: size by `size(Orb)`. The used entries are identical; record it as a new deviation.
S3. `nj_types::orb3b_t`: add `Map(:)` (derived type holding `integer(int32), allocatable :: m(:,:,:)`) and `nlj(:,:,:)`. `N` is unchanged.
    `Dic` (hashmap_t) is removed or kept depending on Q1.
S4. `nj_types::v3b_t`: layout unchanged (`V(P,0:l,0:J,T,j)%v`), but the packed order inside `%v` becomes row-wise
    `Ind = Ket + Bra*(Bra-1)/2`.
S5. `nj_v2b`: `V2b_read(Params,V_NN,Orb_NN)` (Orb argument dropped); `V2b_residual_no2b(...,V_NN [in],V_NNN,V_NN_res [out])`;
    new `V2b_parse`, `V2b_correlation_energy`; `V2b_count` unchanged.
S6. `nj_v3b`: `V3b_no2b_read(Params,V_NNN,Orb_NNN)` (Orb dropped); new `V3b_no2b_lookup`, `V3b_no2b_parse`; public `V3b_no2b_key` deleted.
S7. `nj_tr1b`: `tr1b_t`/`qptr1b_t` get `F, GT`. qp transforms get the internal `T_odd` flag. Public signatures are unchanged.
S8. New `nj_hamiltonian` (core): `H1b`, `H2b`, `H1b_H2b`. It depends on nj_o2b (O2b_pp/pn/nn/index), nj_v3b and nj_types. Add it to the Makefile.
S9. `nj_input`/app: a parse action (`V2b_parse`/`V3b_no2b_parse`) needs `N_max_Out, N_2max_Out, N_3max_Out, NN_Out, NNN_Out`
    (Q2). The `&interaction` namelist is unchanged.

## 6. Port notes: status changes
* Translation table row "`Dict` of 3-body NO2B orbitals -> open-addressing hash map": dev now uses a dense map (`Orb3B.Map`), so the port follows.
  Row "`Dict{UInt64,Int64}` of 2-body orbitals -> dense table": dev now also has a dense `Orb2B.Map[a,b,J+1,T+1]`, so it matches.
* "In-place mutation ... `V2b_residual_no2b`": obsolete. Dev returns a copy, and HF_EFA is removed.
* Deviation 2 (CMS2B `T[a,c]`): **still needed**, since dev is unchanged.
* Deviation 8 (Verbose silencing inside RRPA): RRPA was removed in dev. `H2b` has a `Report` keyword that replaces it.
* New deviation: T1b is sized by `size(Orb)` (S2).
* The Tr1b E1 current fix (3.2) was **fixed upstream**, so the base behaviour must be replaced, not kept.

## 7. Size estimate (Fortran lines, new or changed)
nj_v2b ~230 (parse 110, corr. energy 50, read/residual 70) ; nj_v3b ~300 (map/lookup/index 90, read 40, parse 130, misc) ;
nj_tr1b ~170 ; nj_hamiltonian ~450 ; types/orbitals/T1b ~40 ; callers ~40. **Total ~1250.**
