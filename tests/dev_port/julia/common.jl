# Shared helpers of the dev-port Julia reference drivers (tests/dev_port/julia/<case>.jl).
#
# The drivers run the PATCHED Julia reference (branch port-ref = dev f3f7a33 + fixes F1..F5) found in
# $NUJULIET_REF (default /scratch/local.ssd/knapp/job_25503048.pbs-m1/NuJuliet_dev_ref). They must be started from
# the case work directory (run_julia.sh does this): dev writes everything to io/<Path>/... relative to the cwd.
#
# Dev solvers report most input errors with println + return (no exception), so every step is followed by a check
# that its output exists; a failed check aborts the driver with an error.

const NUJULIET_REF = get(ENV, "NUJULIET_REF", "/scratch/local.ssd/knapp/job_25503048.pbs-m1/NuJuliet_dev_ref")
include(joinpath(NUJULIET_REF, "src", "import_NuJuliet.jl"))
using .NuJuliet
using LinearAlgebra, Printf, Dates

# Live log: when stdout goes to a file Julia buffers it, a task on the interactive thread (julia -t N,1) flushes it.
const FLUSH_INTERVAL = 10.0
const OUT = stdout
const ERR = stderr
if Threads.nthreads(:interactive) > 0
    Threads.@spawn :interactive while true
        sleep(FLUSH_INTERVAL)
        flush(OUT)
        flush(ERR)
    end
else
    println("WARNING: no interactive thread (start Julia with -t N,1) ... the log is flushed only after each step")
end

# BLAS threads = Julia threads (OpenBLAS would otherwise use all cores) ...
BLAS.set_num_threads(Threads.nthreads())

isdir("io") || mkdir("io")

# Bundled small interaction of the dev repository (hw = 16, Nmax = 3, N2max = 6, N3max = 9) ...
const BUNDLED_INT = Interaction_Parameters(
    NN_File  = joinpath(NUJULIET_REF, "io", "NN.bin"),
    NNN_File = joinpath(NUJULIET_REF, "io", "NNN.bin"),
    hw = 16.0, Nmax = 3, N2max = 6, N3max = 9)

# DN2LOGO394 NN + NNN files (file truncation emax = 16, e2max = 32, e3max = 28) ...
const DN2LOGO394_DIR = "/auto/brno12-cerit/nfs4/home/knapp/NuHamil_files/NN+NNN files/DN2LOGO394/"
DN2LOGO394_INT(hw::Int) = Interaction_Parameters(
    NN_File  = DN2LOGO394_DIR * "NN-DN2LOGO394_vbare2.0_hw$(hw)_emax16_e2max32.me2j.bin",
    NNN_File = DN2LOGO394_DIR * "NNN_DN2LOGO394_hw$(hw)_hwt$(hw)_emax16_e2max32_e3max28_jmax3_59.stream.bin",
    hw = Float64(hw), Nmax = 16, N2max = 32, N3max = 28)

const TIMINGS = Tuple{String,Float64}[]

println("Julia ", VERSION, " , threads = ", Threads.nthreads(), " (+ ", Threads.nthreads(:interactive), " interactive)",
        " , BLAS threads = ", BLAS.get_num_threads(), " , reference = ", NUJULIET_REF)
println("cwd = ", pwd(), " , start ", now())
flush(OUT)

# Run one step, time it, check its outputs. `expect` lists paths (relative to the cwd) that must exist afterwards.
function step(f::Function, label::String; expect::Vector{String} = String[])
    println("\n", "="^100, "\nSTEP ", label, "  (", now(), ")\n", "="^100)
    flush(OUT)
    t = @elapsed f()
    push!(TIMINGS, (label, t))
    println(@sprintf("\nSTEP_TIME %-40s %10.2f s", label, t))
    flush(OUT)
    missing_paths = filter(p -> !ispath(p), expect)
    if !isempty(missing_paths)
        error("step " * label * " did not produce: " * join(missing_paths, ", "))
    end
    return
end

# Print & write the timing table (timings_julia.txt in the cwd) ...
function finish(case::String)
    total = sum(last.(TIMINGS); init = 0.0)
    open("timings_julia.txt", "w") do io
        println(io, "# case ", case, " , Julia threads = ", Threads.nthreads(), " , ", now())
        for (l, t) in TIMINGS
            println(io, @sprintf("%-40s %10.2f s", l, t))
        end
        println(io, @sprintf("%-40s %10.2f s", "TOTAL", total))
    end
    println("\n", "="^100)
    for (l, t) in TIMINGS
        println(@sprintf("STEP_TIME %-40s %10.2f s", l, t))
    end
    println(@sprintf("TOTAL_TIME %-39s %10.2f s", case, total))
    println("CASE_DONE ", case, " ", now())
    flush(OUT)
    return
end
