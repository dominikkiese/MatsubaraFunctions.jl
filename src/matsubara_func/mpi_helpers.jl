# set environment variables for mpi usage
ENV["JULIA_EXCLUSIVE"]   = "0"
ENV["UCX_ERROR_SIGNALS"] = "SIGILL,SIGBUS,SIGFPE"

# mpi command for simple loop parallelization (main is busy) 
"""
    function mpi_split(
        r :: UnitRange{Int64}
        ) :: UnitRange{Int64}

Splits UnitRange evenly among available MPI ranks (including main). 
Can, for example, be used to parallelize loops.
"""
function mpi_split(
    r :: UnitRange{Int64}
    ) :: UnitRange{Int64}

    comm     = MPI.COMM_WORLD 
    mpi_size = MPI.Comm_size(comm)
    mpi_rank = MPI.Comm_rank(comm)
    r_length = length(r)

    # if there is only one rank return the full range 
    if mpi_size == 1 
        return r

    # if there are more ranks than possible jobs, assign one job per 
    # available rank and leave the rest idle
    elseif mpi_size >= r_length
        if mpi_rank + 1 <= r_length
            return r[mpi_rank + 1 : mpi_rank + 1]
        else 
            return 1 : 0 
        end 

    # if there are more jobs than ranks, split the jobs as evenly as 
    # possible among the ranks 
    else 
        div, rem = divrem(r_length, mpi_size)
        
        if mpi_rank + 1 <= rem 
            return r[mpi_rank * (div + 1) + 1 : mpi_rank * (div + 1) + (div + 1)]
        else
            return r[end - (mpi_rank + 1 - rem) * div + 1 : end - (mpi_rank + 1 - rem) * div + div]
        end
    end
end

# mpi command for simple inplace reduction 
"""
    function mpi_allreduce!(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: Nothing where {GD, SD, DD, Q <: Number}

Inplace MPI reduction (+) for MatsubaraFunction
"""
function mpi_allreduce!(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Nothing where {GD, SD, DD, Q <: Number}

    MPI.Allreduce!(f.data, +, MPI.COMM_WORLD)

    return nothing
end

# mpi command to check if root 
"""
    mpi_ismain() :: Bool

Returns true for MPI rank 0
"""
mpi_ismain() :: Bool = MPI.Comm_rank(MPI.COMM_WORLD) == 0

# mpi command to print from main 
"""
    function mpi_println(
        s :: String
        ) :: Nothing

Print string s on MPI rank 0
"""
function mpi_println(
    s :: String
    ) :: Nothing

    if mpi_ismain()
        println(s)
        flush(stdout)
    end 

    return nothing 
end

# mpi command to print available resources 
"""
    mpi_info() :: Nothing

Print information about available resources
"""
function mpi_info() :: Nothing
    mpi_println("NOTE: Running in MPI environment with $(MPI.Comm_size(MPI.COMM_WORLD)) rank(s) and $(Threads.nthreads()) thread(s) per rank")
    return nothing 
end