# set environment variables for mpi usage
ENV["JULIA_EXCLUSIVE"]   = "0"
ENV["UCX_ERROR_SIGNALS"] = "SIGILL,SIGBUS,SIGFPE"

"""
    function mpi_comm() :: MPI.Comm

Return the MPI communicator
"""
mpi_comm() :: MPI.Comm = MPI.COMM_WORLD

"""
    function mpi_rank() :: Int64

Return the current MPI rank
"""
mpi_rank() :: Int64 = MPI.Comm_rank(mpi_comm())

"""
    function mpi_size() :: Int64

Return the size of the MPI communicator 
"""
mpi_size() :: Int64 = MPI.Comm_size(mpi_comm())

# mpi command for simple loop parallelization (main is busy) 
"""
    function mpi_split(
        r :: UnitRange{Int64}
        ) :: UnitRange{Int64}

Splits `UnitRange` evenly among available MPI ranks (including main). 
Can, for example, be used to parallelize loops.
"""
function mpi_split(
    r :: UnitRange{Int64}
    ) :: UnitRange{Int64}

    # if there is only one rank return the full range 
    if mpi_size() == 1 
        return r

    # if there are more ranks than possible jobs, assign one job per 
    # available rank and leave the rest idle
    elseif mpi_size() >= length(r)
        if mpi_rank() + 1 <= length(r)
            return r[mpi_rank() + 1 : mpi_rank() + 1]
        else 
            return 1 : 0 
        end 

    # if there are more jobs than ranks, split the jobs as evenly as 
    # possible among the ranks 
    else 
        div, rem = divrem(length(r), mpi_size())
        
        if mpi_rank() + 1 <= rem 
            return r[mpi_rank() * (div + 1) + 1 : mpi_rank() * (div + 1) + (div + 1)]
        else
            return r[end - (mpi_rank() + 1 - rem) * div + 1 : end - (mpi_rank() + 1 - rem) * div + div]
        end
    end
end

"""
    function mpi_allreduce!(
        f :: MeshFunction{{MD, SD, DD, Q, AT}}
        ) :: Nothing where {GD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Inplace MPI reduction (+) for MeshFunction
"""
function mpi_allreduce!(
    f :: MeshFunction{GD, SD, DD, Q, AT}
    ) :: Nothing where {GD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    MPI.Allreduce!(f.data, +, mpi_comm())
    return nothing
end

"""
    mpi_ismain() :: Bool

Returns true for MPI rank 0
"""
mpi_ismain() :: Bool = mpi_rank() == 0

"""
    function mpi_println(
        s :: String
        ) :: Nothing

Print string `s` on MPI rank 0
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

"""
    mpi_info() :: Nothing

Print information about available resources
"""
function mpi_info() :: Nothing
    mpi_println("NOTE: Running in MPI environment with $(mpi_size()) rank(s) and $(Threads.nthreads()) thread(s) per rank")
    return nothing 
end
 
"""
    mpi_barrier() :: Nothing 

Place synchronization barrier for MPI ranks
"""
function mpi_barrier() :: Nothing 
    MPI.Barrier(mpi_comm())
    return nothing 
end

#----------------------------------------------------------------------------------------------#

export 
    mpi_comm,
    mpi_rank,
    mpi_size,
    mpi_split,
    mpi_allreduce!,
    mpi_ismain,
    mpi_println,
    mpi_info,
    mpi_barrier