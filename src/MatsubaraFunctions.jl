module MatsubaraFunctions

    using PrecompileTools

    @recompile_invalidations begin
        using LinearAlgebra 
        using StaticArrays
        using FunctionWrappers
        using MPI
        using HDF5
        using Documenter
    end
    
    include("matsubara_freq/matsubara_freq.jl")
    include("matsubara_grid/matsubara_grid.jl")
    include("matsubara_func/matsubara_func.jl")

    @compile_workload begin
        fg = MatsubaraGrid(1.0, 10, Fermion)
        bg = MatsubaraGrid(1.0, 10, Boson)

        # typical MatsubaraFunction constructors 
        f1D_c = MatsubaraFunction(fg, 1)
        f1D_r = MatsubaraFunction(fg, 1, Float64)

        f2D_c = MatsubaraFunction((bg, fg), 1)
        f2D_r = MatsubaraFunction((bg, fg), 1, Float64)

        f3D_c = MatsubaraFunction((bg, fg, fg), 1)
        f3D_r = MatsubaraFunction((bg, fg, fg), 1, Float64)

        # typical MatsubaraFunction accessors 
        f1D_c[fg[1]]; f1D_c(fg[1]); f1D_c(value(fg[1]))
        f1D_r[fg[1]]; f1D_r(fg[1]); f1D_r(value(fg[1]))

        f2D_c[bg[1], fg[1]]; f2D_c(bg[1], fg[1]); f2D_c(value(bg[1]), value(fg[1]))
        f2D_r[bg[1], fg[1]]; f2D_r(bg[1], fg[1]); f2D_r(value(bg[1]), value(fg[1]))

        f3D_c[bg[1], fg[1], fg[1]]; f3D_c(bg[1], fg[1], fg[1]); f3D_c(value(bg[1]), value(fg[1]), value(fg[1]))
        f3D_r[bg[1], fg[1], fg[1]]; f3D_r(bg[1], fg[1], fg[1]); f3D_r(value(bg[1]), value(fg[1]), value(fg[1]))

        # typical MatsubaraFunction slices and views 
        f1D_c[:]; view(f1D_c, :)
        f1D_r[:]; view(f1D_r, :)

        f2D_c[:, fg[1]]; view(f2D_c, :, fg[1])
        f2D_r[:, fg[1]]; view(f2D_r, :, fg[1])
        f2D_c[bg[1], :]; view(f2D_c, bg[1], :)
        f2D_r[bg[1], :]; view(f2D_r, bg[1], :)

        f3D_c[:, fg[1], fg[1]]; view(f3D_c, :, fg[1], fg[1])
        f3D_r[:, fg[1], fg[1]]; view(f3D_r, :, fg[1], fg[1])
        f3D_c[bg[1], :, fg[1]]; view(f3D_c, bg[1], :, fg[1])
        f3D_r[bg[1], :, fg[1]]; view(f3D_r, bg[1], :, fg[1])
        f3D_c[bg[1], fg[1], :]; view(f3D_c, bg[1], fg[1], :)
        f3D_r[bg[1], fg[1], :]; view(f3D_r, bg[1], fg[1], :)
    end

    export
        # matsubara_freq.jl 
        AbstractParticle,
        Fermion,
        Boson,
        MatsubaraFrequency,
        temperature,
        value, 
        index,
        type,

        # matsubara_grid.jl
        MatsubaraGrid, 
        first_index,
        last_index,
        index_range,
        is_inbounds,
        N,
        first_value,
        last_value,
        value_range,
        indices,
        info,
        save_matsubara_grid!, 
        load_matsubara_grid,

        # matsubara_func.jl
        MatsubaraFunction,
        grids, 
        grids_shape,
        shape, 
        data_shape,
        absmax,
        argmax,
        mpi_comm,
        mpi_rank,
        mpi_size,
        mpi_split,
        mpi_allreduce!,
        mpi_ismain,
        mpi_println,
        mpi_info,
        mpi_barrier,
        add,
        add!, 
        subtract,
        subtract!, 
        mult,
        mult!,
        set!,
        flatten,
        flatten!,
        unflatten!,
        LinearIndex,
        to_Matsubara,
        upper_tail_moments,
        lower_tail_moments,
        sum_me,
        density,
        MatsubaraOperation,
        sgn,
        con,
        MatsubaraSymmetry,
        MatsubaraSymmetryGroup,
        MatsubaraInitFunction,
        save_matsubara_function!, 
        load_matsubara_function,
        save_matsubara_symmetry_group!,
        load_matsubara_symmetry_group
end