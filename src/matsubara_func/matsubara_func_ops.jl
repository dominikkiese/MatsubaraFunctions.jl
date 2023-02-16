#== 
    Basic operations for MatsubaraFunctions:
        -> addition 
        -> subtraction 
        -> multiplication with scalar 
        -> initialization with scalar or MatsubaraFunction
==#

function check_shape_grid!(
    f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2 :: MatsubaraFunction{GD, SD, DD, Q}
    )  :: Nothing where {GD, SD, DD, Q <: Number}

    @assert data_shape(f1) == data_shape(f2) "Data shapes must be equal"

    for i in 1 : GD 
        @assert type(f1.grids[i]) == type(f2.grids[i]) "Grids must have same particle type" 
        @assert temperature(f1.grids[i]) == temperature(f2.grids[i]) "Grids must have same temperature" 
        @assert index_range(f1.grids[i]) == index_range(f2.grids[i]) "Grids must have same index range" 
    end

    return nothing 
end



function add(
    f1     :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, Q}
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    if checks check_shape_grid!(f1, f2) end
    return MatsubaraFunction(f1.grids, f1.shape, f1.data .+ f2.data; checks)
end

function add!(
    f1     :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, Q}
    ;
    checks :: Bool = true
    )      :: Nothing where {GD, SD, DD, Q <: Number}

    if checks check_shape_grid!(f1, f2) end
    f1.data .+= f2.data

    return nothing 
end



function subtract(
    f1     :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, Q}
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    if checks check_shape_grid!(f1, f2) end
    return MatsubaraFunction(f1.grids, f1.shape, f1.data .- f2.data; checks)
end

function subtract!(
    f1     :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, Q}
    ;
    checks :: Bool = true
    )      :: Nothing where {GD, SD, DD, Q <: Number}

    if checks check_shape_grid!(f1, f2) end
    f1.data .-= f2.data

    return nothing 
end



function mult(
    f      :: MatsubaraFunction{GD, SD, DD, Q},
    val    :: Qp
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

    # type promotion checked by Base
    return MatsubaraFunction(f.grids, f.shape, val .* f.data; checks)
end

function mult!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # type promotion checked by Base
    f.data .*= val 

    return nothing
end
 


function set!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # type promotion checked by Base
    f.data .= val

    return nothing
end

function set!(
    f1     :: MatsubaraFunction{GD, SD, DD, Q},
    f2     :: MatsubaraFunction{GD, SD, DD, Q},
    ; 
    checks :: Bool = true
    )      :: Nothing where {GD, SD, DD, Q <: Number}

    if checks check_shape_grid!(f1, f2) end
    f1.data .= f2.data

    return nothing
end