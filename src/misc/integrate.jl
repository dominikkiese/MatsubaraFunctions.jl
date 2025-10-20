function integrate(f :: Function, m :: MT) where {MT <: AbstractMesh}
    return volume_element(m) * sum(f, m)
end

export integrate