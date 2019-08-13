using GeometryTypes
using SparseArrays

function _create_node_func(tess::DelaunayTess2D{T}, pt::T) where T <: Point2
    start_tri = nothing
    for t in tess.faces
        if t.active
            if intriangle(tess, pt)
                start_tri = t
            end
        end
    end
    @assert start_tri != nothing

    tris = [start_tri]
    idxs = []

    while true
        t = tris[end]
        if t.a == pt
            n   = tess.faces[t.nb]
            idx = TRI_NEIGHBOR_A
        elseif t.b == pt
            n   = tess.faces[t.nc]
            idx = TRI_NEIGHBOR_B
        elseif t.c == pt
            n   = tess.faces[t.na]
            idx = TRI_NEIGHBOR_C
        else
            @assert false
        end
        push!(idxs, idx)
        if n != tris[1]
            push!(tris, n)
        else
            break
        end
    end

    # create 'tent' function for the node using barycentric coordinates
    f = (x) -> begin
        for i in 1:length(tris)
            t = tris[i]
            if intriangle(t, x)
                return barycentric(t, x)[idxs[i]]
            end
        end
        return 0.0
    end
end

function assemblePoisson(V::AbstractVector{T}, tess::DelaunayTess2D{T},
    f::Function) where T <: Point2

    # TODO
    n_internal_nodes = 10
    A    = spzeros(n_internal_nodes, n_internal_nodes)
    load = zeros(n_internal_nodes)
end
