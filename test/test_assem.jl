# idk the actual way to do this sorry
push!(LOAD_PATH, string(@__DIR__, "/../src/"))

using GeometryTypes
using FEA
import IterativeSolvers

include("plottess.jl")

function test()
    circdivs = 20
    degs = LinRange(0, 2*π, circdivs)
    pts = Vector{Point2f0}(undef, 0)
    push!(pts, Point2f0(0,0))
    for R in [0.2, 0.4, 0.6, 0.8]
        circ = map(x -> Point2f0(R*cos(x), R*sin(x)), degs)[1:circdivs-1]
        pts = vcat(pts, circ)
    end
    
    outer_circ = map(x -> Point2f0(1.0*cos(x), 1.0*sin(x)), degs)[1:circdivs-1]
    pts = vcat(outer_circ, pts)
    segments = Vector{FEA.IndexedLineSegment}(undef, 0)
    for (idx, _) in enumerate(outer_circ)
        next_idx = (idx % length(outer_circ)) + 1
        push!(segments, FEA.IndexedLineSegment(idx, next_idx, 1))
    end

    for (idx0, _) in enumerate(outer_circ)
        idx = idx0 + length(outer_circ) + 1
        next_idx0 = (idx0 % length(outer_circ)) + 1
        next_idx = next_idx0 + length(outer_circ) + 1
        push!(segments, FEA.IndexedLineSegment(idx, next_idx, 1))
    end

    pslg = FEA.PSLG(segments)

    f = (x) -> 4.0
    f2 = (x) -> begin 
        x1 = x - [0.5 ; 0.5]
        dist1 = sum(x1 .* x1)
        x2 = x - [-0.5 ; -0.5]
        dist2 = sum(x2 .* x2)
        return -exp(-10 * dist1) - exp(-10 * dist2)
    end
    mesh = FEA.createFEAMesh(pts, pslg)
    A, F = FEA.assemblePoisson(mesh, f)
    FEA.deactivate_external(mesh.tess, mesh.pslg)

    U = IterativeSolvers.cg(A, F)
    display(U)
    println("max min ", maximum(U), " ", minimum(U))

    node_colors = Vector{Float64}(undef, length(mesh.tess.verts))

    for (n_idx, m_idx) in enumerate(mesh.vertex_indexing)
        if m_idx > 0
            node_colors[n_idx] = U[m_idx]
        end
    end

    p1 = plottess(mesh.tess, false, node_colors)
    scene = vbox(
        p1,
        colorlegend(
            p1[end-1],            # get Plot object from Scene
            camera = campixel!, # let vbox decide scene limits
            raw = true,          # no axes, other things as well
            width = (            # make the colorlegend longer so it looks nicer
                30,              # the width
                540              # the height
            )
        )
    )
end

test()
