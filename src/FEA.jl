module FEA

include("geo.jl")
include("triangulation.jl")
include("integrate.jl")

export DelaunayTess2D,
        delaunay2D,
        integratetri,
        ACWTriangle

end
