module FEA

include("geo.jl")
include("triangulation.jl")
include("integrate.jl")

export DelaunayTess2D,
        PSLG,
        IndexedLineSegment,
        ACWTriangle,
        delaunay2D,
        conformingDelaunay2D
        integratetri,
        barycentric

end
