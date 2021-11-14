"""
    inpolygon(points, poly)

Check if a point or `Vector` of points is inside a polygon.

This algorithm is very efficient for many points, less so a single point.

# Arguments

- `points`: an `AbstractVector` or a `Tuple` or `Real`, Or a `Vector` of these.
- `poly`: an `AbstractVector` or nested `AbstractVector` with an inner
    `AbstractVector` or `Tuple` of `Real`. It can also be a `GeoInterface.AbstractGeometry`.

Returns a `Bool` or `BitVector{Bool}
"""
function inpolygon end
function inpolygon(point::Union{NTuple{<:Any,At},Pt}, poly::GI.AbstractGeometry; kw...)
    inpolygon(point, GI.coordinates(poly); kw...)
end
function inpolygon(points::AbstractVector, poly::GI.AbstractGeometry; kw...)
    inpolygon(points, GI.coordinates(poly); kw...)
end
inpolygon(point::AbstractVector{<:Real}, poly::AbstractVector; kw...) = inpoly([point], poly; kw...)
inpolygon(point::Tuple, poly::AbstractVector; kw...) = inpoly([point], poly; kw...)
function inpolygon(points::AbstractVector, poly::AbstractVector; kw...)
    edges = Matrix{Int}(undef, 0, 2)
    edgenum = 0
    edges, _ = _get_edges(edges, edgenum, poly)
    nodes = collect(_flat_nodes(poly))
    PolygonInbounds.inpoly2(points, nodes, edges; kw...)
end
