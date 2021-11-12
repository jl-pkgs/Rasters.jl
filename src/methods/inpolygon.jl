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
    nodes = collect(flat_nodes(poly))
    PolygonInbounds.inpoly2(points, nodes, edges; kw...)
end

function _get_edges(edges, edgenum, poly::AbstractVector{<:GI.AbstractGeometry})
    foldl(poly; init=(edges, edgenum)) do (e, en), p
        _get_edges(e, en, GI.coordinates(p))
    end
end
function _get_edges(edges, edgenum, poly::AbstractVector{<:AbstractVector})
    foldl(poly; init=(edges, edgenum)) do (e, en), p
        _get_edges(e, en, p)
    end
end
# Analyse a single polygon
function _get_edges(edges, edgenum, poly::AbstractVector{<:Union{<:NTuple{<:Any,T},<:AbstractVector{T}}}) where T<:Real
    newedges = Matrix{Int}(undef, length(poly), 2)
    for i in eachindex(poly)[1:end-1]
        newedges[i, 1] = i + edgenum
        newedges[i, 2] = i + edgenum + 1
    end
    newedges[end, 1] = length(poly) + edgenum
    newedges[end, 2] = edgenum + 1

    edges = vcat(edges, newedges)
    return edges, edgenum + length(poly)
end

