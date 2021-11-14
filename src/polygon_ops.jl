const Pt{T<:Real} = Union{AbstractVector{T},NTuple{<:Any,T}}
const Poly = AbstractVector{<:Union{NTuple{<:Any,<:Real},AbstractVector{<:Real}}}

function unwrap_point(q::GI.AbstractPoint)
    (q.x, q.y)
end
unwrap_point(q) = q


function _shape_mask(A::AbstractRaster, geom::AbstractVector; shape=:polygon, kw...)
    if shape === :polygon
        _poly_mask!(A, geom; kw...)
    elseif shape == :line
        _line_mask!(A, geom; kw...)
    else
        throw(ArgumentError("shape must be :line or :polygon"))
    end
end

function _poly_mask(A::AbstractRaster, poly::AbstractVector; order=(XDim, YDim))
    missingval isa Nothing && _nomissingerror()
    # We need a tuple of all the dims in `order`
    # We also need the index locus to be the center so we are
    # only selecting cells more than half inside the polygon
    shifted_dims = map(d -> DD.maybeshiftlocus(Center(), d), dims(A))

    poly_bounds = geom_bounds(poly, order)
    array_bounds = bounds(dims(A, order))

    # Only run inpolygon if the polygon bounding box overlaps the array bounding box
    if bounds_overlap(poly_bounds, array_bounds)
        # Get the array as points
        pts = vec(collect(points(shifted_dims; order)))
        if shape === :polygon
            # Use the first column of the output - the points in the polygon,
            # and reshape to match `A`
            inpoly = inpolygon(pts, poly)
            inpoly = BitArray(reshape(view(inpoly, :, 1), size(A)))
        elseif shape === :line
            # Use a tolerance of the average pixel size
            # This is not the most exact metric to use, but we are limited
            # to a single `atol` value.
            meansteps = map(b -> b[2] - b[1], bounds(A)) ./ size(A)
            averagepixel = max(meansteps...)/2
            # Join the line with itself reverse, to form a closed polygon.
            # There must be a better way...
            poly = vcat(poly, reverse(poly))
            inpoly = inpolygon(pts, poly; atol=averagepixel)
            # Take the sedond column of the output - the cells close to the line
            inpoly = BitArray(reshape(view(inpoly, :, 2), size(A)))
        elseif shape === :points

        else
            throw(ArgumentError("`shape` keyword must be :line or :polygon")) 
        end
    else
        inpoly = BitArray(undef, size(A))
        inpoly .= false
    end

    # Rebuild a with the masked values
    return rebuild(A; data=inpoly, missingval=false)
end

function geom_bounds(poly, order)
    nodes = flat_nodes(poly)
    poly_bounds = map(1:length(order)) do i
        extrema((p[i] for p in nodes))
    end
end

function bounds_overlap(a, b)
    axes_cross = map(a, b) do (a_min, a_max), (b_min, b_max)
        if a_max >= b_max
            a_min <= b_max
        else
            a_max >= b_min
        end
    end
    return all(axes_cross)
end

function bbox_overlaps(x, order, poly)
    x_bnds = bounds(x, order)
    p_bbox = GI.bbox(poly)
    # If there is no bbox available just act as if it
    # overlaps and let the checks later on sort it out
    p_bbox isa Nothing && return true
    bbox_dims = length(p_bbox) ÷ 2
    p_bnds = [(p_bbox[i], p_bbox[i+bbox_dims]) for i in 1:bbox_dims]

    return bounds_overlap(p_bnds, x_bnds)
end

function _flat_nodes(A::AbstractVector{<:AbstractVector{<:AbstractVector}})
    Iterators.flatten(map(flat_nodes, A))
end
_flat_nodes(A::AbstractVector{<:AbstractVector{<:AbstractFloat}}) = A
_flat_nodes(A::AbstractVector{<:GI.AbstractGeometry}) = _flat_nodes(map(GI.coordinates, A))


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

# function intersection(l1::Line{T}, l2::Line{T}) where T<:Real
#     a1 = l1.e[1] - l1.s[1]
#     b1 = l1.s.[2] - l1.e[2]
#     c1 = a1 * l1.s[2] + b1 * l1.s[1]

#     a2 = l2.e[1] - l2.s[1]
#     b2 = l2.s[2] - l2.e[2]
#     c2 = a2 * l2.s[2] + b2 * l2.s[1]

#     Δ = a1 * b2 - a2 * b1
#     # If lines are parallel, intersection point will contain infinite values
#     return ((b2 * c1 - b1 * c2) / Δ, (a1 * c2 - a2 * c1) / Δ)
# end

# function _burn_line!(A, linestring)
#     lastpoint = first(linestring)
#     column_num = 1  
#     for point in line[2:end] 
#         if isinside(A, point)

#         end
#         line = Line(lastpoint, point)
#         column_line = Line()
#         crossing = intersection(line, column_line)
#         i = DD.selectindices(dim1, Contains(crossing[1]))
#         ifelse(column_num > 0, column_num)
#     end
# end


struct Point{T}
    x::T
    y::T
end

struct Line{T}
    start::T
    stop::T
end
# All pixels sizes are regular, so we can take shortcuts
function burnline!(A, line, bounds)
    x_scale = (bounds[1][2] - bounds[1][1]) / size(A, 2)
    y_scale = (bounds[2][2] - bounds[2][1]) / size(A, 1)
    raw_x_offset = bounds[1][1]
    raw_y_offset = bounds[2][1]
    # @show raw_x_offset raw_y_offset
    # @show x_scale y_scale
    raw_start, raw_stop = line.start, line.stop # Float
    start = Point((raw_start.x - raw_x_offset)/x_scale, (raw_start.y - raw_y_offset)/y_scale)
    stop = Point((raw_stop.x - raw_x_offset)/x_scale, (raw_stop.y - raw_y_offset)/y_scale)
    # @show start stop raw_start raw_stop
    x, y = floor(Int, start.x) + 1, floor(Int, start.y) + 1 # Int
    # @show x y
    # Grid cells are 1.0 X 1.0.
    diff_x = stop.x - start.x
    diff_y = stop.y - start.y
    # @show diff_x diff_y
    step_x = signbit(diff_x) * -2 + 1
    step_y = signbit(diff_y) * -2 + 1
    # @show step_x step_y
    # Ray/Slope calculations
    # Straight distance to the first vertical/horizontal grid boundaries
    xoffset = stop.x > start.x ? (ceil(start.x) - start.x) : (start.x - floor(start.x))
    yoffset = stop.y > start.y ? (ceil(start.y) - start.y) : (start.y - floor(start.y))
    # @show xoffset yoffset
    # Angle of ray/slope.
    angle = atan(-diff_y, diff_x)
    # How far to move along the ray to cross the first cell boundary.
    max_x = xoffset / cos(angle)
    max_y = yoffset / sin(angle)
    # @show angle max_x max_y
    # How far to move along the ray to move 1 grid cell.
    delta_x = 1.0 / cos(angle)
    delta_y = 1.0 / sin(angle)
    # @show delta_x delta_y
    # Travel one grid cell at a time.
    manhattan_distance = floor(Int, abs(floor(start.x) - floor(stop.x)) + abs(floor(start.y) - floor(stop.y)))
    # @show manhattan_distance
    for t in 0:manhattan_distance
        # checkbounds(Bool, A, x, y) || return A
        @inbounds A[x, y] = true
        # Only move in either X or Y coordinates, not both.
        if abs(max_x) <= abs(max_y)
            max_x += delta_x
            x += step_x
        else
            max_y += delta_y
            y += step_y
        end
        # @show x y
    end
    return A
end

# using ProfileView, Plots
# points = [(253.0, 600.1), (755.0, 2111.0)]
# line = Line(Point(points[1]...), Point(points[2]...))
# cursor = (2, 1)
# A = falses(110, 110)
# @btime burnline!(A, line, cursor, bounds)
# heatmap(300:20:2300, 150:10:1150, A; aspect_ratio=2)
# plot!(reverse.(points))

using BenchmarkTools
A = falses(1100, 1100)
bounds = (100.0, 1200.0), (200.0, 2400.0)
points = [(rand(100.0:0.0001:1200.0), rand(200.0:0.0002:2400.0)) for i in 1:10_000]
lines = [Line(Point(points[i]...), Point(points[i+1]...)) for i in 1:length(points)-1]
@btime Threads.@spawn for line in lines 
    burnline!(A, line, bounds);
end
@btime for line in lines 
    burnline!(A, line, bounds);
end
using Polyester

heatmap(A)
