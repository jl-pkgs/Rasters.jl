const Pt{T<:Real} = Union{AbstractVector{T},NTuple{<:Any,T}}
const Poly = AbstractVector{<:Union{NTuple{<:Any,<:Real},AbstractVector{<:Real}}}

function unwrap_point(q::GI.AbstractPoint)
    (q.x, q.y)
end
unwrap_point(q) = q

function _shape_mask(A::AbstractRaster, geom::AbstractVector; 
    shape=:polygon, order=(XDim, YDim), kw...
)
    gbounds = geom_bounds(geom, order)
    abounds = bounds(dims(A, order))
    missingval(A) isa Nothing && _nomissingerror()
    # Only mask if the gemoetry bounding box overlaps the array bounding box
    mask = if bounds_overlap(gbounds, abounds)
        if shape === :polygon
            _poly_mask(A, geom; order, kw...)
        elseif shape === :line
            _line_mask(A, geom; order, kw...)
        elseif shape === :point
            _point_mask(A, geom; order, kw...)
        else
            throw(ArgumentError("shape must be :line or :polygon"))
        end
    else
        # Otherwise return all false
        falses(size(A))
    end
    # Rebuild a with the masked values
    return rebuild(A; data=mask, missingval=false)
end

function _poly_mask(A::AbstractRaster, poly::AbstractVector; order)
    # We need a tuple of all the dims in `order`
    # We also need the index locus to be the center so we are
    # only selecting cells more than half inside the polygon
    shifted_dims = map(d -> DD.maybeshiftlocus(Center(), d), dims(A))
    # Get the array as points
    pts = vec(collect(points(shifted_dims; order)))
    # Use the first column of the output - the points in the polygon,
    # and reshape to match `A`
    inpoly = inpolygon(pts, poly)
    return BitArray(reshape(view(inpoly, :, 1), size(A)))
end

function _line_mask(A::AbstractRaster, lines::AbstractVector; order)
    # Use a tolerance of the average pixel size
    # This is not the most exact metric to use, but we are limited
    # to a single `atol` value.
    # meansteps = map(b -> b[2] - b[1], bounds(A)) ./ size(A)
    # averagepixel = max(meansteps...)/2
    # # Join the line with itself reverse, to form a closed polygon.
    # # There must be a better way...
    # poly = vcat(poly, reverse(poly))
    # inpoly = inpolygon(pts, poly; atol=averagepixel)
    # # Take the sedond column of the output - the cells close to the line
    # return BitArray(reshape(view(inpoly, :, 2), size(A)))
    B = falses(size(A))
    P = NamedTuple{(:x,:y)}
    for i in eachindex(lines)[1:end-1]
        start = ntuple(n -> lines[i][n], length(order)) 
        stop = ntuple(n -> lines[i + 1][n], length(order)) 
        line = Line(P(start), P(stop))
        _burn_line!(B, line, bounds(A))
    end
    return B
end

function _point_mask(A::AbstractRaster, points::AbstractVector; order)
    # Just find which pixels contian the points, and set them to true
    data = falses(size(A))
    for point in points 
        selectors = map(dims(A, order), ntuple(i -> i, length(order))) do d, i
            rebuild(d, Contains(point[i]))
        end
        I = selectindices(A, selectors...)
        data[I...] = true
    end
end

function geom_bounds(poly, order)
    nodes = _flat_nodes(poly)
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
    Iterators.flatten(map(_flat_nodes, A))
end
_flat_nodes(A::AbstractVector{<:AbstractVector{<:AbstractFloat}}) = A
_flat_nodes(A::AbstractVector{<:Tuple}) = A
_flat_nodes(A::AbstractVector{<:GI.AbstractGeometry}) = _flat_nodes(map(GI.coordinates, A))


function _to_edges(poly)
    edgenum = 0
    edges = Matrix{Int}(undef, 0, 2)
    _to_edges!(edges, edgenum, poly)
end

function _to_edges!(edges, edgenum, poly::AbstractVector{<:GI.AbstractGeometry})
    foldl(poly; init=(edges, edgenum)) do (e, en), p
        _to_edges!(e, en, GI.coordinates(p))
    end
end
function _to_edges!(edges, edgenum, poly::AbstractVector{<:AbstractVector})
    foldl(poly; init=(edges, edgenum)) do (e, en), p
        _to_edges!(e, en, p)
    end
end
# Analyse a single polygon
function _to_edges!(edges, edgenum, poly::AbstractVector{<:Union{<:NTuple{<:Any,T},<:AbstractVector{T}}}) where T<:Real
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

struct Line{T}
    start::T
    stop::T
end
# All pixels sizes are regular, so we can take shortcuts
function _burn_line!(A, line, bounds)
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

# using BenchmarkTools, Plots
# A = falses(1100, 1100)
# bounds = (100.0, 1200.0), (200.0, 2400.0)
# points = [(rand(100.0:0.0001:1200.0), rand(200.0:0.0002:2400.0)) for i in 1:1000_000]
# lines = [Line(Point(points[i]...), Point(points[i+1]...)) for i in 1:length(points)-1]
# @time for line in lines 
#     burnline!(A, line, bounds);
# end
# heatmap(A)

# using Polyester, ProfileView
# @time for line in lines 
#     burnline!(A, line, bounds);
# end
# @btime @batch minbatch=5000 for line in lines 
#     burnline!(A, line, bounds)
# end

