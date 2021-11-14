"""
    rasterize(data; kw...)
    rasterize(points, values; kw...)

Rasterize the points and values in `data`, or the `points` and `values` objects,
into the [`Raster`](@ref) or [`RasterStack`](@ref) `x`. 

# Arguments

- `data`: a Tables.jl compatible object containing points and values or a
    polygon - an GeoInterface.jl `AbstractGeometry`, or a nested `Vector` of `Vectors`.
- `points`: A `Vector` or nested `Vectors` holding `Vector` or `Tuple` of `Real`
- `values` A `Vector` of values to be written to a `Raster`, or a Vector of `NamedTupled`
    to write to a `RasterStack`.

# Keywords

These are detected automatically from `A` and `data` where possible.

- `to`: a `Raster` or `RasterStack` to rasterize to.
- `order`: A `Tuple` of pairs `Dim => Symbol` for the keys in the data that match
    the dimension.
- `value`: A `Tuple` of `Symbol` for the keys in the data that provide
    values to add to `A`.
- `fill`: the value to fill a polygon with, if `data` is a polygon. 
- `atol`: an absolute tolerance for rasterizing to dimensions with `Points` sampling.

# Example

```jldoctest
using Rasters, Plots, Dates, Shapefile, GeoInterface, Downloads
using Rasters.LookupArrays

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "boundary_lines.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# Load the shapes for denmark
indonesia_border = Shapefile.Handle(shapefile_name).shapes[1]

# Make an empty EPSG 4326 projected Raster of the area of Indonesia
dimz = Y(-15.0:0.1:10.9; mode=Projected(; sampling=Intervals(Start()), crs=EPSG(4326))), 
       X(90.0:0.1:145; mode=Projected(; sampling=Intervals(Start()), crs=EPSG(4326)))
A = Raster(zeros(UInt16, dimz); missingval=0)

# Rasterize each island with a different number
for (i, shp) in enumerate(coordinates(indonesia_border))
    rasterize!(A, shp; fill=i, order=(X, Y))
end

# And plot
p = plot(A; color=:spring)
plot!(p, indonesia_border; fillalpha=0, linewidth=0.7)
savefig("build/indonesia_rasterized.png")

# output

```

![rasterize](indonesia_rasterized.png)

$EXPERIMENTAL
"""
rasterize(args...; to, kw...) = _rasterize(to, args...; kw...)

function _rasterize(to::DimTuple, args...; kw...)
    A = Raster(zeros(Union{Float64,Missing}, to))
    return rasterize!(A, args...; kw...)
end
function _rasterize(to::AbstractRaster, args...; kw...)
    A = similar(to) .= missingval(to)
    return rasterize!(A, args...; kw...)
end
function _rasterize(to::AbstractRasterStack, args...; kw...)
    st = map(to) do A
        similar(A) .= missingval(A)
    end
    return rasterize!(st, args...; kw...)
end


"""
    rasterize!(x, data; order, name, atol)
    rasterize!(x, points, values; order, atol)

Rasterize the points and values in `data`, or the `points` and `values` objects,
into the [`Raster`](@ref) or [`RasterStack`](@ref) `x`.

# Arguments

- `x`: a `Raster` or `RasterStack` to rasterize to.
- `data`: a Tables.jl compatible object containing points and values or a
    polygon - an GeoInterface.jl `AbstractGeometry`, or a nested `Vector` of `Vectors`.
- `points`: A `Vector` or nested `Vector` holding `Vector` or `Tuple` of `Real`
- `values` A `Vector` of values to be written when `x` is a `Raster`, or a Vector of
    `NamedTupled` to write when `x` is a `RasterStack`.

# Keywords

These are detected automatically from `A` and `data` where possible.

- `point`: A `Tuple` of pairs `Dim => Symbol` for the keys in the data that match
    the dimension.
- `value`: A `Tuple` of `Symbol` for the keys in the data that provide
    values to add to `A`.
- `fill`: the value to fill a polygon with, if `data` is a polygon. 
- `atol`: an absolute tolerance for rasterizing to dimensions with `Points` sampling.

# Example

Rasterize a shapefile for denmark and plot, with a border.

```jldoctest
using Rasters, Plots, Dates, Shapefile, Downloads
using Rasters.LookupArrays

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "boundary_lines.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# Loade the shapes for china
china_border = Shapefile.Handle(shapefile_name).shapes[10]

# Make an empty EPSG 4326 projected Raster of the China area
dimz = Y(Projected(15.0:0.1:55.0; sampling=Intervals(Start()), crs=EPSG(4326))), 
       X(Projected(70.0:0.1:140; sampling=Intervals(Start()), crs=EPSG(4326)))
A = Raster(zeros(UInt8, dimz); missingval=0)

# Rasterize the border polygon 
rasterize!(A, china_border; fill=1, order=(X, Y))

# And plot
p = plot(A; color=:spring)
plot!(p, china_border; fillalpha=0, linewidth=0.6)
savefig("build/china_rasterized.png")

# output

```

![rasterize](china_rasterized.png)

$EXPERIMENTAL
"""
function rasterize!(A::AbstractRaster, points, values;
    order=(XDim, YDim, ZDim), atol=nothing
)
    isdisk(A) && _warn_disk(rasterize)
    ordered_dims = dims(A, ntuple(i -> order[i], length(first(points))))
    _without_mapped_crs(A) do A1
        map(points, values) do p, v
            any(map(ismissing, p)) && return nothing
            selectors = map((d, x) -> _at_or_contains(d, x, atol), ordered_dims, p)
            if length(selectors) == length(dims(A))
                A1[selectors...] = v
            else
                A1[selectors...] .= v
            end
            return nothing
        end
    end
    return A
end
function rasterize!(st::AbstractRasterStack, points, values;
    order=(XDim, YDim, ZDim), atol=nothing
)
    isdisk(first(st)) && _warn_disk(rasterize!)
    ordered_dims = dims(st, order)
    _without_mapped_crs(st) do st1
        map(points, values) do p, v
            any(map(ismissing, p)) && return nothing
            selectors = map((d, x) -> _at_or_contains(d, x, atol), ordered_dims, p)
            map(Base.values(st1), v) do A, v_n
                if length(selectors) == length(dims(A))
                    A[selectors...] = v_n
                else
                    A[selectors...] .= v_n
                end
            end
            return nothing
        end
    end
    return st
end
function rasterize!(A::AbstractRaster, data; kw...)
    if Tables.istable(data)
        _rasterize_table!(A, data; kw...)
    else
        _rasterize_geometry!(A, data; kw...)
    end
end

function _rasterize_table!(A::AbstractRaster, data;
    order=_auto_pointcols(A, data),
    name=nothing, kw...
)
    istable(data)

    if name isa Nothing 
        name = first(_not_a_dimcol(data, order))
    end
    isdisk(data) && _warn_disk(rasterize)
    ordered_dims = map(p -> DD.basetypeof(p[1])(p[2]), order)
    ordered_keys = map(last, order)
    points = (map(k -> r[k], ordered_keys) for r in Tables.rows(data))
    if name isa Symbol
        values = (r[name] for r in Tables.rows(data))
    elseif value isa Tuple
        values = (r[first(name)] for r in Tables.rows(data))
    end
    return rasterize!(A, points, values; order=ordered_dims, kw...)
end
function _rasterize_table!(st::AbstractRasterStack, data;
    point=_auto_pointcols(st, data), name=_not_a_dimcol(data, point), kw...
)
    isdisk(data) && _warn_disk(rasterize!)
    point_dims = map(p ->  DD.basetypeof(p[1])(p[2]), point)
    order = map(last, point)
    points = (map(pk -> r[pk], order) for r in Tables.rows(data))
    if name isa Symbol
        values = (r[name] for r in Tables.rows(data))
    elseif name isa Tuple
        values = (map(vk -> r[vk], name) for r in Tables.rows(data))
    end
    return rasterize!(st, points, values; order=point_dims, kw...)
end

function _rasterize_geometry!(x, coll::GI.AbstractFeatureCollection; kw...)
    foreach(coll.features) do feature
        _rasterize_geometry!(x, feature; kw...)
    end
    return x
end
function _rasterize_geometry!(x, feature::GI.AbstractFeature; kw...)
    _rasterize_geometry!(x, feature.geometry; kw...)
end
function _rasterize_geometry!(x, poly::GI.AbstractGeometry;
    order=(XDim, YDim, ZDim), kw...
)
    shape = if poly isa GI.AbstractPolygon
        :polygon
    elseif poly isa GI.AbstractLineString
        :line
    else
        throw(ArgumentError("`shape` keyword must be :line or :polygon")) 
    end

    if bbox_overlaps(x, order, poly)
        _rasterize_geometry!(x, GI.coordinates(poly); order, shape, kw...)
    end
    return x
end
function _rasterize_geometry!(x::RasterStackOrArray, poly::AbstractVector{<:AbstractVector};
    fill, order=(XDim, YDim, ZDim), kw... 
)
    ordered_dims = dims(st, order)
    B = _poly_mask(first(st), poly; order=ordered_dims)
    return _fill!(x, B, poly, fill)
end

function _fill!(A::AbstractRasterStack, B, poly, fill)
    map((a, f) -> _fill!(a, B, poly, f), st, fill)
    return A
end
function _fill!(A, B, poly, fill)
    broadcast!(A, A, B) do a, b
        b ? (fill isa Function ? fill(a) : fill) : a
    end
end

function _at_or_contains(d, v, atol)
    selector = sampling(d) isa Intervals ? Contains(v) : At(v; atol=atol)
    DD.basetypeof(d)(selector)
end

function _auto_pointcols(A, data)
    names = Tables.columnnames(data)
    if names == ()
        names = keys(first(Tables.rows(data)))
    end
    Tuple(DD.basedims(d) => DD.dim2key(d) for d in dims(A) if DD.dim2key(d) in names)
end
