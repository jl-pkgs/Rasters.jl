"""
    crop(x; to)
    crop(xs...; to)

Crop one or multiple [`AbstractRaster`](@ref) or [`AbstractRasterStack`](@ref) `x`
to match the size of the object `to`, or smallest of any dimensions that are shared.

Otherwise crop to the size of the keyword argument `to`. This can be a
`Tuple` of `Dimension` or any object that will return one from `dims(to)`.

# Keywords

- `to`: the array to crop to. If `to` keyword is passed, the smallest shared
    area of all `x` is used.
- `atol`: the absolute tolerance value to use when comparing the index of x and `to`.
    If `atol` isnt set, `Near` will be used.

# Example

```jldoctest
using Rasters, Plots
evenness = Raster(EarthEnv{HabitatHeterogeneity}, :evenness)
rnge = Raster(EarthEnv{HabitatHeterogeneity}, :range)

# Roughly cut out New Zealand from the evenness raster
nz_bounds = X(Between(165, 180)), Y(Between(-32, -50))
nz_evenness = evenness[nz_bounds...]

# Crop range to match evenness
nz_range = crop(rnge; to=nz_evenness, atol=1e-7)
plot(nz_range)

savefig("build/crop_example.png")
# output
```

![crop]/crop_example.png)

$EXPERIMENTAL
"""
function crop end
function crop(l1::RasterStackOrArray, l2::RasterStackOrArray, ls::RasterStackOrArray...; kw...)
    crop((l1, l2, ls); kw...)
end
function crop(xs::Union{Tuple,NamedTuple}; to=_subsetdims(_shortest, xs), kw...)
    map(l -> crop(l; to, kw...), xs)
end
crop(x::RasterStackOrArray; to, kw...) = _crop_to(x, to; kw...)

# crop `A` to values of dims of `to`
_crop_to(A::RasterStackOrArray, to; kw...) = _crop_to(A, dims(to); kw...)
function _crop_to(x::RasterStackOrArray, to::DimTuple; atol=maybe_eps(to))
    # Take a view over the dim ranges
    _without_mapped_crs(x) do x1
        dimranges = map(to, atol) do d, atol_n
            dx = dims(x1, d)
            l = lookup(dx)
            fi = DD.selectindices(l, At(first(d); atol=atol_n))
            li = DD.selectindices(l, At(last(d); atol=atol_n))
            newindex = fi <= li ? (fi:li) : (li:fi)
            rebuild(dx, newindex)
        end
        view(x1, dimranges...)
    end
end
function _crop_to(x::RasterStackOrArray, geom::GI.AbstractGeometry; 
    order=(XDim, YDim), atol=maybe_eps(to)
)
    # Wrap the bounds in dims so we can reorder them
    wrapped_bounds = map(rebuild, dims(x, order), geom_bounds(geom))
    ordered_bounds = dims(wrapped_bounds, dims(x))
    # Take a view over the bounds
    _without_mapped_crs(x) do x1
        dimranges = map(ordered_bounds) do d, b 
            dx = dims(x1, b)
            l = lookup(dx)
            selector = sampling(l) isa Intervals ? Contains : Near
            fi = DD.selectindices(l, selector(first(d)))
            li = DD.selectindices(l, selector(last(d)))
            newindex = fi <= li ? (fi:li) : (li:fi)
            rebuild(dx, newindex)
        end
        view(x1, dimranges...)
    end
end

"""
    extend(layers::AbstractRaster...)
    extend(layers::Union{NamedTuple,Tuple})
    extend(A::Union{AbstractRaster,AbstractRasterStack}; to)

Extend multiple [`AbstractRaster`](@ref) to match the area covered by all.
A single `AbstractRaster` can be extended by passing the new `dims` tuple
as the second argument.

```jldoctest
using Rasters, Plots
evenness = Raster(EarthEnv{HabitatHeterogeneity}, :evenness)
rnge = Raster(EarthEnv{HabitatHeterogeneity}, :range)

# Roughly cut out South America
sa_bounds = X(Between(-88, -32)), Y(Between(-57, 13))
sa_evenness = evenness[sa_bounds...]

# Extend range to match the whole-world raster
sa_range = extend(sa_evenness; to=rnge)
plot(sa_range)

savefig("build/extend_example.png")
# output

```

![extend](extend_example.png)

$EXPERIMENTAL
"""
function extend end
function extend(l1::RasterStackOrArray, l2::RasterStackOrArray, ls::RasterStackOrArray...; kw...)
    extend((l1, l2, ls...); kw...)
end
function extend(xs::Union{NamedTuple,Tuple}; to=_subsetdims(_longest, xs))
    # Extend all layers to `to`, by default the _largestdims
    return map(l -> extend(l; to), xs)
end
extend(x::RasterStackOrArray; to=dims(x)) = _extend_to(x, to)

_extend_to(x::RasterStackOrArray, to) = _extend_to(x, dims(to))
function _extend_to(A::AbstractRaster, to::GI.AbstractGeometry; order=(XDim, YDim))
    all(map(s -> s isa Regular, span(A, order))) || throw(ArgumentError("All dims must have `Regular` span to be extended with a polygon"))
    wrapped_bounds = map(rebuild, dims(x, order), geom_bounds(geom))
    ordered_bounds = dims(wrapped_bounds, dims(x))
    newdims = map(ordered_bounds) do b
        d = dims(A, b)
        l = lookup(d)
        if order(l) isa ForwardOrdered
            # Use ranges for math because they have TwicePrecision magic
            lowerrange = if first(b) < first(bounds(l))
                # Define a range down to the lowest value, but anchored at the existing value
                first(l):-step(l):first(b)-step(l)
            else
                first(l):step(l):first(l)
            end
            upperrange = if last(b) > last(bounds(l))
                last(l):step(l):last(b)+step(l)
            else
                last(l):step(l):last(l)
            end
            newrange = first(lowerrange):step(l):last(upperrange)
        elseif order(d) isa ReverseOrdered
            lowerrange = if first(b) < first(bounds(l))
                last(l):step(l):first(b)+step(l)
            else
                last(l):step(l):last(l)
            end
            upperrange = if first(b) > last(bounds(l))
                first(l):-step(l):first(b)-step(l)
            else
                first(l):step(l):first(l)
            end
            newrange = last(upperrange):step(l):first(lowerrange)
            newlookup = rebuild(l; data=newrange)
            return rebuild(d, newlookup)
        end
    end
    return _extend_to(A, newdims)
end
function _extend_to(A::AbstractRaster, to::DimTuple)
    sze = map(length, to)
    T = eltype(A)
    # Create a new extended array
    newdata = similar(parent(A), T, sze)
    # Fill it with missing/nodata values
    newdata .= missingval(A)
    # Rebuild the original object with larger data and dims.
    newA = rebuild(A; data=newdata, dims=to)
    # Calculate the range of the old array in the extended array
    ranges = map(dims(A), to) do d, nd
        # TODO use open Interval here
        l = lookup(nd)
        start = DD.selectindices(l, Near(first(d)))
        stop = DD.selectindices(l, Near(last(d)))
        start <= stop ? (start:stop) : (stop:start)
    end
    # Copy the original data to the new array
    # Somehow this is slow from disk?
    newA[ranges...] .= read(A)
    return newA
end
_extend_to(st::AbstractRasterStack, to::Tuple) = map(A -> _extend_to(A, to), st)


# Shared utils

# Get the largest or smallest dimensions in a tuple of AbstractRaster
function _subsetdims(f, layers)
    # Combine the dimensions of all layers
    dims = DD.combinedims(layers...; check=false)
    # Search through all the dimensions choosing the shortest
    alldims = map(DD.dims, layers)
    return map(dims) do d
        matchingdims = map(ds -> DD.dims(ds, (d,)), alldims)
        reduce(matchingdims) do a, b
            _choose(f, a, b)
        end |> first
    end
end

# Choose a dimension from either missing dimension
# (empty Tuple) or a comparison between two 1-Tuples
_choose(f, ::Tuple{}, ::Tuple{}) = ()
_choose(f, ::Tuple{}, (b,)::Tuple) = (b,)
_choose(f, (a,)::Tuple, ::Tuple{}) = (a,)
_choose(f, (a,)::Tuple, (b,)::Tuple) = (f(a, b) ? a : b,)

# Choose the shortest or longest dimension
_shortest(a, b) = length(a) <= length(b)
_longest(a, b) = length(a) >= length(b)
