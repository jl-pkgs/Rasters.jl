"""
    mosaic(f, regions...; dims, missingval, atol)
    mosaic(f, regions::Tuple; dims, missingval, atol)

Combine `layer`s using the function `f`, (e.g. `mean`, `sum`,
`first` or `last`) where values from `regions` overlap.

# Keywords

- `dims`: The dimesions of `regions` to mosaic over, `(XDim, YDim)` by default.
    If dims contains an index it will be ignored, but this may change in future.
- `missingval`: Fills empty areas, and defualts to the
    `missingval` of the first layer.
- `atol`: Absolute tolerance for comparison between index values.
    This is often required due to minor differences in range values
    due to floating point error. It is not applied to non-float dimensions.
    A tuple of tolerances may be passed, matching the dimension order.

If your mosaic has has apparent line errors, increase the `atol` value.

# Example

Here we cut out australia and africa from a stack, and join them with `mosaic`.

```jldoctest
using Rasters, Plots
st = RasterStack(WorldClim{Climate}; month=1);

africa = st[X(Between(-20.0, 60.0)), Y(Between(35.0, -40.0))]
a = plot(africa)

aus = st[X(Between(100.0, 160.0)), Y(Between(-10.0, -50.0))]
b = plot(aus)

# Combine with mosaic
mos = mosaic(first, aus, africa)
c = plot(mos)

savefig(a, "build/mosaic_example_africa.png")
savefig(b, "build/mosaic_example_aus.png")
savefig(c, "build/mosaic_example_combined.png")
# output

```

### Individual continents

![arica](mosaic_example_africa.png)

![aus](mosaic_example_aus.png)

### Mosaic of continents

![mosaic](mosaic_example_combined.png)

$EXPERIMENTAL
"""
mosaic(f::Function, regions...; kw...) = mosaic(f, regions; kw...)
function mosaic(f::Function, regions::Tuple{<:AbstractRaster,Vararg};
    missingval=missingval(first(regions)), filename=nothing, kw...
)
    missingval isa Nothing && throw(ArgumentError("Layers have no missingval, so pass a `missingval` keyword explicitly"))
    T = Base.promote_type(typeof(missingval), Base.promote_eltype(regions...))
    dims = _mosaic(map(DD.dims, regions))
    data = if filename isa Nothing
        Array{T,length(dims)}(undef, map(length, dims))
    else
        l1 = first(regions)
        create(filename, T, dims; name=name(l1), missingval, metadata=metadata(l1))
        parent(Raster(filename))
    end
    A = rebuild(first(regions); data, dims, missingval)
    open(A; write=true) do a
        mosaic!(f, a, regions; missingval, kw...)
    end
    return A
end
function mosaic(f::Function, regions::Tuple{<:AbstractRasterStack,Vararg}; kw...)
    map(regions...) do A...
        mosaic(f, A...; kw...)
    end
end

"""
    mosaic!(f, x, regions...; missingval, atol)
    mosaic!(f, x, regions::Tuple; missingval, atol)

Combine `regions`s in `x` using the function `f`.

# Arguments

- `f` a function (e.g. `mean`, `sum`, `first` or `last`) that is applied to
    values where `regions` overlap.
- `x`: A `Raster` or `RasterStack`. May be a an opened disk-based `Raster`,
    the result will be written to disk.
    slow read speed with the current algorithm
- `regions`: source objects to be joined. These should be memory-backed
    (use `read` first), or may experience poor performance. If all objects have
    the same extent, `mosaic` is simply a merge.

# Keywords

- `missingval`: Fills empty areas, and defualts to the `missingval/
    of the first layer.
- `atol`: Absolute tolerance for comparison between index values.
    This is often required due to minor differences in range values
    due to floating point error. It is not applied to non-float dimensions.
    A tuple of tolerances may be passed, matching the dimension order.

# Example

Cut out Australia and Africa stacks, then combined them
into a single stack.

```jldoctest
using Rasters, Statistics, Plots
st = read(RasterStack(WorldClim{Climate}; month=1))
aus = st[X(Between(100.0, 160.0)), Y(Between(-10.0, -50.0))]
africa = st[X(Between(-20.0, 60.0)), Y(Between(35.0, -40.0))]
mosaic!(first, st, aus, africa)
plot(st)
savefig("build/mosaic_bang_example.png")
# output

```

![mosaic](mosaic_bang_example.png)

$EXPERIMENTAL
"""
function mosaic!(f::Function, A::AbstractRaster{T}, regions;
    missingval=missingval(A), atol=_default_atol(T)
) where T
    _without_mapped_crs(A) do A1
        broadcast!(A1, DimKeys(A1; atol)) do ds
            # Get all the regions that have this point
            ls = foldl(regions; init=()) do acc, l
                if DD.hasselection(l, ds)
                    v = l[ds...]
                    (acc..., l)
                else
                    acc
                end
            end
            values = foldl(ls; init=()) do acc, l
                v = l[ds...]
                if isnothing(Rasters.missingval(l))
                    (acc..., v)
                elseif ismissing(Rasters.missingval(l))
                    ismissing(v) ? acc : (acc..., v)
                else
                    v === Rasters.missingval(l) ? acc : (acc..., v)
                end
            end
            if length(values) === 0
                missingval
            else
                f(values)
            end
        end
    end
    return A
end
function mosaic!(f::Function, st::AbstractRasterStack, regions; kw...)
    map(st, regions...) do A, r...
        mosaic!(f, A, r; kw...)
    end
end
mosaic!(f::Function, x, regions...; kw...) = mosaic!(f, x, regions; kw...)

_mosaic(alldims::Tuple{<:DimTuple,Vararg{<:DimTuple}}) = map(_mosaic, alldims...)
function _mosaic(dims::Dimension...)
    map(dims) do d
        DD.comparedims(first(dims), d; val=false, length=false, lookup=true)
    end
    return rebuild(first(dims), _mosaic(lookup(dims)))
end
_mosaic(lookups::LookupArrayTuple) = _mosaic(first(lookups), lookups)
function _mosaic(lookup::Categorical, lookups::LookupArrayTuple)
    newindex = union(lookups...)
    if order isa ForwardOrdered
        newindex = sort(newindex; order=LA.ordering(order(lookup)))
    end
    return rebuild(lookup; data=newindex)
end
function _mosaic(lookup::AbstractSampled, lookups::LookupArrayTuple)
    order(lookup) isa Unordered && throw(ArgumentError("Cant mozaic an Unordered lookup"))
    return _mosaic(span(lookup), lookup, lookups)
end
function _mosaic(span::Regular, lookup::AbstractSampled, lookups::LookupArrayTuple)
    newindex = if order(lookup) isa ForwardOrdered
        mi = minimum(map(first, lookups))
        ma = maximum(map(last, lookups))
        if mi isa AbstractFloat
            # Handle slight range erorrs to make sure
            # we dont drop one step of the range
            mi:step(span):ma + 2eps(ma)
        else
            mi:step(span):ma
        end
    else
        mi = minimum(map(last, lookups))
        ma = maximum(map(first, lookups))
        if mi isa AbstractFloat
            ma:step(span):mi - 2eps(mi)
        else
            ma:step(span):mi
        end
    end
    return rebuild(lookup; data=newindex)
end

function _mosaic(::Irregular, lookup::AbstractSampled, lookups::LookupArrayTuple)
    newindex = sort(union(map(parent, lookups)...); order=LA.ordering(order(lookup)))
    return rebuild(lookup; data=newindex)
end
function _mosaic(span::Explicit, lookup::AbstractSampled, lookups::LookupArrayTuple)
    # TODO make this less fragile to floating point innaccuracy
    newindex = sort(union(map(parent, lookups)...); order=LA.ordering(order(lookup)))
    bounds = map(val ∘ DD.span, lookups)
    lower = map(b -> view(b, 1, :), bounds)
    upper = map(b -> view(b, 2, :), bounds)
    newlower = sort(union(lower...); order=LA.ordering(order(lookup)))
    newupper = sort(union(upper...); order=LA.ordering(order(lookup)))
    newbounds = vcat(permutedims(newlower), permutedims(newupper))
    return rebuild(lookup; data=newindex, span=Explicit(newbounds))
end

_without_mapped_crs(f, A) = _without_mapped_crs(f, A, mappedcrs(A))
_without_mapped_crs(f, A::AbstractRaster, ::Nothing) = f(A)
function _without_mapped_crs(f, A::AbstractRaster, mappedcrs)
    A = setmappedcrs(A, nothing)
    x = f(A)
    if x isa AbstractRaster
        x = setmappedcrs(x, mappedcrs)
    end
    return x
end
_without_mapped_crs(f, A::AbstractRasterStack, ::Nothing) = f(A)
function _without_mapped_crs(f, A::AbstractRasterStack, mappedcrs) 
    st1 = map(A -> setmappedcrs(A, nothing), st)
    x = f(st1)
    if x isa AbstractRasterStack
        x = map(A -> setmappedcrs(A, mappedcrs), x)
    end
    return x
end

# These are pretty random default, but seem to work
_default_atol(T::Type{<:Float32}) = 100eps(T)
_default_atol(T::Type{<:Float64}) = 1000eps(T)
_default_atol(T::Type{<:Integer}) = T(1)
_default_atol(::Type) = nothing

