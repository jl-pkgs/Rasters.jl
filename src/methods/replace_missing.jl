"""
    replace_missing(a::AbstractRaster, newmissingval)
    replace_missing(a::AbstractRasterStack, newmissingval)

Replace missing values in the array or stack with a new missing value,
also updating the `missingval` field/s.

# Example

```jldoctest
using Rasters
A = Raster(WorldClim{Climate}, :prec; month=1) |> replace_missing
missingval(A)
# output
missing
```

"""
replace_missing(x; missingval=missing) = replace_missing(x, missingval)
function replace_missing(A::AbstractRaster{T}, missingval::MV=missing) where {T,MV}
    missingval = convert(promote_type(T, MV), missingval)
    newdata = if ismissing(Rasters.missingval(A))
        if ismissing(missingval)
            copy(parent(read(A)))
        else
            collect(Missings.replace(parent(A), missingval))
        end
    else
        replace(parent(A), Rasters.missingval(A) => missingval)
    end
    return rebuild(A; data=newdata, missingval=missingval)
end
replace_missing(x::RasterSeriesOrStack, args...) = map(A -> replace_missing(A, args...), x)

