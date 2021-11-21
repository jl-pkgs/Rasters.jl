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
    MT = promote_type(T, MV)
    mv = convert(MT, missingval)
    repmissing(x) = isequal(x, Rasters.missingval(A)) ? mv : x
    # Disk-backed arrays need to be lazy, memory-backed don't.
    # But in both cases we make sure we return an array with the missingval
    # in the eltype, even if there are no missing values in the array.
    if isdisk(A)
        data = repmissing.(parent(A))
        if missingval isa Missing
            MT = promote_type(eltype(data), MV)
            data = MissingDiskArray(MT, data)
        end
    else
        data = similar(parent(A), MT)
        data .= repmissing.(parent(A))
    end
    return rebuild(A; data, missingval)
end
replace_missing(x::RasterSeriesOrStack, args...) = map(A -> replace_missing(A, args...), x)
