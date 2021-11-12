"""
    points(A::AbstractRaster; dims=(YDim, XDim), ignore_missing) => Array{Tuple}

Returns a generator of the points in `A` for dimensions in `dims`,
where points are a tuple of the values in each specified dimension
index.

# Keywords

- `dims` the dimensions to return points from. The first slice of other
    layers will be used.
- `ignore_missing`: wether to ignore missing values in the array when considering
    points. If `true`, all points in the dimensions will be returned, if `false`
    only the points that are not `=== missingval(A)` will be returned.

The order of `dims` determines the order of the points.

$EXPERIMENTAL
"""
function points(A::AbstractRaster; ignore_missing=false, order=(XDim, YDim, ZDim))
    ignore_missing ? _points(A; order) : _points_missing(A; order)
end
function points(dims::DimTuple; order=(XDim, YDim, ZDim))
    indices = DimIndices(dims)
    ordered_dims = DD.dims(dims, order)
    # Lazily reorder the pionts and index into the dims in the generator
    ordered_point(I) = map(ordered_dims, DD.dims(I, ordered_dims)) do d, i
        d[val(i)] 
    end
    return (ordered_point(I) for I in indices)
end

_points(A::AbstractRaster; kw...) = points(dims(A); kw...)
function _points_missing(A::AbstractRaster; order)
    indices = DimIndices(A)
    ordered_dims = dims(A, order)
    # Lazily reorder the points and index into the dims in the generator
    # or return missing if the matching array value is missing
    function ordered_point_or_missing(I) 
        if A[I...] === missingval(A)
            missing
        else
            map((d, i) -> d[val(i)], ordered_dims, DD.dims(I, ordered_dims))
        end
    end
    return (ordered_point_or_missing(I) for I in indices)
end

