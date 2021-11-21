# File extensions. GDAL is the catch-all for everything else
const EXT = Dict(
    GRDraster => (".grd", ".gri"), 
    NCDraster => (".nc",), 
    SMAPraster => (".h5",),
    HDF5raster => (".h5",),
)
const REV_EXT = Dict(
    ".grd" => GRDraster, 
    ".gri" => GRDraster, 
    ".nc" => NCDraster, 
    ".h5" => HDF5raster,
)

# Get the source backend for a file extension, falling back to GDALraster
_sourcetype(filenames::NamedTuple) = _sourcetype(first(filenames))
function _sourcetype(filename::AbstractString)
    source = get(REV_EXT, splitext(filename)[2], GDALraster)
    if source === HDF5raster
        source = h5open(filename) do ds
            haskey(ds, SMAPGEODATA) ? SMAPraster : NCDraster
        end
    end
    return source
end

# Internal read method
function _open(f, filename::AbstractString; source=_sourcetype(filename), kw...)
    _open(f, source, filename; kw...)
end

function create(filename::AbstractString, T::Type, dims::Tuple; 
    source=_sourcetype(filename), kw...
)
    create(filename, source, T, dims; kw...)
end
