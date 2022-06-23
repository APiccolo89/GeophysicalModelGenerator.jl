# This is data_types.jl
# contains type definitions to be used in GeophysicalModelGenerator

import Base: show

export  GeoData, ParaviewData, UTMData, CartData,
        LonLatDepthGrid, XYZGrid, Velocity_SphericalToCartesian!,
        Convert2UTMzone, Convert2CartData, ProjectionPoint

"""
    struct ProjectionPoint
        Lon     :: Float64
        Lat     :: Float64
        EW      :: Float64
        NS      :: Float64
        zone    :: Integer
        isnorth :: Bool
    end

Structure that holds the coordinates of a point that is used to project a data set from Lon/Lat to a Cartesian grid and vice-versa.
"""
struct ProjectionPoint
    Lat     :: Float64 
    Lon     :: Float64
    EW      :: Float64
    NS      :: Float64
    zone    :: Integer
    isnorth :: Bool
end

"""
    ProjectionPoint(; Lat=49.9929, Lon=8.2473)

Defines a projection point used for map projections, by specifying latitude and longitude
"""
function ProjectionPoint(; Lat=49.9929, Lon=8.2473)
    # Default = Mainz (center of universe)
    x_lla = LLA(Lat, Lon, 0.0);    # Lat/Lon/Alt of geodesy package 
    x_utmz = UTMZ(x_lla, wgs84)    # UTMZ of 

    ProjectionPoint(Lat, Lon, x_utmz.x, x_utmz.y, Int64(x_utmz.zone), x_utmz.isnorth)
end

"""
    ProjectionPoint(EW::Float64, NS::Float64, Zone::Int64, isnorth::Bool)

Defines a projection point used for map projections, by specifying UTM coordinates (EW/NS), UTM Zone and whether you are on the northern hemisphere

"""
function ProjectionPoint(EW::Float64, NS::Float64, Zone::Int64, isnorth::Bool)
    
    x_utmz = UTMZ(EW,NS,0.0,Zone, isnorth)    # UTMZ of 
    x_lla = LLA(x_utmz, wgs84);    # Lat/Lon/Alt of geodesy package 
    
    ProjectionPoint(x_lla.lat, x_lla.lon, EW, NS, Zone, isnorth)
end


# data structure for a list of values - TO BE REMOVED
mutable struct ValueList
    name::String
    unit::String
    values::Vector{Float64}
end

""" 
    GeoData(lon::Any, lat:Any, depth::GeoUnit, fields::NamedTuple)
    
Data structure that holds one or several fields with longitude, latitude and depth information.

- `depth` can have units of meter, kilometer or be unitless; it will be converted to km.
- `fields` should ideally be a NamedTuple which allows you to specify the names of each of the fields. 
- In case you only pass one array we will convert it to a NamedTuple with default name.
- A single field should be added as `(DataFieldName=Data,)` (don't forget the comma at the end).
- Multiple fields  can be added as well. `lon`,`lat`,`depth` should all have the same size as each of the `fields`.
- In case you want to display a vector field in paraview, add it as a tuple: `(Velocity=(Veast,Vnorth,Vup), Veast=Veast, Vnorth=Vnorth, Vup=Vup)`; we automatically apply a vector transformation when transforming this to a `ParaviewData` structure from which we generate Paraview output. As this changes the magnitude of the arrows, you will no longer see the `[Veast,Vnorth,Vup]` components in Paraview which is why it is a good ideas to store them as separate Fields.
- Yet, there is one exception: if the name of the 3-component field is `colors`, we do not apply this vector transformation as this field is regarded to contain RGB colors. 
- `Lat`,`Lon`,`Depth` should have the same size as the `Data` array. The ordering of the arrays is important. If they are 3D arrays, as in the example below, we assume that the first dimension corresponds to `lon`, second dimension to `lat` and third dimension to `depth` (which should be in km). See below for an example.

# Example     
```julia-repl
julia> Lat         =   1.0:3:10.0;
julia> Lon         =   11.0:4:20.0;
julia> Depth       =   (-20:5:-10)*km;
julia> Lon3D,Lat3D,Depth3D = LonLatDepthGrid(Lon, Lat, Depth);
julia> Lon3D
3×4×3 Array{Float64, 3}:
[:, :, 1] =
 11.0  11.0  11.0  11.0
 15.0  15.0  15.0  15.0
 19.0  19.0  19.0  19.0

[:, :, 2] =
 11.0  11.0  11.0  11.0
 15.0  15.0  15.0  15.0
 19.0  19.0  19.0  19.0

[:, :, 3] =
 11.0  11.0  11.0  11.0
 15.0  15.0  15.0  15.0
 19.0  19.0  19.0  19.0
julia> Lat3D
 3×4×3 Array{Float64, 3}:
 [:, :, 1] =
  1.0  4.0  7.0  10.0
  1.0  4.0  7.0  10.0
  1.0  4.0  7.0  10.0
 
 [:, :, 2] =
  1.0  4.0  7.0  10.0
  1.0  4.0  7.0  10.0
  1.0  4.0  7.0  10.0
 
 [:, :, 3] =
  1.0  4.0  7.0  10.0
  1.0  4.0  7.0  10.0
  1.0  4.0  7.0  10.0
julia> Depth3D
  3×4×3 Array{Unitful.Quantity{Float64, 𝐋, Unitful.FreeUnits{(km,), 𝐋, nothing}}, 3}:
  [:, :, 1] =
   -20.0 km  -20.0 km  -20.0 km  -20.0 km
   -20.0 km  -20.0 km  -20.0 km  -20.0 km
   -20.0 km  -20.0 km  -20.0 km  -20.0 km
  
  [:, :, 2] =
   -15.0 km  -15.0 km  -15.0 km  -15.0 km
   -15.0 km  -15.0 km  -15.0 km  -15.0 km
   -15.0 km  -15.0 km  -15.0 km  -15.0 km
  
  [:, :, 3] =
   -10.0 km  -10.0 km  -10.0 km  -10.0 km
   -10.0 km  -10.0 km  -10.0 km  -10.0 km
   -10.0 km  -10.0 km  -10.0 km  -10.0 km
julia> Data        =   zeros(size(Lon3D));
julia> Data_set    =   GeophysicalModelGenerator.GeoData(Lon3D,Lat3D,Depth3D,(DataFieldName=Data,))   
GeoData 
  size      : (3, 4, 3)
  lon       ϵ [ 11.0 : 19.0]
  lat       ϵ [ 1.0 : 10.0]
  depth     ϵ [ -20.0 km : -10.0 km]
  fields    : (:DataFieldName,)
  attributes: ["note"]
```
"""
struct GeoData
    lon     ::  GeoUnit
    lat     ::  GeoUnit 
    depth   ::  GeoUnit
    fields  ::  NamedTuple
    atts    ::  Dict
    
    # Ensure that the data is of the correct format
    function GeoData(lon,lat,depth,fields,atts=nothing)
        
        # check depth & convert it to units of km in case no units are given or it has different length units
        if unit.(depth[1])==NoUnits 
            depth = depth*km                # in case depth has no dimensions
        end
        depth = uconvert.(km,depth)         # convert to km
        depth = GeoUnit(depth)              # convert to GeoUnit structure with units of km

        if isa(lat, StepRangeLen)
            lat = Vector(lat);
        end

        if isa(lon, StepRangeLen)
            lon = Vector(lon);
        end
        
        # Check ordering of the arrays in case of 3D
        if sum(size(lon).>1)==3
            if maximum(abs.(diff(lon,dims=2)))>maximum(abs.(diff(lon,dims=1))) || maximum(abs.(diff(lon,dims=3)))>maximum(abs.(diff(lon,dims=1)))
                error("It appears that the lon array has a wrong ordering")
            end
            if maximum(abs.(diff(lat,dims=1)))>maximum(abs.(diff(lat,dims=2))) || maximum(abs.(diff(lat,dims=3)))>maximum(abs.(diff(lat,dims=2)))
                error("It appears that the lat array has a wrong ordering")
            end
        end

        # fields should be a NamedTuple. In case we simply provide an array, lets transfer it accordingly
        if !(typeof(fields)<: NamedTuple)
            if (typeof(fields)<: Tuple)
                if length(fields)==1
                    fields = (DataSet1=first(fields),)  # The field is a tuple; create a NamedTuple from it
                else
                    error("Please employ a NamedTuple as input, rather than  a Tuple")  # out of luck
                end
            else
                fields = (DataSet1=fields,)
            end
        end

        DataField = fields[1];
        if typeof(DataField)<: Tuple
            DataField = DataField[1];           # in case we have velocity vectors as input
        end

        if !(size(lon)==size(lat)==size(depth)==size(DataField))    
            error("The size of Lon/Lat/Depth and the Fields should all be the same!")
        end
        
        if isnothing(atts)
            # if nothing is given as attributes, then we note that in GeoData
            atts = Dict("note" => "No attributes were given to this dataset")
        else
            # check if a dict was given
            if !(typeof(atts)<: Dict)
                error("Attributes should be given as Dict!")
            end
        end 
        
        return new(lon,lat,depth,fields,atts)

     end

end

# Print an overview of the Geodata struct:
function Base.show(io::IO, d::GeoData)
    println(io,"GeoData ")
    println(io,"  size      : $(size(d.lon))")
    println(io,"  lon       ϵ [ $(minimum(d.lon.val)) : $(maximum(d.lon.val))]")
    println(io,"  lat       ϵ [ $(minimum(d.lat.val)) : $(maximum(d.lat.val))]")
    println(io,"  depth     ϵ [ $(minimum(d.depth.val)) : $(maximum(d.depth.val))]")
    println(io,"  fields    : $(keys(d.fields))")
    if any( propertynames(d) .== :atts)
        println(io,"  attributes: $(keys(d.atts))")
    end
end



"""
    ParaviewData(x::GeoUnit, y::GeoUnit, z::GeoUnit, values::NamedTuple)

Cartesian data in `x/y/z` coordinates to be used with Paraview.
This is usually generated automatically from the `GeoData` structure, but you can also invoke do this manually:

```julia-repl
julia> Data_set    =   GeophysicalModelGenerator.GeoData(1.0:10.0,11.0:20.0,(-20:-11)*km,(DataFieldName=(-20:-11),))   
julia> Data_cart = convert(ParaviewData, Data_set)
```
"""
mutable struct ParaviewData
    x       ::  GeoUnit
    y       ::  GeoUnit
    z       ::  GeoUnit
    fields  ::  NamedTuple
end

# Print an overview of the ParaviewData struct:
function Base.show(io::IO, d::ParaviewData)
    println(io,"ParaviewData ")
    println(io,"  size  : $(size(d.x))")
    println(io,"  x     ϵ [ $(minimum(d.x.val)) : $(maximum(d.x.val))]")
    println(io,"  y     ϵ [ $(minimum(d.y.val)) : $(maximum(d.y.val))]")
    println(io,"  z     ϵ [ $(minimum(d.z.val)) : $(maximum(d.z.val))]")
    println(io,"  fields: $(keys(d.fields))")
end

# conversion function from GeoData -> ParaviewData
function Base.convert(::Type{ParaviewData}, d::GeoData)  
    
    # Utilize the Geodesy.jl package & use the Cartesian Earth-Centered-Earth-Fixed (ECEF) coordinate system
    lon         =   Array(ustrip.(d.lon.val));
    lat         =   Array(ustrip.(d.lat.val));
    LLA_Data    =   LLA.(lat,lon, Array(ustrip.(d.depth.val))*1000);            # convert to LLA from Geodesy package
    X,Y,Z       =   zeros(size(lon)), zeros(size(lon)), zeros(size(lon));
    
    # convert to cartesian ECEF reference frame. Note that we use kilometers and the wgs84
    for i in eachindex(X)
        data_xyz = ECEF(LLA_Data[i], wgs84)        
        X[i] = data_xyz.x/1e3;
        Y[i] = data_xyz.y/1e3;
        Z[i] = data_xyz.z/1e3;
    end
    

    # This is the 'old' implementation, which does not employ a reference ellipsoid 
    # X = R .* cosd.( lon ) .* cosd.( lat );
    # Y = R .* sind.( lon ) .* cosd.( lat );
    # Z = R .* sind.( lat );

    # In case any of the fields in the tuple has length 3, it is assumed to be a vector, so transfer it
    field_names = keys(d.fields)
    for i=1:length(d.fields)
        if typeof(d.fields[i]) <: Tuple
            if length(d.fields[i]) == 3
                # the tuple has length 3, which is therefore assumed to be a velocity vector
                
                # If the field name contains the string "color" we do not apply a vector transformation as it is supposed to contain RGB colors
                if !occursin("color", string(field_names[i]))
                    println("Applying a vector transformation to field: $(field_names[i])")
                    Velocity_SphericalToCartesian!(d, d.fields[i])  # Transfer it to x/y/z format
                end
            end
        end
    end



    return ParaviewData(GeoUnit(X),GeoUnit(Y),GeoUnit(Z),d.fields)
end



""" 
    UTMData(EW::Any, NS:Any, depth::GeoUnit, UTMZone::Int, NorthernHemisphere=true, fields::NamedTuple)
    
Data structure that holds one or several fields with UTM coordinates (east-west), (north-south) and depth information.

- `depth` can have units of meters, kilometer or be unitless; it will be converted to meters (as UTMZ is usually in meters)
- `fields` should ideally be a NamedTuple which allows you to specify the names of each of the fields. 
- In case you only pass one array we will convert it to a NamedTuple with default name.
- A single field should be added as `(DataFieldName=Data,)` (don't forget the comma at the end).
- Multiple fields  can be added as well.
- In case you want to display a vector field in paraview, add it as a tuple: `(Velocity=(Veast,Vnorth,Vup), Veast=Veast, Vnorth=Vnorth, Vup=Vup)`; we automatically apply a vector transformation when transforming this to a `ParaviewData` structure from which we generate Paraview output. As this changes the magnitude of the arrows, you will no longer see the `[Veast,Vnorth,Vup]` components in Paraview which is why it is a good ideas to store them as separate Fields.
- Yet, there is one exception: if the name of the 3-component field is `colors`, we do not apply this vector transformation as this field is regarded to contain RGB colors. 
- `Lat`,`Lon`,`Depth` should have the same size as the `Data` array. The ordering of the arrays is important. If they are 3D arrays, as in the example below, we assume that the first dimension corresponds to `lon`, second dimension to `lat` and third dimension to `depth` (which should be in km). See below for an example.

# Example     
```julia-repl
julia> ew          =   422123.0:100:433623.0
julia> ns          =   4.514137e6:100:4.523637e6
julia> depth       =   -5400:250:600
julia> EW,NS,Depth =   XYZGrid(ew, ns, depth);
julia> Data        =   ustrip.(Depth);
julia> Data_set    =   UTMData(EW,NS,Depth,33, true, (FakeData=Data,Data2=Data.+1.))  
UTMData 
  UTM zone : 33-33 North
    size    : (116, 96, 25)
    EW      ϵ [ 422123.0 : 433623.0]
    NS      ϵ [ 4.514137e6 : 4.523637e6]
    depth   ϵ [ -5400.0 m : 600.0 m]
    fields  : (:FakeData, :Data2)
  attributes: ["note"]
```
If you wish, you can convert this from `UTMData` to `GeoData` with
```julia-repl
julia> Data_set1 =  convert(GeoData, Data_set)
GeoData 
  size      : (116, 96, 25)
  lon       ϵ [ 14.075969111533457 : 14.213417764154963]
  lat       ϵ [ 40.77452227533946 : 40.86110443583479]
  depth     ϵ [ -5.4 km : 0.6 km]
  fields    : (:FakeData, :Data2)
  attributes: ["note"]
```
which would allow visualizing this in paraview in the usual manner:
```julia-repl
julia> Write_Paraview(Data_set1, "Data_set1")
1-element Vector{String}:
 "Data_set1.vts"
```
"""
struct UTMData
    EW       ::  GeoUnit
    NS       ::  GeoUnit 
    depth    ::  GeoUnit
    zone     ::  Any
    northern ::  Any
    fields   ::  NamedTuple 
    atts     ::  Dict
    
    # Ensure that the data is of the correct format
    function UTMData(EW,NS,depth,zone,northern,fields,atts=nothing)
        
        # check depth & convert it to units of km in case no units are given or it has different length units
        if unit.(depth)[1]==NoUnits 
            depth = depth*m                # in case depth has no dimensions
        end
        depth = uconvert.(m,depth)         # convert to meters
        depth = GeoUnit(depth)             # convert to GeoUnit structure with units of meters

        # Check ordering of the arrays in case of 3D
        if sum(size(EW).>1)==3
            if maximum(abs.(diff(EW,dims=2)))>maximum(abs.(diff(EW,dims=1))) || maximum(abs.(diff(EW,dims=3)))>maximum(abs.(diff(EW,dims=1)))
                error("It appears that the EW array has a wrong ordering")
            end
            if maximum(abs.(diff(NS,dims=1)))>maximum(abs.(diff(NS,dims=2))) || maximum(abs.(diff(NS,dims=3)))>maximum(abs.(diff(NS,dims=2)))
                error("It appears that the NS array has a wrong ordering")
            end
        end

        # fields should be a NamedTuple. In case we simply provide an array, lets transfer it accordingly
        if !(typeof(fields)<: NamedTuple)
            if (typeof(fields)<: Tuple)
                if length(fields)==1
                    fields = (DataSet1=first(fields),)  # The field is a tuple; create a NamedTuple from it
                else
                    error("Please employ a NamedTuple as input, rather than  a Tuple")  # out of luck
                end
            else
                fields = (DataSet1=fields,)
            end
        end

        DataField = fields[1];
        if typeof(DataField)<: Tuple
            DataField = DataField[1];           # in case we have velocity vectors as input
        end

        if !(size(EW)==size(NS)==size(depth)==size(DataField))    
            error("The size of EW/NS/Depth and the Fields should all be the same!")
        end

        if length(zone)==1
            zone = ones(Int64,size(EW))*zone
            northern = ones(Bool,size(EW))*northern
        end
        
        # take care of attributes
        if isnothing(atts)
            # if nothing is given as attributes, then we note that in GeoData
            atts = Dict("note" => "No attributes were given to this dataset")
        else
            # check if a dict was given
            if !(typeof(atts)<: Dict)
                error("Attributes should be given as Dict!")
            end
        end

        return new(EW,NS,depth,zone,northern, fields,atts)
        
     end

end

# Print an overview of the UTMData struct:
function Base.show(io::IO, d::UTMData)
    println(io,"UTMData ")
    if d.northern[1]
        println(io,"  UTM zone : $(minimum(d.zone))-$(maximum(d.zone)) North")
    else
        println(io,"  UTM zone : $(minimum(d.zone))-$(maximum(d.zone)) South")
    end
    println(io,"    size    : $(size(d.EW))")
    println(io,"    EW      ϵ [ $(minimum(d.EW.val)) : $(maximum(d.EW.val))]")
    println(io,"    NS      ϵ [ $(minimum(d.NS.val)) : $(maximum(d.NS.val))]")
    println(io,"    depth   ϵ [ $(minimum(d.depth.val)) : $(maximum(d.depth.val))]")
    println(io,"    fields  : $(keys(d.fields))")
    if any( propertynames(d) .== :atts)
        println(io,"  attributes: $(keys(d.atts))")
    end
end

"""
Converts a `UTMData` structure to a `GeoData` structure
"""
function Base.convert(::Type{GeoData}, d::UTMData)  

    Lat = zeros(size(d.EW));
    Lon = zeros(size(d.EW));
    for i in eachindex(d.EW.val)

        # Use functions of the Geodesy package to convert to LLA
        utmz_i  = UTMZ(d.EW.val[i],d.NS.val[i],Float64(ustrip.(d.depth.val[i])),d.zone[i],d.northern[i])
        lla_i   = LLA(utmz_i,wgs84)
        lon = lla_i.lon;
       # if lon<0; lon = 360+lon; end # as GMT expects this

        Lat[i] = lla_i.lat
        Lon[i] = lon
    end 

    # handle the case where an old GeoData structure is converted
    if any( propertynames(d) .== :atts)
        atts = d.atts;
    else
        atts = Dict("note" => "No attributes were given to this dataset") # assign the default
    end

    depth = d.depth.val
    if GeophysicalModelGenerator.Unit(d.depth[1])==m
        depth = depth/1000
    end

    return GeoData(Lon,Lat,depth,d.fields,atts)

end

"""
Converts a `GeoData` structure to a `UTMData` structure
"""
function Base.convert(::Type{UTMData}, d::GeoData)  

    EW = zeros(size(d.lon));
    NS = zeros(size(d.lon));
    depth = zeros(size(d.lon));
    zone = zeros(Int64,size(d.lon));
    northern = zeros(Bool,size(d.lon));
    for i in eachindex(d.lon.val)

        # Use functions of the Geodesy package to convert to LLA
        lla_i   = LLA(d.lat.val[i],d.lon.val[i],Float64(ustrip.(d.depth.val[i])*1e3))
        utmz_i  = UTMZ(lla_i, wgs84)
        
        EW[i] = utmz_i.x
        NS[i] = utmz_i.y
        depth[i] = utmz_i.z
        zone[i] = utmz_i.zone;
        northern[i] = utmz_i.isnorth
    end 

    # handle the case where an old GeoData structure is converted
    if any( propertynames(d) .== :atts)
        atts = d.atts;
    else
        atts = Dict("note" => "No attributes were given to this dataset") # assign the default
    end

    return UTMData(EW,NS,depth,zone, northern, d.fields, atts)

end


"""
    Convert2UTMzone(d::GeoData, p::ProjectionPoint)  

Converts a `GeoData` structure to fixed UTM zone, around a given `ProjectionPoint`  
    This useful to use real data as input for a cartesian geodynamic model setup (such as in LaMEM). In that case, we need to project map coordinates to cartesian coordinates.
    One way to do this is by using UTM coordinates. Close to the `ProjectionPoint` the resulting coordinates will be rectilinear and distance in meters. The map distortion becomes larger the further you are away from the center.
      
"""
function Convert2UTMzone(d::GeoData, proj::ProjectionPoint)  

    EW = zeros(size(d.lon));
    NS  = zeros(size(d.lon));
    zone        = zeros(Int64,size(d.lon));
    northern    = zeros(Bool,size(d.lon));
    trans       = UTMfromLLA(proj.zone, proj.isnorth, wgs84) 
    for i in eachindex(d.lon.val)

        # Use functions of the Geodesy package to convert to LLA
        lla_i  =   LLA(d.lat.val[i],d.lon.val[i],Float64(ustrip.(d.depth.val[i])*1e3))
        utm_i  =   trans(lla_i)

        EW[i] = utm_i.x
        NS[i] = utm_i.y
        zone[i] = proj.zone;
        northern[i] = proj.isnorth
    end 

    # handle the case where an old GeoData structure is converted
    if any( propertynames(d) .== :atts)
        atts = d.atts;
    else
        atts = Dict("note" => "No attributes were given to this dataset") # assign the default
    end

    return UTMData(EW,NS,d.depth.val,zone, northern, d.fields,atts)

end



""" 
    CartData(x::Any, y::Any, z::GeoUnit, fields::NamedTuple)
    
Data structure that holds one or several fields with with Cartesian x/y/z coordinates. Distances are in kilometers

- `x`,`y`,`z` can have units of meters, kilometer or be unitless; they will be converted to kilometers
- `fields` should ideally be a NamedTuple which allows you to specify the names of each of the fields. 
- In case you only pass one array we will convert it to a NamedTuple with default name.
- A single field should be added as `(DataFieldName=Data,)` (don't forget the comma at the end).
- Multiple fields  can be added as well.
- In case you want to display a vector field in paraview, add it as a tuple: `(Velocity=(Vx,Vnorth,Vup), Veast=Veast, Vnorth=Vnorth, Vup=Vup)`; we automatically apply a vector transformation when transforming this to a `ParaviewData` structure from which we generate Paraview output. As this changes the magnitude of the arrows, you will no longer see the `[Veast,Vnorth,Vup]` components in Paraview which is why it is a good ideas to store them as separate Fields.
- Yet, there is one exception: if the name of the 3-component field is `colors`, we do not apply this vector transformation as this field is regarded to contain RGB colors. 
- `x`,`y`,`z` should have the same size as the `Data` array. The ordering of the arrays is important. If they are 3D arrays, as in the example below, we assume that the first dimension corresponds to `x`, second dimension to `y` and third dimension to `z` (which should be in km). See below for an example.

# Example     
```julia-repl
julia> x        =   0:2:10
julia> y        =   -5:5
julia> z        =   -10:2:2
julia> X,Y,Z    =   XYZGrid(x, y, z);
julia> Data     =   Z
julia> Data_set =   CartData(X,Y,Z, (FakeData=Data,Data2=Data.+1.))
CartData 
    size    : (6, 11, 7)
    x       ϵ [ 0.0 km : 10.0 km]
    y       ϵ [ -5.0 km : 5.0 km]
    z       ϵ [ -10.0 km : 2.0 km]
    fields  : (:FakeData, :Data2)
  attributes: ["note"]
```
`CartData` is particularly useful in combination with cartesian geodynamic codes, such as LaMEM, which require cartesian grids.
You can directly save your data to Paraview with
```julia-repl
julia> Write_Paraview(Data_set, "Data_set")
1-element Vector{String}:
 "Data_set.vts"
```

If you wish, you can convert this to `UTMData` (which will simply convert the )
```julia-repl
julia> Data_set1 =  convert(GeoData, Data_set)
GeoData 
  size  : (116, 96, 25)
  lon   ϵ [ 14.075969111533457 : 14.213417764154963]
  lat   ϵ [ 40.77452227533946 : 40.86110443583479]
  depth ϵ [ -5.4 km : 0.6 km]
  fields: (:FakeData, :Data2)
```
which would allow visualizing this in paraview in the usual manner:

"""
struct CartData
    x       ::  GeoUnit
    y       ::  GeoUnit 
    z       ::  GeoUnit
    fields  ::  NamedTuple
    atts    ::  Dict 
    
    # Ensure that the data is of the correct format
    function CartData(x,y,z,fields,atts=nothing)
       
        # Check ordering of the arrays in case of 3D
        if sum(size(x).>1)==3
            if maximum(abs.(diff(x,dims=2)))>maximum(abs.(diff(x,dims=1))) || maximum(abs.(diff(x,dims=3)))>maximum(abs.(diff(x,dims=1)))
                error("It appears that the x-array has a wrong ordering")
            end
            if maximum(abs.(diff(y,dims=1)))>maximum(abs.(diff(y,dims=2))) || maximum(abs.(diff(y,dims=3)))>maximum(abs.(diff(y,dims=2)))
                error("It appears that the y-array has a wrong ordering")
            end
        end

        # check depth & convert it to units of km in case no units are given or it has different length units
        x = Convert!(x,km)
        y = Convert!(y,km)
        z = Convert!(z,km)
        
        # fields should be a NamedTuple. In case we simply provide an array, lets transfer it accordingly
        if !(typeof(fields)<: NamedTuple)
            if (typeof(fields)<: Tuple)
                if length(fields)==1
                    fields = (DataSet1=first(fields),)  # The field is a tuple; create a NamedTuple from it
                else
                    error("Please employ a NamedTuple as input, rather than a Tuple")  # out of luck
                end
            else
                fields = (DataSet1=fields,)
            end
        end

        DataField = fields[1];
        if typeof(DataField)<: Tuple
            DataField = DataField[1];           # in case we have velocity vectors as input
        end

        if !(size(x)==size(y)==size(z)==size(DataField))    
            error("The size of x/y/z and the Fields should all be the same!")
        end

        # take care of attributes
        if isnothing(atts)
            # if nothing is given as attributes, then we note that
            atts = Dict("note" => "No attributes were given to this dataset")
        else
            # check if a dict was given
            if !(typeof(atts)<: Dict)
                error("Attributes should be given as Dict!")
            end
        end

        return new(x,y,z,fields,atts)
        
     end

end

# Print an overview of the UTMData struct:
function Base.show(io::IO, d::CartData)
    println(io,"CartData ")
    println(io,"    size    : $(size(d.x))")
    println(io,"    x       ϵ [ $(minimum(d.x.val)) : $(maximum(d.x.val))]")
    println(io,"    y       ϵ [ $(minimum(d.y.val)) : $(maximum(d.y.val))]")
    println(io,"    z       ϵ [ $(minimum(d.z.val)) : $(maximum(d.z.val))]")
    println(io,"    fields  : $(keys(d.fields))")
    if any( propertynames(d) .== :atts)
        println(io,"  attributes: $(keys(d.atts))")
    end
end

"""
    CartData(xyz::Tuple{Array,Array,Array})

This creates a `CartData` struct if you have a Tuple with 3D coordinates as input.
# Example 
```julia
julia> data = CartData(XYZGrid(-10:10,-5:5,0))
CartData 
    size    : (21, 11, 1)
    x       ϵ [ -10.0 km : 10.0 km]
    y       ϵ [ -5.0 km : 5.0 km]
    z       ϵ [ 0.0 km : 0.0 km]
    fields  : (:Z,)
  attributes: ["note"]
```
"""
CartData(xyz::Tuple) = CartData(xyz[1],xyz[2],xyz[3],(Z=xyz[3],))



"""
    Convert2UTMzone(d::CartData, proj::ProjectionPoint)  

This transfers a `CartData` dataset to a `UTMData` dataset, that has a single UTM zone. The point around which we project is `ProjectionPoint`
"""
function Convert2UTMzone(d::CartData, proj::ProjectionPoint)  

    return UTMData(ustrip.(d.x.val).*1e3 .+ proj.EW,ustrip.(d.y.val).*1e3 .+ proj.NS,
                   ustrip.(d.z.val).*1e3,proj.zone, proj.isnorth, d.fields, d.atts)

end

"""
    Convert2CartData(d::UTMData, proj::ProjectionPoint)
Converts a `UTMData` structure to a `CartData` structure, which essentially transfers the dimensions to km
"""
function Convert2CartData(d::UTMData, proj::ProjectionPoint)  

        # handle the case where an old structure is converted
        if any( propertynames(d) .== :atts)
            atts = d.atts;
        else
            atts = Dict("note" => "No attributes were given to this dataset") # assign the default
        end

    return CartData( (ustrip.(d.EW.val) .- proj.EW)./1e3, (ustrip.(d.NS.val) .- proj.NS)./1e3,
                     ustrip.(d.depth.val)./1e3, d.fields,atts)
end


"""
    Convert2CartData(d::GeoData, proj::ProjectionPoint)
Converts a `GeoData` structure to a `CartData` structure, which essentially transfers the dimensions to km
"""
function Convert2CartData(d::GeoData, proj::ProjectionPoint)  

    d_UTM = Convert2UTMzone(d,proj)
    return CartData( (ustrip.(d_UTM.EW.val) .- proj.EW)./1e3, (ustrip.(d_UTM.NS.val) .- proj.NS)./1e3,
                     ustrip.(d_UTM.depth.val)./1e3, d_UTM.fields,d_UTM.atts)
end

"""
    Lon, Lat, Depth = LonLatDepthGrid(Lon::Any, Lat::Any, Depth:Any)

Creates 3D arrays of `Lon`, `Lat`, `Depth` from 1D vectors or numbers

# Example 1: Create 3D grid
```julia-repl
julia> Lon,Lat,Depth =  LonLatDepthGrid(10:20,30:40,(-10:-1)km);
julia> size(Lon)
(11, 11, 10)
```

# Example 2: Create 2D lon/lat grid @ a given depth
```julia-repl
julia> Lon,Lat,Depth =  LonLatDepthGrid(10:20,30:40,-50km);
julia> size(Lon)
(11, 11)
```

# Example 3: Create 2D lon/depth grid @ a given lat
```julia-repl
julia> Lon,Lat,Depth =  LonLatDepthGrid(10:20,30,(-10:-1)km);
julia> size(Lon)
(11, 11)
```
# Example 4: Create 1D vertical line @ a given lon/lat point
```julia-repl
julia> Lon,Lat,Depth =  LonLatDepthGrid(10,30,(-10:-1)km);
julia> size(Lon)
(10, )
```

"""
function LonLatDepthGrid(Lon::Any, Lat::Any, Depth::Any)

    nLon    = length(Lon)
    nLat    = length(Lat)
    nDepth  = length(Depth)

    if nLon==nLat==nDepth==1
        error("Cannot use this routine for a 3D point (no need to create a grid in that case")
    end 
    if maximum([length(size(Lon)), length(size(Lat)), length(size(Depth))])>1
        error("You can only give 1D vectors or numbers as input")
    end

    Lon3D   =   zeros(nLon,nLat,nDepth);
    Lat3D   =   zeros(nLon,nLat,nDepth);
    Depth3D =   zeros(nLon,nLat,nDepth);

    for i=1:nLon
        for j=1:nLat
            for k=1:nDepth
                Lon3D[i,j,k]    =   ustrip.(Lon[i]);
                Lat3D[i,j,k]    =   ustrip.(Lat[j]);
                Depth3D[i,j,k]  =   ustrip.(Depth[k]);
            end
        end
    end

    # Add dimensions back
    Lon3D   =   Lon3D*unit(  Lon[1])
    Lat3D   =   Lat3D*unit(  Lat[1])
    Depth3D = Depth3D*unit(Depth[1])

    return Lon3D, Lat3D, Depth3D
end

"""
    X,Y,Z = XYZGrid(X_vec::Any, Y_vec::Any, Z_vec::Any)

Creates a `X,Y,Z` grid. It works just as `LonLatDepthGrid` apart from the better suited name.

# Example 1: Create 3D grid
```julia-repl
julia> X,Y,Z =  XYZGrid(10:20,30:40,(-10:-1)km);
julia> size(X)
(11, 11, 10)
```

See `LonLatDepthGrid` for more examples.
"""
function XYZGrid(X_vec::Any, Y_vec::Any, Z_vec::Any)
    return X,Y,Z = LonLatDepthGrid(X_vec,Y_vec,Z_vec)
end


"""
    Velocity_SphericalToCartesian!(Data::GeoData, Velocity::Tuple)

In-place conversion of velocities in spherical velocities `[Veast, Vnorth, Vup]` to cartesian coordinates (for use in paraview).

NOTE: the magnitude of the vector will be the same, but the individual `[Veast, Vnorth, Vup]` components
will not be retained correctly (as a different `[x,y,z]` coordinate system is used in paraview). 
Therefore, if you want to display or color that correctly in Paraview, you need to store these magnitudes as separate fields

"""
function Velocity_SphericalToCartesian!(Data::GeoData, Velocity::Tuple)
    # Note: This is partly based on scripts originally written by Tobias Baumann, Uni Mainz 

    for i in eachindex(Data.lat.val)
        az  =   Data.lon.val[i];
        el  =   Data.lat.val[i];

        R           = [-sind(az) -sind(el)*cosd(az) cosd(el)*cosd(az);
                        cosd(az) -sind(el)*sind(az) cosd(el)*sind(az); 
                        0.0       cosd(el)          sind(el)            ];
        
        V_sph       =   [Velocity[1][i]; Velocity[2][i]; Velocity[3][i] ];
       
        # Normalize spherical velocity
        V_mag       =  sum(sqrt.(V_sph.^2));        # magnitude
        V_norm      =  V_sph/V_mag                  

        V_xyz_norm  =  R*V_norm;
        V_xyz       =  V_xyz_norm.*V_mag;          # scale with magnitude

        # in-place saving of rotated velocity    
        Velocity[1][i] = V_xyz[1];  
        Velocity[2][i] = V_xyz[2];
        Velocity[3][i] = V_xyz[3];
    end
end

# Internal function that converts arrays to a GeoUnit with certain units
function Convert!(d,u)
    if unit.(d)[1]==NoUnits 
        d = d*u                # in case it has no dimensions
    end
    d = uconvert.(u,d)         # convert to u
    d = GeoUnit(d)             # convert to GeoUnit structure with units of u

    return d
end