using GeometricalPredicates
using Base: Int64, Float64, NamedTuple
using Printf
using Parameters        # helps setting default parameters in structures
using SpecialFunctions: erfc
using GeophysicalModelGenerator
using ScatteredInterpolation

abstract type trench_slab end
abstract type AbstractThermalStructure end


# Basic folder where to contain the function of trench
# Trench structure: structure that contains all the information concerning the trench boundary 
@with_kw_noshow mutable struct Trench <: trench_slab
    n_seg_xy::Int64 =  1     # Number of segment of the trench plane view (for now, 1 segment)
    A::Array{Float64}       =  (0.0,0.0)  # Coordinate 1 {set of coordinates}
    B::Array{Float64}     =  (0.0,1.0)  # Coordinate 2 {set of coordinates}
    theta_max::Float64 = 45 # max bending angle, (must be converted into radians)
    type_bending::Symbol = :Ribe     # Mode Ribe | Linear | Costumize 
    n_seg::Int64       = 50         # definition segment 
    L0:: Float64       = 400        # length of the slab
    D0:: Float64      = 100        # thickness of the slab 
    WZ:: Float64       = 50        # thickness of the weak zone 
    Lb:: Float64    = 200       # Length at which all the bending is happening (Lb<=L0)
    d_decoupling:: Float64 = 100       # decoupling depth of the slab
end


@with_kw_noshow mutable struct McKenzie_subducting_slab <: AbstractThermalStructure
    Tsurface::Float64 = 20.0       # top T
    Tmantle::Float64  = 1350.0     # bottom T
    Age::Float64     = 60.0          # thermal age of plate [in Myrs]
    Adiabat::Float64  = 0.4        # Adiabatic gradient in K/km
    v_s::Float64      = 2.0      # velocity of subduction  [cm/yrs]
    Cp::Float64       = 1050.0     # Heat capacity   []
    k::Float64        = 3.0        # Heat conductivity 
    rho::Float64      = 3300.0     # denisty of the mantle [km/m3]
    it::Int64       = 36       # number of harmonic summation (look Mckenzie formula)
end
"""
    Compute_ThermalStructureSlab(Temp, X, l, d, s::McKenzie_subducting_slab,ldc,t::Trench)
Compute the thermal structure for the slab: first compute the half-space cooling model. Then compute the McKenzie solution
then as a function of l (length of the slab) it does a weigthed average between the two field. At ldc {length of decoupling} the temperature is dominated by the Mckenzie solution
"""
function Compute_ThermalStructureSlab(Temp, Z, l, d, s::McKenzie_subducting_slab,ldc,t::Trench)

    @unpack Tsurface, Tmantle, Adiabat, Age, v_s, Cp, k, rho, it = s
    @unpack L0, D0, d_decoupling = t

    # calculate halfspace cooling temperature
    kappa       =   1e-6
    SecYear     =   3600*24*365
    ThermalAge  =   Age*1e6*SecYear

    # Compute the + Adiabat*abs.(Z)
    for i in eachindex(Temp)
        Temp[i] = (Tmantle) +  (Tsurface .-Tmantle)*erfc((abs.(d[i])*1e3)./(2*sqrt(kappa*ThermalAge)))
    end

    # Convert the velocity 
    v_s = v_s/(100*365.25*60*60*24);
    
    # calculate the Reynolds number
    Re = (rho*Cp*v_s*D0*1000)/2/k;

    # McKenzie model
    sc = 1/D0
    σ  = zeros(size(Temp));
    Temp_McK = zeros(size(Temp));
    # Dividi et impera
    for i=1:it
        a   = (-1).^(i)./(i.*pi)
        b   = (Re .- (Re.^2 .+ i .^2. *pi .^2) .^(0.5)) .*l .*sc;
        c   = sin.(i .*pi .*(1 .- abs.(d .*sc))) ;
        e   = exp.(b);
        σ .+= a.*e.*c
    end

    Temp_McK           .= (Tmantle) .+2 .*(Tmantle-Tsurface).*σ;
    ind_neg             = findall(Temp_McK .< Tsurface);
    Temp_McK[ind_neg]  .= Tsurface; 


    weight = 0.1 .+(0.9-0.1) ./(ldc*sc) .*(l*sc)
    ind_1 = findall(weight .>1.0);
    ind_2 =findall(weight .<0.1);
    weight[ind_1] .= 1.0; 
    weight[ind_2] .= 0.1;


    Temp .= weight .*Temp_McK+(1 .-weight) .*Temp;
    Temp .= Temp + Adiabat*abs.(Z);

    return Temp
end

"""
    compute_slab_surface!(D0::Float64,L0::Float64,Lb::Float64,WZ::Float64,n_seg::Int64,theta_max::Float64,type_bending::Symbol)

Compute the coordinate of the slab top, bottom surface using the mid surface of the slab as reference. It computes it by discretizing the slab surface
in n segments, and computing the average bending angle {which is a function of the current length of the slab}. Then compute the coordinate 
assuming that the trench is at 0.0, and assuming a positive theta_max angle. 

"""
function compute_slab_surface!(D0::Float64,L0::Float64,Lb::Float64,WZ::Float64,n_seg::Int64,theta_max::Float64,type_bending::Symbol)
    # Convert theta_max into radians
    theta_max = theta_max*pi/180;

    # Allocate the top,mid and bottom surface, and the weakzone 
    Top = zeros(n_seg+1,2);

    Bottom = zeros(n_seg+1,2);
    Bottom[1,2] = -D0; 

    MidS    = zeros(n_seg+1,2);
    MidS[1,2] = -D0./2;

    WZ_surf     = zeros(n_seg+1,2);
    WZ_surf[1,2] = WZ;
    # Initialize the length. 
    l = 0.0;   # initial length 

    it = 1;  # iteration 

    dl = L0/n_seg; # dl 

    while l<L0
        ln = l+dl
        # Compute the mean angle within the segment
        theta_mean = (compute_bending_angle!(theta_max,Lb,l,type_bending)+compute_bending_angle!(theta_max,Lb,ln,type_bending))./2;
        # Compute the mid surface coordinate
        MidS[it+1,1] = MidS[it,1]+dl*cos(theta_mean);

        MidS[it+1,2] = MidS[it,2]-dl.*sin(theta_mean);
        #Compute the top surface coordinate

        Top[it+1,1] = MidS[it+1,1]+0.5.*D0.*abs(sin(theta_mean));

        Top[it+1,2] = MidS[it+1,2]+0.5.*D0.*abs(cos(theta_mean));
        #Compute the bottom surface coordinate

        Bottom[it+1,1] = MidS[it+1,1]-0.5.*D0.*abs(sin(theta_mean));

        Bottom[it+1,2] = MidS[it+1,2]-0.5.*D0.*abs(cos(theta_mean));
        # Compute the top surface for the weak zone 

        WZ_surf[it+1,1] = MidS[it+1,1]+(0.5.*D0+WZ).*abs(sin(theta_mean));

        WZ_surf[it+1,2] = MidS[it+1,2]+(0.5.*D0+WZ).*abs(cos(theta_mean));
        # update l and it
        l = ln;
        it = it+1;
    end

    return Top,Bottom,WZ_surf; #{Filling the structure?}

    end

"""
    compute_bending_angle!(θ_max,Lb,l,type)
θ_max = maximum bending angle 
Lb    = length at which the function of bending is applied (Lb<=L0)
l     = current position within the slab
type  = type of bending {Ribe}{Linear}{Customize} 
function that computes the θ(l). 
"""
function compute_bending_angle!(theta_max::Float64,Lb::Float64,l::Float64,type::Symbol)
    if l>Lb
        return theta_max
    elseif type === :Ribe    
        # Compute theta
        return theta_max*l^2*((3*Lb-2*l))/(Lb^3);
    elseif type === :Linear
        # Compute the actual angle
        return l*(theta_max-0)/(Lb);
    end
end

"""
       transform_coordinate!(X,Y,Z,XT,YT,A,B,direction)
Transform the coordinate such that the new x axis (XT) is parallel to the segment A-B of the slab. The rotation is
anticlockwise. If θ_max is negative, it multiplies YT with the sign of the angle, changing the dip of the subduction.
It returns Bn -> which is the point B coordinate in the new transformed system. 
"""

function transform_coordinate!(X,Y,Z,XT,YT,A,B,direction)
    
    # find the slope of AB points
    s = (B[2]-A[2])/(B[1]-A[1])

    angle_rot = -atand(s); 

    # Shift the origin 
    XT .= X .-A[1]; 

    YT .= Y .-A[2]; 

    Bn = zeros(3);
  
    Bn[1] = B[1]-A[1];
  
    Bn[2] = B[2]-A[2]; 
  
    Bn[3] = 0.0;

    # Transform the coordinates
    Rot3D!(XT,YT,Z, angle_rot, 0);

    YT .= YT*direction; 


    #Find Point B in the new coordinate system 
    roty = [cosd(-0) 0 sind(-0) ; 0 1 0 ; -sind(-0) 0  cosd(-0)];
   
    rotz = [cosd(angle_rot) -sind(angle_rot) 0 ; sind(angle_rot) cosd(angle_rot) 0 ; 0 0 1]

    Bn = rotz* Bn;
    
    Bn = roty*Bn; 

    return Bn

end

"""

"""
function find_slab!(X,Y,Z,d,ls,theta_max,A,B,Top,Bottom,seg_slab,D0,L0)

    # Create the XT,YT 
    XT = zeros(size(X));

    YT = zeros(size(Y)); 

    # Function to transform the coordinate 
    @show sign(theta_max)
    xb = transform_coordinate!(X,Y,Z,XT,YT,A,B,sign(theta_max)); 

    # dl 
    dl = L0/seg_slab; 

    l = 0  # length at the trench position

    # Construct the slab 
    for i = 1:(seg_slab-1)
        
        ln = l+dl; 
        
        pa = (Top[i,1],Top[i,2]); # D = 0 | L = l 

        pb = (Bottom[i,1],Bottom[i,2]); # D = -D0 | L=l 

        pc = (Bottom[i+1,1],Bottom[i+1,2]); # D = -D0 |L=L+dl

        pd = (Top[i+1,1],Top[i+1,2]) # D = 0| L = L+dl

        # Create the polygon 
        poly_y = [pa[1],pb[1],pc[1],pd[1]];

        poly_z = [pa[2],pb[2],pc[2],pd[2]];

        # find a sub set of particles
        ymin = minimum(poly_y);
        
        ymax = maximum(poly_y);
        
        zmin = minimum(poly_z);
        
        zmax = maximum(poly_z);

        ind_s = findall(0.0.<= XT.<= xb[1] .&& ymin .<= YT .<= ymax .&& zmin .<= Z .<= zmax);

        # Find the particles 
        yp = YT[ind_s];
        
        zp = Z[ind_s];
        
        # Initialize the ind that are going to be used by inpoly
        ind = zeros(Bool,size(zp));
        
        # inPoly! [Written by Arne Spang, must be updated with the new version]
        inPolygon!(ind,poly_y,poly_z,yp,zp); 
        
        # indexes of the segment
        ind_seg = ind_s[ind]
        
        # Prepare the variable to interpolate {I put here because to allow also a variation of thickness of the slab}
        D = [0.0,-D0,-D0,0.0];

        L = [l,l,ln,ln];
        
        # Interpolations
        points = [pa[1] pa[2];pb[1] pb[2];pc[1] pc[2];pd[1] pd[2]]'
        
        itp1 = interpolate(Shepard(), points, D);
        
        itp2 = interpolate(Shepard(), points, L);
        
        # Loop over the chosen particles and interpolate the current value of L and D. 
        particle_n = length(ind_seg)

        for ip = 1:particle_n

            point_ = [YT[ind_seg[ip]],Z[ind_seg[ip]]];

            d[ind_seg[ip]] = evaluate(itp1,point_)[1];

            ls[ind_seg[ip]] = evaluate(itp2,point_)[1];
        end

        #Update l
        l = ln; 
    end
end


"""
    create_slab!(X::Array{Float64},Y::Array{Float64},Z::Array{Float64},Ph::Array{Int32},T::Array{Float64},t::Trench,strat,temp)

Main function that creates the slab. Unpack the variable from the structure, create two arrays that contains the information of l and d. 
    l and d are the array containing the length of the slab per each coordinate belonging to the slab, and the distance from the surface. 
    In this function compute_slab_surface and find_slab are called. And after d and ls are computed, it fill up the temperature and phase arrays


"""
function create_slab!(X::Array{Float64},Y::Array{Float64},Z::Array{Float64},Ph::Array{Int32},T::Array{Float64},t::Trench,strat,temp)

    d = ones(size(X)).*NaN64;

    # -> l = length from the trench along the slab 
    ls = ones(size(X)).*NaN64;


    D0 = t.D0; 

    L0 = t.L0;

    n_seg = t.n_seg;

    Lb    = t.Lb; 

    theta_max = t.theta_max;

    n_seg_xy = t.n_seg_xy; 

    WZ = t.WZ; 

    # Allocate d-l array and A,B

    A = zeros(Float64,2,1);
    B = zeros(Float64,2,1);

    #Loop over the segment of the slab
    # In theory this loop would loop all the segment of the trench, but if I introduce a n_seg_xy = 1, it throws me an error concerning the iteration. I tried to fix, but the the error message is mysterious 

    #for is =1:1
        is = 1;
        A[1] = t.A[is];
        A[2] = t.A[is+1];
        B[1] = t.B[is];
        B[2] = t.B[is+1];


        # Compute Top-Bottom surface 
        # Or loop over the segment of the top/bottom surface and inpolygon each element or 

        Top,Bottom,WZ_surf =compute_slab_surface!(D0,L0,Lb,WZ,n_seg,abs(theta_max),t.type_bending);

        find_slab!(X,Y,Z,d,ls,t.theta_max,A,B,Top,Bottom,t.n_seg,t.D0,t.L0);
        
        l_decouplingind = findall(Top[:,2].<=-t.d_decoupling);

        l_decoupling = Top[l_decouplingind[1],1];
        

        # Function to fill up the temperature and the phase. I need to personalize addbox! 

        ind = findall((-D0 .<= d .<= 0.0));

        # Compute thermal structure accordingly. See routines below for different options {Future: introducing the length along the trench for having lateral varying properties along the trench}
        T[ind] = Compute_ThermalStructureSlab(T[ind], Z[ind], ls[ind], d[ind],temp,l_decoupling,t);

        # Set the phase. Different routines are available for that - see below.
        Ph[ind] = Compute_Phase(Ph[ind], T[ind], XT[ind], ls[ind], d[ind], strat)

        # Place holder of the weak zone: it is simply using the find slab routine, and cutting it at d_decoupling. 

    #end

end
"""
I was not able to export this function. When I port this function on the Setup_geometry, this part is completely useless. 
"""

# Internal function that rotates the coordinates
function Rot3D!(X,Y,Z, StrikeAngle, DipAngle)

    # rotation matrixes
    roty = [cosd(-DipAngle) 0 sind(-DipAngle) ; 0 1 0 ; -sind(-DipAngle) 0  cosd(-DipAngle)];
    rotz = [cosd(StrikeAngle) -sind(StrikeAngle) 0 ; sind(StrikeAngle) cosd(StrikeAngle) 0 ; 0 0 1]

    for i in eachindex(X)
        CoordVec = [X[i], Y[i], Z[i]]
        CoordRot =  rotz*CoordVec;
        CoordRot =  roty*CoordRot;
        X[i] = CoordRot[1];
        Y[i] = CoordRot[2];
        Z[i] = CoordRot[3];
    end

    return nothing
end


"""
    inPolygon!(INSIDE::Matrix, PolyX::Vector, PolyY::Vector, X::Matrix, Y::Matrix; fast=false)

Checks if points given by matrices `X` and `Y` are in or on (both cases return true) a polygon given by `PolyX` and `PolyY`. Boolean `fast` will trigger faster version that may miss points that are exactly on the edge of the polygon. Speedup is a factor of 3.

"""
function inPolygon!(INSIDE::Matrix{Bool}, PolyX::Vector{T}, PolyY::Vector{T}, X::Matrix{T}, Y::Matrix{T}; fast=false) where T <: Real
    iSteps = collect(eachindex(PolyX))
    jSteps = [length(PolyX); collect(1:length(PolyX)-1)]

    if fast
        for j = 1 : size(X, 2)
            for i = 1 : size(X, 1)
                INSIDE[i,j] = inPolyPointF(PolyX, PolyY, X[i,j], Y[i,j], iSteps, jSteps)
            end
        end
    else
        for j = 1 : size(X, 2)
            for i = 1 : size(X, 1)
                INSIDE[i,j] = (inPolyPoint(PolyX, PolyY, X[i,j], Y[i,j], iSteps, jSteps) || inPolyPoint(PolyY, PolyX, Y[i,j], X[i,j], iSteps, jSteps))
            end
        end
    end
end

"""
    inPolygon!(inside::Vector, PolyX::Vector, PolyY::Vector, x::Vector, y::Vector; fast=false)

Same as above but `inside`, `X` and `Y` and are vectors.

"""
function inPolygon!(inside::Vector{Bool}, PolyX::Vector{T}, PolyY::Vector{T}, x::Vector{T}, y::Vector{T}; fast=false) where T <: Real
    iSteps = collect(eachindex(PolyX))
    jSteps = [length(PolyX); collect(1:length(PolyX)-1)]

    if fast
        for i = eachindex(x)
            inside[i] = inPolyPointF(PolyX, PolyY, x[i], y[i], iSteps, jSteps)
        end
    else
        for i = eachindex(x)
            inside[i] = (inPolyPoint(PolyX, PolyY, x[i], y[i], iSteps, jSteps) || inPolyPoint(PolyY, PolyX, y[i], x[i], iSteps, jSteps))
        end
    end
end

"""
    inPolyPoint(PolyX::Vector, PolyY::Vector, x::Number, y::Number, iSteps::Vector, jSteps::)

Checks if a point given by x and y is in or on (both cases return true) a polygon given by PolyX and PolyY, iSteps and jSteps provide the connectivity between the polygon edges. This function should be used through inPolygon!().

"""
function inPolyPoint(PolyX::Vector{T}, PolyY::Vector{T}, x::T, y::T, iSteps::Vector{Int64}, jSteps::Vector{Int64}) where T <: Real
    inside1, inside2, inside3, inside4 = false, false, false, false
    for (i,j) in zip(iSteps, jSteps)
        xi = PolyX[i]
        yi = PolyY[i]
        xj = PolyX[j]
        yj = PolyY[j]

        con1 = ((yi > y) != (yj > y))
        con2 = ((yi >= y) != (yj >= y))
        if con1 && (x > (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)
            inside1 = !inside1
        end

        if con1 && (x >= (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)
            inside2 = !inside2
        end

        if con2 && (x > (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)
            inside3 = !inside3
        end

        if con2 && (x >= (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)
            inside4 = !inside4
        end
    end
    return ((inside1 || inside2) || (inside3 || inside4))
end

"""
    inPolyPointF(PolyX::Vector, PolyY::Vector, x::Number, y::Number, iSteps::Vector, jSteps::)

Faster version of inPolyPoint() but will miss some points that are on the edge of the polygon.

"""
function inPolyPointF(PolyX::Vector{T}, PolyY::Vector{T}, x::T, y::T, iSteps::Vector{Int64}, jSteps::Vector{Int64}) where T <: Real
    inside = false
    for (i,j) in zip(iSteps, jSteps)
        xi = PolyX[i]
        yi = PolyY[i]
        xj = PolyX[j]
        yj = PolyY[j]

        if ((yi > y) != (yj > y)) && (x > (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)
            inside = !inside
        end
    end
    return inside
end