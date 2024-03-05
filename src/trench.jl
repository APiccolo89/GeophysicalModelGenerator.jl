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
    n_seg_xy =  1     # Number of segment of the trench plane view (for now, 1 segment)
    A        =  (0.0,0.0)  # Coordinate 1 {set of coordinates}
    B        =  (0.0,1.0)  # Coordinate 2 {set of coordinates}
    # Place holder for segment transformation(i.e., for a given segment you can define a function such that curve and bend the coordinate)
    # Place holder: changing properties along strike: -function to compute the length of the trench
    # function that compute how much a certain property is changing (i.e., amount of subducted sediment/continent/temperature)
    # function that encode how much change a certain properties along the strike) 
    #(in the future should be possible to change the theta max along strike, or the velocity of convergence or the lenght)
    # = Slab portion 
    theta_max = 45 # max bending angle, (must be converted into radians)
    type_bending  = "Ribe"     # Mode Ribe | Linear | Costumize 
    n_seg         = 50         # definition segment 
    L0            = 400        # length of the slab
    D0            = 100        # thickness of the slab 
    WZ            = 50        # thickness of the weak zone 
    Lb            = 200       # Length at which all the bending is happening (Lb<=L0)
    d_decoupling  = 100       # decoupling depth of the slab
end



@with_kw_noshow mutable struct McKenzie_subducting_slab <: AbstractThermalStructure
    Tsurface = 20       # top T
    Tmantle  = 1350     # bottom T
    Age      = 60          # thermal age of plate [in Myrs]
    Adiabat  = 0        # Adiabatic gradient in K/km
    v_s      = 2.0      # velocity of subduction  [cm/yrs]
    Cp       = 1050     # Heat capacity   []
    k        = 3        # Heat conductivity 
    rho      = 3300     # denisty of the mantle [km/m3]
    it       = 36       # number of harmonic summation (look Mckenzie formula)
end

function Compute_ThermalStructureSlab(Temp, X, l, d, s::McKenzie_subducting_slab,ldc,t::Trench)

    @unpack Tsurface, Tmantle, Adiabat, Age, v_s, Cp, k, rho, it = s
    @unpack L0, D0, d_decoupling = t

    # calculate halfspace cooling temperature
    kappa       =   1e-6
    SecYear     =   3600*24*365
    ThermalAge  =   Age*1e6*SecYear

    # Compute the 
    for i in eachindex(Temp)
        Temp[i] = Tmantle+  (Tsurface .-Tmantle)*erfc((abs.(d[i])*1e3)./(2*sqrt(kappa*ThermalAge)))
    end

    @show it ;

    # Convert the velocity 
    v_s = v_s/(100*365.25*60*60*24);
    @show v_s; 
    # calculate the Reynolds number
    Re = (rho*Cp*v_s*D0*1000)/2/k;
    @show Re
    @show k 
    @show Cp


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


    weight = 0.1 .+(0.8-0.1) ./(ldc*sc) .*(l*sc)
    ind_1 = findall(weight .>1.0);
    ind_2 =findall(weight .<0.1);
    weight[ind_1] .= 1.0; 
    weight[ind_2] .= 0.1;

    #@show weight


    Temp .= weight .*Temp_McK+(1 .-weight) .*Temp;

    return Temp
end

# function to compute top surface bottom surface 
# What do I have in mind? 
# => For a given segment of any orientation: 
# 1. Transform the coordinate system such that xT // to the current segment
# 2. On top of that apply any additional transformation (i.e., circular boundary or whatever)
#   e.g., a.] Given A-B rotate the coordinate such that x^T// the segment. 
#         b.] If the A-B segment are the tip of a more complex curve transform the transform the coordinate (yTT=yT-xT^2=0.0)
#         c.] Then select all the particles that belongs to each of the segment of the slab. 

function compute_slab_surface!(D0::Float64,L0::Float64,Lb::Float64,WZ::Float64,n_seg::Int64,theta_max::Float64,type_bending::String)
    # Convert theta_max into radians

    theta_max = theta_max*pi/180;


    # Allocate the top,mid and bottom surface, and the weakzone 
    Top = zeros(n_seg+1,2);

    theta_mean = zeros(n_seg,1); 

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
        theta_mean[it] = (compute_bending_angle!(theta_max,Lb,l,type_bending)+compute_bending_angle!(theta_max,Lb,ln,type_bending))./2;
        # Compute the mid surface coordinate
        MidS[it+1,1] = MidS[it,1]+dl*cos(theta_mean[it]);

        MidS[it+1,2] = MidS[it,2]-dl.*sin(theta_mean[it]);
        #Compute the top surface coordinate

        Top[it+1,1] = MidS[it+1,1]+0.5.*D0.*abs(sin(theta_mean[it]));

        Top[it+1,2] = MidS[it+1,2]+0.5.*D0.*abs(cos(theta_mean[it]));
        #Compute the bottom surface coordinate

        Bottom[it+1,1] = MidS[it+1,1]-0.5.*D0.*abs(sin(theta_mean[it]));

        Bottom[it+1,2] = MidS[it+1,2]-0.5.*D0.*abs(cos(theta_mean[it]));
        # Compute the top surface for the weak zone 

        WZ_surf[it+1,1] = MidS[it+1,1]+(0.5.*D0+WZ).*abs(sin(theta_mean[it]));

        WZ_surf[it+1,2] = MidS[it+1,2]+(0.5.*D0+WZ).*abs(cos(theta_mean[it]));
        # update l and it
        l = ln;
        it = it+1;

    end

    return Top,MidS,Bottom,WZ_surf,theta_mean; #{Filling the structure?}

    end

function compute_bending_angle!(theta_max::Float64,Lb::Float64,l::Float64,type::String)
    # Input argument: 
    # theta_max -> maximum bending angle in radians 
    # Lb        -> the lenght at which the bending of the slab become effectively constant 
    # l         -> the actual length 

    if type == "Ribe"    
        # Compute theta
        theta = theta_max*l^2*((3*Lb-2*l))/(Lb^3);

    elseif type == "Linear"

        # Compute the slope assuming that the minumum angle is 0.0 
        s = (theta_max-0)/(L0);

        # Compute the actual angle
        theta= l*s;
    end
    # If l>L0 -> theta = theta_max
    if l>Lb

        theta=theta_max;
    end

    return theta 
end

function create_slab!(X,Y,Z,Ph,T,t,strat,temp)
    # Spell out the trench structure
    # Loop over the segment avaiable in the structure
    # -> transform the coordinate 
    # -> create XT,YT,ZT such that XT//AB segment and YT is perpendicular
    # -> See if there are curved boundaries -> compute dy to compute YTT | trench == 0.0. 
    #1. Spell out the structure

    D0 = t.D0; 

    L0 = t.L0;

    n_seg = t.n_seg;

    Lb    = t.Lb; 

    theta_max = t.theta_max;

    A = t.A;

    B = t.B;

    n_seg_xy = t.n_seg_xy; 

    WZ = t.WZ; 

    # Allocate d-l structure 

    # -> d = distance from the top surface
    d = ones(size(X)).*NaN64;

    # -> l = length from the trench along the slab 
    ls = ones(size(X)).*NaN64;

    #Loop over the segment of the slab
    for is =1:n_seg_xy

        Point_A = t.A[is];

        Point_B = t.B[is];

        # Compute Top-Bottom surface 
        # Or loop over the segment of the top/bottom surface and inpolygon each element or 
        Top,MidS,Bottom,WZ_surf =compute_slab_surface!(D0,L0,Lb,WZ,n_seg,abs(theta_max),t.type_bending);

        XT,d,ls,xb=find_slab!(X,Y,Z,d,ls,t.theta_max,A,B,Top,Bottom,t.n_seg,t.D0,t.L0);

        l_decouplingind = findall(Top[:,2].<=-t_.d_decoupling);

        l_decoupling = Top[l_decouplingind[1],1];
        

        # Function to fill up the temperature and the phase. I need to personalize addbox! 

        ind = findall((-D0 .<= d .<= 0.0));

        # Compute thermal structure accordingly. See routines below for different options
        T[ind] = Compute_ThermalStructureSlab(T[ind], XT[ind], ls[ind], d[ind],temp,l_decoupling,t);

        # Set the phase. Different routines are available for that - see below.
        Ph[ind] = Compute_Phase(Ph[ind], T[ind], XT[ind], ls[ind], d[ind], strat)

    end

    return d, ls
end

function find_slab!(X,Y,Z,d,ls,theta_max,A,B,Top,Bottom,seg_slab,D0,L0)

    # Create the XT,YT 
    XT = zeros(size(X));

    YT = zeros(size(Y)); 

    # Function to transform the coordinate 
    @show theta_max
    XT,YT, xb = transform_coordinate!(X,Y,XT,YT,A,B,sign(theta_max)); 

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
        inPoly!(poly_y,poly_z,yp,zp,ind); 
        
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

    # The required variable {XT the x vector transformed, d, the layout of d, l, the layout of length, xb the limit of the slab along x}
    return XT, d, ls, xb
end



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

function transform_coordinate!(X,Y,XT,YT,A,B,direction)
    
    # find the slope of AB points
    s = (B[2]-A[2])/(B[1]-A[1])

    @show s 

    angle_rot = -atand(s); 

    @show angle_rot

    # Shift the origin 
    XT = X .-A[1]; 

    YT = Y .-A[2]; 

    Bn = zeros(3);
    Bn[1] = B[1]-A[1];
    Bn[2] = B[2]-A[2]; 
    Bn[3] = 0.0;

    # Transform the coordinates
    Rot3D!(XT,YT,Z, angle_rot, 0);
    @show direction
    YT = YT*direction; 

    roty = [cosd(-0) 0 sind(-0) ; 0 1 0 ; -sind(-0) 0  cosd(-0)];
    rotz = [cosd(angle_rot) -sind(angle_rot) 0 ; sind(angle_rot) cosd(angle_rot) 0 ; 0 0 1]

    Bn = rotz* Bn;
    Bn = roty*Bn; 
    @show Bn 

    return XT,YT,Bn

end


function inPolyPoint(PolyX::Vector{T}, PolyY::Vector{T}, x::T, y::T, inside::Bool, inside2::Bool, iSteps::Vector{Int64}, jSteps::Vector{Int64}) where T <: Real
    for (i,j) in zip(iSteps, jSteps)
        xi = PolyX[i]
        yi = PolyY[i]
        xj = PolyX[j]
        yj = PolyY[j]

        if ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)
            inside = !inside
        end

        if ((yi > y) != (yj > y)) && (x <= (xj - xi) * (y - yi) / (yj - yi + eps()) + xi)
            inside2 = !inside2
        end
    end
    return (inside || inside2)
end

function inPoly!(PolyX::Vector{T}, PolyY::Vector{T}, x::Vector{T}, y::Vector{T}, inside::Vector{Bool}) where T <: Real
    yN  = false
    yN2 = false
    iSteps = collect(eachindex(PolyX))
    jSteps = [length(PolyX); collect(1:length(PolyX)-1)]

    for i = eachindex(x)
        inside[i] = inPolyPoint(PolyX, PolyY, x[i], y[i], yN, yN2, iSteps, jSteps)
    end
end

function inPoly!(PolyX::Vector{T}, PolyY::Vector{T}, X::Matrix{T}, Y::Matrix{T}, INSIDE::Matrix{Bool}) where T <: Real
    yN  = false
    yN2 = false
    iSteps = collect(eachindex(PolyX))
    jSteps = [length(PolyX); collect(1:length(PolyX)-1)]

    for j = 1 : size(X, 2)
        for i = 1 : size(X, 1)
            INSIDE[i,j] = inPolyPoint(PolyX, PolyY, X[i,j], Y[i,j], yN, yN2, iSteps, jSteps)
        end
    end
end