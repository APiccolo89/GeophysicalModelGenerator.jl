using GeometricalPredicates
using Base: Int64, Float64, NamedTuple
using Printf
using Parameters        # helps setting default parameters in structures
using SpecialFunctions: erfc
abstract type trench_slab end


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
    WZ            = 50         # thickness of the weak zone 
    Lb            = 200        # Length at which all the bending is happening (Lb<=L0)
    d_decoupling  = 100        # decoupling depth of the slab
    # -> Something to tell what is the phase of the weak zone
    # -> Something to tell what are the phases of the slab. 
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
    print("$theta_max")

    theta_max = theta_max*pi/180;

    print("$theta_max")

    # Allocate the top,mid and bottom surface, and the weakzone 
    Top = zeros(n_seg+1,2);
    theta_mean = zeros(n_seg,1); 
    Bottom = zeros(n_seg+1,2);
    Bottom[1,2] = -D0; 

    MidS    = zeros(n_seg+1,2);
    MidS[1,2] = -D0./2;

    WZ_surf     = zeros(n_seg+1,2);

    # Initialize the length. 
    l = 0.0;   # initial length 

    it = 1;  # iteration 

    dl = L0/n_seg; # dl 

    while l<L0
        ln = l+dl
        # Compute the mean angle within the segment

        theta_mean[it] = (compute_bending_angle!(theta_max,Lb,l,type_bending)+compute_bending_angle!(theta_max,Lb,ln,type_bending))./2;
        # Compute the mid surface coordinate
        @show theta_mean[it] 
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
        @show l
        l = ln;
        @show l

        it = it+1;

    end

    return Top,MidS,Bottom,WZ_surf,theta_mean #{Filling the structure?}

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

function create_slab!(X,Y,Z,Ph,T,t::Trench)
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

    # Allocate d-l structure 

    # -> d = distance from the top surface
    d = fill!(array(Float64,size(X),NaN))

    # -> l = length from the trench along the slab 
    l = fill!(array(Float64,size(X),NaN))

    #Loop over the segment of the slab
    for is =1:n_seg

        Point_A = t.A[is];

        Point_B = t.B[is];

        # Compute Top-Bottom surface 
        # Or loop over the segment of the top/bottom surface and inpolygon each element or 

        Top,MidS,Bottom,WZ_surf =compute_slab_surface!(D0,L0,Lb,WZ,n_seg,abs(theta_max),t.type_bending)

        # Compute Top-Bottom surface 
        # Or loop over the segment of the top/bottom surface and inpolygon each element or 
        # interpolate top-bottom surface onto particles -> transform ZT and use inbox 
        #-> Top-bottom surface are giving maximum depth+max y and x1-x2. 
        # 

    end

end

function find_slab!(X,Y,Z,d,l,theta_max,A,B,Top,Bottom,seg_slab,D0,L0)

    # Create the XT,YT 
    XT = zeros(size(X));

    YT = zeros(size(Y)); 

    # Function to transform the coordinate 
    XT,YT, xb = transform_coordinate!(X,Y,XT,YT,A,B); 

    dl = L0/seg_slab; 

    l = 0 
    # Construct the slab 
    for i = 1:(n_seg-1)
        ln = l+dl; 
        
        y_lc = [Top[i,1],Top[i+1,1],Bottom[i,1],Bottom[i+1,1]];
        
        z_lc = [Top[i,2],Top[i+1,2],Bottom[i,2],Bottom[i+1,2]];
        
        
        pa = (Top[i,1],Top[i,2]); # D = 0 | L = l 

        pb = (Bottom[i,1],Bottom[i,2]); # D = -D0 | L=l 

        pc = (Bottom[i+1,1],Bottom[i+1,2]); # D = -D0 |L=L+dl

        pd = (Top[i+1,1],Top[i+1,1]) # D = 0| L = L+dl

        # Create the polygon 

        scaling_x = 1/((maximum(x_lc))-minimum(x_lc))


        # Use inpolygon 
        # find all the points that are likely to be within the polygon 
        ind_part = findall(minimum(z_lc) .<= Z .<= maximum(z_lc) .&& minumum(y_lc) .<= YT .<= maximum(y_lc));
        
        for in = 1:length(ind_part)
            y = YT[ind[in]];
            z = Z[ind[in]];

        end


        # Interpolate into the marker the l/d 
        # Procede
        # -
        # find all the particles within this poligon 
        
        l = ln; 
    end
    # Use box function of Boris using as coordinate XT,l,d 




end

function transform_coordinate!(X,Y,XT,YT,A,B)
    
    # find the slope of AB points
    s = (B[2]-A[2])/(B[2]-B[1])

    angle_rot = tand(s); 

    # Shift the origin 
    XT = X .+A[1]; 

    YT = Y .+A[2]; 

    # Transform the coordinates
    XT = XT.*cos(angle_rot)-YT.*sin(angle_rot);

    YT = XT.*sin(angle_rot)+YT.*cos(angle_rot); 

    xb[1] = B[1]*cos(angle_rot)-B[2]*sin(angle_rot);

    xb[2] = B[1]*sin(angle_rot)+B[2]*cos(angle_rot); 

    return XT,YT,xb

end







