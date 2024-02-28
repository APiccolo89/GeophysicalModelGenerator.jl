# Basic folder where to contain the function of trench
# Trench structure: structure that contains all the information concerning the trench boundary 
@with_kw_noshow mutable struct Trench 
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


function _compute_slab_surface!(t:Trench)
    # Spell out all the component of the surface 
    D0        = t.D0;

    L0        = t.L0;

    Lb        = t.Lb; 

    WZ        = t.WZ; 

    n_seg     = t.n_seg; 

    theta_max = t.theta_max;

    # Convert theta_max into radians

    theta_max = theta_max*pi/180;

    if t.type_bending =="Ribe"

        f_theta(x) = _compute_ribe_bending_angle(theta_max,Lb,x)
    elseif t.type_bending =="Linear"

        f_theta(x) = _compute_linear_bending_angle(theta_max,Lb,x)
    end

    # Allocate the top,mid and bottom surface, and the weakzone 
    Top = zeros(n_seg+1,2);

    Bottom = zeros(n_seg+1,2);

    MidS    = zeros(n_seg+1,2);

    WZ_surf     = zeros(n_seg+1,2);

    # Initialize the length. 
    l = 0.0;   # initial length 

    iter = 1;  # iteration 

    dl = L0/n_seg; # dl 

    while l<L0
        ln = l+dl;
        # Compute the mean angle within the segment

        theta_mean[it] = (ftheta(l)+ftheta(ln))./2;
        # Compute the mid surface coordinate

        MidS[it+1,1] = MidS[it,1]+dl*cos(theta_mean(it));

        MidS[it+1,2] = MidS[it,2]-dl.*sin(theta_mean(it));
        #Compute the top surface coordinate

        Top[it+1,1] = MidS[it+1,1]+0.5.*D0.*abs(sin(theta_mean(it)));

        Top[it+1,2] = MidS[it+1,2]+0.5.*D0.*abs(cos(theta_mean(it)));
        #Compute the bottom surface coordinate

        Bottom[it+1,1] = MidS[it+1,1]-0.5.*D0.*abs(sin(theta_mean(it)));

        Bottom[it+1,2] = MidS[it+1,2]-0.5.*D0.*abs(cos(theta_mean(it)));
        # Compute the top surface for the weak zone 

        WS_surf[it+1,1] = MidS[it+1,1]+(0.5.*D0+WZ_tk).*abs(sin(theta_mean(it)));

        WS_surf[it+1,2] = MidS[it+1,2]+(0.5.*D0+WZ_tk).*abs(cos(theta_mean(it)));
        # update l and it

        l = ln;

        it = it+1;

    end

    return Top,MidS,Bottom,WZ_surf
    end



function _compute_ribe_bending_angle(theta_max::Float64,Lb::Float64,l::float64)
    # Input argument: 
    # theta_max -> maximum bending angle in radians 
    # Lb        -> the lenght at which the bending of the slab become effectively constant 
    # l         -> the actual length 

    
    # Compute theta
    theta = theta*l^2*((3*Lb-2*l))/(Lb^3);

    if l>Lb 
        theta = theta_max; 
    end

    return theta
end

function _compute_linear_bending_angle(theta_max::Float64,Lb::Float64,l::float64)

    # Compute the slope assuming that the minumum angle is 0.0 
    s = (theta_max-0)/(L0);

    # Compute the actual angle
    theta= l*s;

    # If l>L0 -> theta = theta_max
    if l>Lb
        theta=theta;
    end

    return theta 
end












