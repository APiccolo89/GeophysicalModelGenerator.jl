# Basic folder where to contain the function of trench
# Trench structure: structure that contains all the information concerning the trench boundary 
@with_kw_noshow mutable struct Trench 
    n_seg_xy =  1     # Number of segment of the trench plane view (for now, 1 segment)
    A        =  (0.0,0.0)  # Coordinate 1 {set of coordinates}
    B        =  (0.0,0.0)  # Coordinate 2 {set of coordinates}
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
end





