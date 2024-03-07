using GeophysicalModelGenerator
include("../src/trench.jl")
using ScatteredInterpolation 
using BenchmarkTools

"""
1. Creation of subduction: 
This example is creating two subductions zone that has different dip angle and common trench. 
The slab are created using the trench structure. 
Trench structure is composed of the following member:
n_seg_xy => the number of segment through which the x-y trench is discretized. Each segment has a point A and B that defines it
Point A : the first point of each segment array Point A = [[coord 1x, coord1y ]....]
Point B : the second point of each segment
------- For now I did not work on introducing this feature. n_seg_xy is a place holder, and Point A and B are the coordinate of two points that define the slab trench
type: Describe the type of bending: for now we have only :Linear and :Ribe 
L0  : the total length of the slab; 
D0  : the total thickness of the slab; 
n_seg: the total amount of segment in which the slab is discretized along its length 
WZ  : the thickness of the weakzone {The option for making the weak zone is not there yet. }
Lb : an additional parameter concerning the bending angle. It tells to the function at which length the angle of subduction becomes constant. 
d_decoupling: it represents the depth at which both of the surface of the slab are submerged in the upper mantle. 
1. -> The function starts by extracting the relevant information from the structure. 
    a. Compute the top, bottom and weak zone surface using: θ_max, L0, n_seg, D0, Lb 
    b. Find the slab: 
        b1. Transform the coordinate such that the new x coordinate is parallel to the segment A-B. {translate the origin to the coordinate of A, and rotate the coordinate anticlockwise}
           -> The sign of θ_max determine the dip direction of the slab w.r.t. the new coordinate system. 
        b2. Loop over the segment of the slab for creating a small poligon and to find the points using the transformed coordinate system
        b3. Interpolate d {distance from the top surface} and l {lenght} of the slab
    c. Find the length at which the slab reaches the decoupling depth
    d. Compute the temperature field 
    e. Compute the phase field

"""
    # number of cells in every direction
    nx = 256;

    ny = nx;

    nz = nx;


    # define domain size
    x        = LinRange(0.0,1200.0,nx);

    y        = LinRange(0.0,1200.0,ny);

    z        = LinRange(-660,50,nz);

    X,Y,Z    = XYZGrid(x, y, z);

    Cart     = CartData(X,Y,Z, (Data=Z,));

    # initialize phase and temperature matrix
    Phase   = ones(Int32,size(X));

    Temp    = ones(Float64,size(X))*1350;

    # add different phases: crust->2, Mantle Lithosphere->3 Mantle->1
    AddBox!(Phase, Temp, Cart; xlim=(0.0,1200.0),ylim=(0.0,1200.0), zlim=(-660.0,0.0), phase = LithosphericPhases(Layers=[30 400 600], Phases=[2 3 1], Tlab=1300 ), T=HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4) )

    # add air phase 0
    AddBox!(Phase, Temp, Cart; xlim=(0.0,1200.0),ylim=(0.0,1200.0), zlim=(0.0,50.0), phase = ConstantPhase(0), T=ConstantTemp(20.0))

    # Create the first trench structre
    t_1 = Trench(1,[200.0,400.0],[500.0,700.0],90.0,:Ribe,50,500.0,100.0,20.0,200.0,100.0);

    # Create the second trench structre
    t_2 = Trench(1,[200.0,400.0],[500.0,700.0],-30.0,:Ribe,50,500.0,100.0,20.0,400.0,100.0);

    # Stratigraphy of the slab
    stratigraphy_slab = LithosphericPhases(Layers=[10 100], Phases=[5 5 5], Tlab=1300 );

    #Temperature field of the first slab
    temperature_slab=McKenzie_subducting_slab(20.0,1350.0,30.0,0.4,2.5,1050.0,3.0,3300.0,36);

    #Temperature field of the second slab
    temperature_slab2=McKenzie_subducting_slab(20.0,1350.0,30.0,0.4,10.0,1050.0,3.0,3300.0,36);

    #Create the first slab
    create_slab!(X,Y,Z,Phase,Temp,t_1,stratigraphy_slab,temperature_slab);
    # Create the second slab 
    create_slab!(X,Y,Z,Phase,Temp,t_2,stratigraphy_slab,temperature_slab2);

    # Save data to paraview:
    Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp)) 

    Write_Paraview(Data_Final, "Double_Subduction")


#function Benchmark_trench()
#
#    nx = 128;
#
#    ny = nx;
#
#    nz = nx;
#
#    # define domain size
#    x        = LinRange(0.0,1200.0,nx);
#
#    y        = LinRange(0.0,1200.0,ny);
#
#    z        = LinRange(-660,50,nz);
#
#    X,Y,Z    = XYZGrid(x, y, z);
#
#    Cart     = CartData(X,Y,Z, (Data=Z,));
#
#    # initialize phase and temperature matrix
#    Phase   = ones(Int32,size(X));
#
#    Temp    = ones(Float64,size(X))*1350;
#    
#    t_ = Trench(1,[200.0,400.0],[500.0,700.0],90.0,:Ribe,50,500.0,100.0,20.0,200.0,100.0);
#
#    # Benchmarking each function of the routine
#
#    # -> d = distance from the top surface
#    d = ones(size(X)).*NaN64;
#
#    # -> l = length from the trench along the slab 
#    ls = ones(size(X)).*NaN64;
#
#    stratigraphy_slab = LithosphericPhases(Layers=[10 90], Phases=[2 3 1], Tlab=1300 );
#
#    temperature_slab=McKenzie_subducting_slab(20,1350,30,0.4,2.5,1050,3.0,3300,36);
#
#    XT = zeros(size(X));
#
#    YT = zeros(size(Y)); 
#
#    # Function to transform the coordinate 
#    # Most Problematic function 
#    xb = transform_coordinate!($X,$Y,$XT,$YT,$t_.A,$t_.B,sign($t_.theta_max)); 
#    
#    @code_warntype (transform_coordinate!(X,Y,Z,XT,YT,t_.A,t_.B,sign(t_.theta_max)))
#    @btime (transform_coordinate!($X,$Y,$Z,$XT,$YT,$t_.A,$t_.B,sign($t_.theta_max)))
#    # Rotate coordinate: 
#    Rot3D!(XT,YT,Z, 40, 0);
#    @btime (Rot3D!($XT,$YT,$Z, $40, $0))
#    # dl 
#
#    @code_warntype (compute_slab_surface!(t_.D0,t_.L0,t_.Lb,t_.WZ,t_.n_seg,abs(t_.theta_max),t_.type_bending))
#    @btime (compute_slab_surface!($t_.D0,$t_.L0,$t_.Lb,$t_.WZ,$t_.n_seg,abs($t_.theta_max),$t_.type_bending))
#
#    @btime (compute_bending_angle!(abs($t_.theta_max),$t_.Lb,10.0,$t_.type_bending))
#    @code_warntype (compute_bending_angle!(abs(t_.theta_max),t_.Lb,10.0,t_.type_bending))
#
#    Top,Bottom = compute_slab_surface!(t_.D0,t_.L0,t_.Lb,t_.WZ,t_.n_seg,abs(t_.theta_max),t_.type_bending)
#    @code_warntype(find_slab!(X,Y,Z,d,ls,t_.theta_max,t_.A,t_.B,Top,Bottom,t_.n_seg,t_.D0,t_.L0))
#    @btime(find_slab!($X,$Y,$Z,$d,$ls,$t_.theta_max,$t_.A,$t_.B,$Top,$Bottom,$t_.n_seg,$t_.D0,$t_.L0))
#
#
#    @btime (create_slab!($X,$Y,$Z,$Phase,$Temp,$t_,$d,$ls,$stratigraphy_slab,$temperature_slab))
#    @code_warntype(create_slab!(X,Y,Z,Phase,Temp,t_,d,ls,stratigraphy_slab,temperature_slab))
#
#    Test_Double_subduction()
#    @btime (Test_Double_subduction())
#    
#end
#
#
#