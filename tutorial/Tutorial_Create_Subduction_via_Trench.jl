include("../src/trench.jl")
using GeophysicalModelGenerator


# number of cells in every direction
nx = 64;
ny = nx+1;
nz = nx+2;


# define domain size
x        = LinRange(0.0,1200.0,nx);
y        = LinRange(0.0,1200.0,ny);
z        = LinRange(-660,50,nz);

X,Y,Z    = XYZGrid(x, y, z);

Cart     = CartData(X,Y,Z, (Data=Z,));

# initialize phase and temperature matrix
Phase   = ones(Int32,size(X));

Temp    = ones(Float64,size(X))*1350;

t_ = Trench(1,(200.0,400.0),(500.0,700.0),45.0,"Ribe",50,300.0,100.0,20.0,200.0,100.0)

A,B,C,D,E = compute_slab_surface!(t_.D0,t_.L0,t_.Lb,t_.WZ,t_.n_seg,abs(t_.theta_max),t_.type_bending);

XT = zeros(size(X));

YT = zeros(size(Y));




# add different phases: crust->2, Mantle Lithosphere->3 Mantle->1
AddBox!(Phase, Temp, Cart; xlim=(0.0,800.0),ylim=(0.0,800.0), zlim=(-660.0,0.0), phase = LithosphericPhases(Layers=[30 400 600], Phases=[2 3 1], Tlab=1300 ), T=HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4) )

# add air phase 0
AddBox!(Phase, Temp, Cart; xlim=(0.0,800.0),ylim=(0.0,800.0), zlim=(0.0,50.0), phase = ConstantPhase(0), T=ConstantTemp(20.0))

# Save data to paraview:
Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp)) 
Write_Paraview(Data_Final, "Subduction_Setup")





