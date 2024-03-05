using GeophysicalModelGenerator
include("../src/trench.jl")
using Plots
using ScatteredInterpolation 


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

t_ = Trench(1,(200.0,400.0),(500.0,700.0),90.0,"Ribe",50,500.0,100.0,20.0,200.0,100.0)

t_2 = Trench(1,(200.0,400.0),(500.0,700.0),-30.0,"Ribe",50,500.0,100.0,20.0,200.0,100.0)
t_2.Lb = 400.0

stratigraphy_slab = LithosphericPhases(Layers=[10 90], Phases=[2 3 1], Tlab=1300 )

temperature_slab=McKenzie_subducting_slab(20,1350,30,0.4,2.5,1050,3.0,3300,36);

temperature_slab2=McKenzie_subducting_slab(20,1350,30,0.4,10.0,1050,3.0,3300,36);



d,l = create_slab!(X,Y,Z,Phase,Temp,t_,stratigraphy_slab,temperature_slab);
d2,l2 = create_slab!(X,Y,Z,Phase,Temp,t_2,stratigraphy_slab,temperature_slab2);


Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp,d,l,d2,l2)) 
Write_Paraview(Data_Final, "Subduction_example")

# add different phases: crust->2, Mantle Lithosphere->3 Mantle->1
AddBox!(Phase, Temp, Cart; xlim=(0.0,800.0),ylim=(0.0,800.0), zlim=(-660.0,0.0), phase = LithosphericPhases(Layers=[30 400 600], Phases=[2 3 1], Tlab=1300 ), T=HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4) )

# add air phase 0
AddBox!(Phase, Temp, Cart; xlim=(0.0,800.0),ylim=(0.0,800.0), zlim=(0.0,50.0), phase = ConstantPhase(0), T=ConstantTemp(20.0))

# Save data to paraview:
Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp)) 
Write_Paraview(Data_Final, "Subduction_Setup")





