using GeophysicalModelGenerator

# number of cells in every direction
nx = 64
ny = 64
nz = 64

# define domain size
x        = LinRange(0.0,800.0,nx)
y        = LinRange(0.0,800.0,ny)
z        = LinRange(-660,50,nz)
X,Y,Z    = XYZGrid(x, y, z);
Cart     = CartData(X,Y,Z, (Data=Z,))

# initialize phase and temperature matrix
Phase   = ones(Int32,size(X))
Temp    = ones(Float64,size(X))*1350

# add different phases: crust->2, Mantle Lithosphere->3 Mantle->1
AddBox!(Phase, Temp, Cart; xlim=(0.0,800.0),ylim=(0.0,800.0), zlim=(-660.0,0.0), phase = LithosphericPhases(Layers=[30 400 600], Phases=[2 3 1], Tlab=1300 ), T=HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4) )

# add air phase 0
AddBox!(Phase, Temp, Cart; xlim=(0.0,800.0),ylim=(0.0,800.0), zlim=(0.0,50.0), phase = ConstantPhase(0), T=ConstantTemp(20.0))

# Save data to paraview:
Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp)) 
Write_Paraview(Data_Final, "Subduction_Setup")





