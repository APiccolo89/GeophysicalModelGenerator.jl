

using GeophysicalModelGenerator


# number of cells in every direction
nx = 100
ny = 100
nz = 200

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
# Layers present thickness and get add up
AddBox!(Phase, Temp, Cart; xlim=(0.0,800.0),ylim=(0.0,800.0), zlim=(-800.0,0.0), phase = LithosphericPhases(Layers=[15 30 100 800], Phases=[2 3 1 5], Tlab=1300 ), T=LinearTemp(Ttop=20, Tbot=1600) )#T=HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4)


# add sediment basin # depending on resolution show and angle it is not shown in paraview
AddPolygon!(Phase, Temp, Cart; xlim=(0.0,0.0, 160.0, 200.0),ylim=(0.0,800.0), zlim=(0.0,-10.0,-20.0,0.0), phase = ConstantPhase(8), T=LinearTemp(Ttop=20, Tbot=30))

# add thinning of the continental crust attached to the slab and its thickness for a smooth transition
AddPolygon!(Phase, Temp, Cart; xlim=(0.0, 200.0, 0.0),ylim=(0.0,800.0), zlim=(-100.0,-150.0,-150.0), phase = ConstantPhase(5), T=LinearTemp(Ttop=1000, Tbot=1100))

# add accretionary prism 
AddPolygon!(Phase, Temp, Cart; xlim=(800.1, 600.0, 800.1),ylim=(0.0,800.0), zlim=(0.0,0.0,-60.0), phase = ConstantPhase(8), T=LinearTemp(Ttop=20, Tbot=30))


# add air phase 0
AddBox!(Phase, Temp, Cart; xlim=(0.0,800.0),ylim=(0.0,800.0), zlim=(0.0,50.0), phase = ConstantPhase(0), T=ConstantTemp(20.0))

# # Save data to paraview:
Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp)) 
Write_Paraview(Data_Final, "Sedimentary_basin")


