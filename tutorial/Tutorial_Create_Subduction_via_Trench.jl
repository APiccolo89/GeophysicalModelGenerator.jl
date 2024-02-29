include("../src/trench.jl")
using GeophysicalModelGenerator
using Plots


# number of cells in every direction
nx = 512;
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
t_2 = Trench(1,(200.0,400.0),(500.0,700.0),-90.0,"Ribe",50,500.0,100.0,20.0,200.0,100.0)


Top,MidS,Bottom,WZ_surf,theta_mean = compute_slab_surface!(t_.D0,t_.L0,t_.Lb,t_.WZ,t_.n_seg,abs(t_.theta_max),t_.type_bending);

d = zeros(size(X));

l = zeros(size(X));


XT,YT,Z,d,l=find_slab!(X,Y,Z,d,l,t_.theta_max,t_.A,t_.B,Top,Bottom,t_.n_seg,t_.D0,t_.L0)

#XT = zeros(size(X));

#YT = zeros(size(Y)); 

# Function to transform the coordinate 
#XT,YT, xb = transform_coordinate!(X,Y,XT,YT,t_.A,t_.B,sign(t_.theta_max)); 
#p = plot()
#for i= 1:t_.n_seg;
#    @show i
#pa = (Top[i,1],Top[i,2]); # D = 0 | L = l 
#
#pb = (Bottom[i,1],Bottom[i,2]); # D = -D0 | L=l 
#
#pc = (Bottom[i+1,1],Bottom[i+1,2]); # D = -D0 |L=L+dl
#
#pd = (Top[i+1,1],Top[i+1,2]); # D = 0| L = L+dl
#
## Create the polygon 
#poly_y = [pa[1],pb[1],pc[1],pd[1]];
#poly_z = [pa[2],pb[2],pc[2],pd[2]];
#@show poly_y
#@show poly_z
#
#ymin = minimum(poly_y);
#ymax = maximum(poly_y);
#zmin = minimum(poly_z);
#zmax = maximum(poly_z);
#ind_chosen = findall(0.0.<= XT.<= xb[1] .&& ymin .<= YT .<= ymax .&& zmin .<= Z .<= zmax);
#ind_seg = []
#yp = YT[ind_chosen];
#zp = Z[ind_chosen];
#
#ind = zeros(Bool,size(zp));
#
#inPoly!(poly_y,poly_z,yp,zp,ind); 
#ind_prophet = ind_chosen[ind]
#d[ind_prophet] .=1.0;
#@show
#
#if i ==1 
#    scatter!(p,Y[ind_chosen][ind],Z[ind_chosen][ind],label=:none)
#    plot!(p,poly_y,poly_z,label=:none)
#else
#    scatter!(p,Y[ind_chosen][ind],Z[ind_chosen][ind],label=:none)
#    plot!(p,poly_y,poly_z,label=:none)
#
#end
#display(p)
#
#end
#
Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp,d)) 
Write_Paraview(Data_Final, "Bla_bla_car")

# add different phases: crust->2, Mantle Lithosphere->3 Mantle->1
AddBox!(Phase, Temp, Cart; xlim=(0.0,800.0),ylim=(0.0,800.0), zlim=(-660.0,0.0), phase = LithosphericPhases(Layers=[30 400 600], Phases=[2 3 1], Tlab=1300 ), T=HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4) )

# add air phase 0
AddBox!(Phase, Temp, Cart; xlim=(0.0,800.0),ylim=(0.0,800.0), zlim=(0.0,50.0), phase = ConstantPhase(0), T=ConstantTemp(20.0))

# Save data to paraview:
Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp)) 
Write_Paraview(Data_Final, "Subduction_Setup")





