using GeophysicalModelGenerator
using ScatteredInterpolation 
using BenchmarkTools


    # number of cells in every direction
    nx = 128;

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
    AddBox!(Phase, Temp, Cart; xlim=(0.0,1200.0),ylim=(0.0,1200.0), zlim=(-660.0,0.0), phase = LithosphericPhases(Layers=[30 400 600], Phases=[2 3 1], Tlab=1300 ), T=HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4) );

    # add air phase 0
    AddBox!(Phase, Temp, Cart; xlim=(0.0,1200.0),ylim=(0.0,1200.0), zlim=(0.0,50.0), phase = ConstantPhase(0), T=ConstantTemp(20.0));

    TsHC = HalfspaceCoolingTemp(Tsurface=20.0, Tmantle=1350, Age=120, Adiabat=0.4)
    TsMK = McKenzie_subducting_slab(Tsurface = 20.0, Tmantle = 1350.0,v_s = 4.0, adiabat = 0.4,Cp = 1050.0, k = 3.0, it = 3.0)
    T_slab = LinearWeightedTemperature(0.1,1.0,100.0,:X,TsMK,TsHC);

    AddBox!(Phase, Temp, Cart; xlim=(0.0,600.0),ylim=(0.0,600.0), zlim=(0.0,-80.0),StrikeAngle=0, DipAngle=45, phase = ConstantPhase(5), T=McKenzie(20.0));

    Data_Final      =   CartData(X,Y,Z,(Phase=Phase,Temp=Temp)) 

    Write_Paraview(Data_Final, "Double_Subduction")
