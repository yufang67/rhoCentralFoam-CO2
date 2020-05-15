// Input OP conditions
Info<< "\nSet ejector Operating conditions 9MPa\n" << endl;
// Initialisation
// Inlet
   scalar e_inlet1 = -254224.9;
   scalar rho_inlet1  = 784.6;

   scalar e_inlet2 = -122636.6;
   scalar rho_inlet2  = 177.32;

   scalar e_mid  = e_inlet2;
   scalar rho_mid  = rho_inlet2;

// Outlet
   scalar e_outlet= -125.7982e3;
   scalar rho_outlet = 187.46;
// total length
   scalar ll1      = 19.302e-3;
   scalar ll2      = 159.698e-3;
   scalar yl      = 3.0e-3;
//
//   
// BC condition
   scalar p_inlet1  = 9.0e6;
   scalar T_inlet1  = 299.65;

   scalar p_inlet2  = 5.5523e6;
   scalar T_inlet2  = 293.13;
// taget pressure 
   scalar p_out =  5.6e6;
//// characteristic length
   scalar lx    = ll2;
//// relaxation factor
   scalar relax = 50;
