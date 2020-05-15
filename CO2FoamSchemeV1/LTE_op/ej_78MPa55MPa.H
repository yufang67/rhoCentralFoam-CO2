// Input OP conditions
Info<< "\nSet ejector Operating conditions 7.8MPa_5.5MPa\n" << endl;
// Initialisation
// Inlet
   scalar e_inlet1 = -218868.6;
   scalar rho_inlet1  = 631.0;

   scalar e_inlet2 = -107586.9;
   scalar rho_inlet2  = 151.9;

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
   scalar p_inlet1  = 7.8e6;
   scalar T_inlet1  = 305.0;

   scalar p_inlet2  = 5.5e6;
   scalar T_inlet2  = 300.0;
// taget pressure 
   scalar p_out =  5.6e6;
//// characteristic length
   scalar lx    = ll2;
//// relaxation factor
   scalar relax = 50;
