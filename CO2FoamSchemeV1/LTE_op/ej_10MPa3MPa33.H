// Input OP conditions
Info<< "\nSet ejector Operating conditions 10MPa_3MPa_4MPa\n" << endl;
// Initialisation
// Inlet
   scalar e_inlet1 = -257459.1;
   scalar rho_inlet1  = 801.6;

   scalar e_inlet2 = -109116.5;
   scalar rho_inlet2  = 81.4;

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
   scalar p_inlet1  = 10.0e6;
   scalar T_inlet1  = 300.0;

   scalar p_inlet2  = 3.0e6;
   scalar T_inlet2  = 268.15;
// taget pressure 
   scalar p_out =  3.3e6;
//// characteristic length
   scalar lx    = ll2;
//// relaxation factor
   scalar relax = 10;
