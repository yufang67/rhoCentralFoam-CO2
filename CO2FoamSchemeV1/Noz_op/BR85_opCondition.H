// Input OP conditions
Info<< "\nSet Bera Operating conditions\n" << endl;
// Initialisation
// Inlet
//
   scalar e_inlet1 = -152235.21;
   scalar rho_inlet1 = 353.9;
   scalar e_inlet  = e_inlet1;
   scalar rho_max  = rho_inlet1;

   scalar e_inlet2   = 0.0;
   scalar rho_inlet2 = 0.0;
   scalar p_inlet2   = 0.0;
   scalar T_inlet2   = 0.0;
// Outlet
   scalar e_outlet= 328.78e+03; 
   scalar rho_min = 389.85;     
// total length
   scalar ll      = 35e-3;//26e-3;//35e-3;
//
//   
// BC condition
   scalar p_inlet1  = 8.5e6;
   scalar T_inlet1  = 313.15; 
// taget pressure 
   scalar p_out = 5.5e6; 
//// characteristic length
   scalar lx    = ll;
//// relaxation factor
   scalar relax = 10;
