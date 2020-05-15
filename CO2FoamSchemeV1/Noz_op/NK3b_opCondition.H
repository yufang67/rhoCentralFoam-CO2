// Input OP conditions
Info<< "\nSet Naka3b_7MPa Operating conditions\n" << endl;
// Initialisation
// Inlet
   scalar e_inlet1 = -240311.4;
   scalar rho_inlet1 = 716.0;
   scalar e_inlet  = e_inlet1;
   scalar rho_max  = rho_inlet1;

   scalar e_inlet2   = 0.0;
   scalar rho_inlet2 = 0.0;
   scalar p_inlet2   = 0.0;
   scalar T_inlet2   = 0.0;
// Outlet
   scalar e_outlet= 393.31e+03; 
   scalar rho_min = 51.94;
// total length
   scalar ll      = 100.0e-3;
//
//   
// BC condition
   scalar p_inlet1  = 7.10001e6;
   scalar T_inlet1  = 299.85;
//   scalar rho_in= rho_inlet1;//485.5;//621.475;
//   scalar e_in  = e_inlet;//-181535.542;//-211417.86;
// taget pressure 
   scalar p_out = 0.9e6;
//// characteristic length
   scalar lx    = ll;
//// relaxation factor
   scalar relax = 100;
