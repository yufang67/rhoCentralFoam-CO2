// Input OP conditions
Info<< "\nSet Naka3a_9MPa Operating conditions\n" << endl;
// Initialisation
// Inlet
   scalar e_inlet1 = -218491.9;
   scalar rho_inlet1 = 651.74;
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
   scalar ll      = 50.0e-3;
//
//   
// BC condition
   scalar p_inlet1  = 9.10001e6;
   scalar T_inlet1  = 308.95;
//   scalar rho_in= rho_inlet1;//485.5;//621.475;
//   scalar e_in  = e_inlet;//-181535.542;//-211417.86;
// taget pressure 
   scalar p_out = 1.0e6;
//// characteristic length
   scalar lx    = ll;
//// relaxation factor
   scalar relax = 1000;
