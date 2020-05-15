// Input OP conditions
Info<< "\nSet JW Operating conditions\n" << endl;
// Initialisation
// Inlet
   scalar e_inlet1 = -54785.1;
   scalar rho_inlet1 = 16.53;
   scalar e_inlet  = e_inlet1;
   scalar rho_max  = rho_inlet1;

   scalar e_inlet2   = 0.0;
   scalar rho_inlet2 = 0.0;
   scalar p_inlet2   = 0.0;
   scalar T_inlet2   = 0.0;
// Outlet
   scalar e_outlet= 439.0e+03; 
   scalar rho_min = 13.747;
// total length
   scalar ll      = 30.0e-3;
//
//   
// BC condition
   scalar p_inlet1  = 0.925e6;
   scalar T_inlet1  = 309.15;
//   scalar rho_in= rho_inlet1;//485.5;//621.475;
//   scalar e_in  = e_inlet;//-181535.542;//-211417.86;
// taget pressure 
   scalar p_out = 0.615e6;
//// characteristic length
   scalar lx    = ll;
//// relaxation factor
   scalar relax = 1;
