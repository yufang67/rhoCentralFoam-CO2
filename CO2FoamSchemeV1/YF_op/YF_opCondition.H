// Input OP conditions
Info<< "\nSet YF Operating conditions\n" << endl;
// Initialisation
// Inlet
   scalar e_inlet1 = -218360.5;
   scalar rho_inlet1 = 663.64;
   scalar e_inlet  = e_inlet1;
   scalar rho_max  = rho_inlet1;

   scalar e_inlet2   = 0.0;
   scalar rho_inlet2 = 0.0;
   scalar p_inlet2   = 0.0;
   scalar T_inlet2   = 0.0;
// Outlet
   scalar e_outlet= 384.42e+03; 
   scalar rho_min = 148.41;
// total length
   scalar ll      = 35.0e-3;
//
//   
// BC condition
   scalar p_inlet1  = 10.00001e6;
   scalar T_inlet1  = 311.33;
//   scalar rho_in= rho_inlet1;//485.5;//621.475;
//   scalar e_in  = e_inlet;//-181535.542;//-211417.86;
// taget pressure 
   scalar p_out = 5.0001e6;
//// characteristic length
   scalar lx    = ll;
//// relaxation factor
   scalar relax = 1;
