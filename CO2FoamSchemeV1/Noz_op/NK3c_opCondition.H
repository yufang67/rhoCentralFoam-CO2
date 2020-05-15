// Input OP conditions
Info<< "\nSet Naka3c_6MPa Operating conditions\n" << endl;
// Initialisation
// Inlet
   scalar e_inlet1 = -253593.08;
   scalar rho_inlet1 = 758.667;
   scalar e_inlet  = e_inlet1;
   scalar rho_max  = rho_inlet1;
   scalar ux_min  =5.0;

   scalar e_inlet2   = 0.0;
   scalar rho_inlet2 = 0.0;
   scalar p_inlet2   = 0.0;
   scalar T_inlet2   = 0.0;
// Outlet
   scalar e_outlet= 244.45e+03; //438.235e+03; 
   scalar rho_min = 52.92;      //18.579;
   scalar ux_max  =100.0;
// total length
   scalar ll      = 100.0e-3;
//
//   
// BC condition
   scalar p_inlet1  = 6.10001e6;
   scalar T_inlet1  = 294.95; //21.8degreeC
// taget pressure 
   scalar p_out = 1e6;
//// characteristic length
   scalar lx    = ll;
//// relaxation factor
   scalar relax = 2000;
