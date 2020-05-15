// Input OP conditions
// Inlet
   scalar e_inlet1 = -218360.5;
   scalar rho_inlet1  = 663.64;

   scalar e_inlet2 = -90814.3;
   scalar rho_inlet2  = 92.18;

   scalar e_mid  = e_inlet2; //-125798.2;
   scalar rho_mid  = rho_inlet2;//187.46;

// Outlet
   scalar e_outlet= -111.66e3;
   scalar rho_outlet = 113.08;
// total length
   scalar ll1      = 29.36e-3;
   scalar ll2      = 125.25e-3;
   scalar yl      = 4.0e-3;
//
//   
// BC condition
   scalar p_inlet1  = 10.0001e6;
   scalar T_inlet1  = 311.33;

   scalar p_inlet2  = 3.91e6;
   scalar T_inlet2  = 295.82;
// taget pressure 
   scalar p_out = 4.3401e6;
//// characteristic length
   scalar lx    = ll2;
//// relaxation factor
   scalar relax = 10;

//
Info<< "\n***************YF_ejector Operating Condition******************\n" << endl;
Info<< "\nPrimary inlet: p= "<<p_inlet1<<" T= "<<T_inlet1<<"K\n" << endl;
Info<< "\nSecondary inlet: p= "<<p_inlet2<<" T= "<<T_inlet2<<"K\n" << endl;
Info<< "\nOutlet: p= "<<p_out<<" relax= "<<relax<< endl;
Info<< "\n****************************************************************\n" << endl;
