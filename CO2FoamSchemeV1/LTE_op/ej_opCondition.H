// Input OP conditions
Info<< "\nSet ejector Operating conditions 21_03_2018 \n" << endl;
// Initialisation
// Inlet
   scalar e_inlet1 = -125399.7;  //-122740.2;
   scalar rho_inlet1  = 246.97;  //235.21; 

   scalar e_inlet2 = -118024.07; //-118546.9;
   scalar rho_inlet2  = 140.55;  //141.04; 

   scalar e_mid  = e_inlet2; //-125798.2;
   scalar rho_mid  = rho_inlet2;//187.46;

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
   scalar p_inlet1  = 7.58755e6; //7.4001e6;
   scalar T_inlet1  = 311.88;    //311.18;

   scalar p_inlet2  = 4.63749e6; //4.575e6;
   scalar T_inlet2  = 284.42;    //284.37;
// taget pressure 
   scalar p_out     = 5.6788e6;     //5.6e6;
//// characteristic length
   scalar lx    = ll2;
//// relaxation factor
   scalar relax = 50;
