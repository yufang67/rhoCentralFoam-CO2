// Input OP conditions
Info<< "\nSet nozzle Operating conditions\n" << endl;
// Initialisation
// Inlet
   scalar e_inlet1 = -171254.9; //-195763.1;//-194096.0; // -194485.5; //-194096.0;//-191859.7;//-193626.97;
   scalar rho_inlet1 = 435.88;  //556.24; //547.7;     //550.54;    //547.7;//537.7;//545.57;
   scalar e_inlet  = e_inlet1;
   scalar rho_max  = rho_inlet1;

   scalar e_inlet2   = 0.0;
   scalar rho_inlet2 = 0.0;
   scalar p_inlet2   = 0.0;
   scalar T_inlet2   = 0.0;
// Outlet
   scalar e_outlet= 328.78e+03; //382.47e3;
   scalar rho_min = 389.85;     //182.31;//148.41;//182.31;
// total length
   scalar ll      = 35e-3;//26e-3;//35e-3;
//
//   
// BC condition
   scalar p_inlet1  = 8.801e6; //9.301e6;  //9.25e6;//9.0e6;
   scalar T_inlet1  = 313.0;   //313.15; //313.0;//310.15;
//   scalar rho_in= rho_inlet1;//485.5;//621.475;
//   scalar e_in  = e_inlet;//-181535.542;//-211417.86;
// taget pressure 
   scalar p_out = 7.0e6; //7.5e6;//5.9e6;
//// characteristic length
   scalar lx    = ll; 37e-3;//25e-3;//100e-3;//80e-3;
//// relaxation factor
   scalar relax = 1;
