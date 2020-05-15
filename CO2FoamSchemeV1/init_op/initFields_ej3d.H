Info<< "\nCreating CO2 table ...\n" << endl;
//    //    Info<< "\n---------pp=" << c << endl;
    __grid_MOD_make_grid();
//    grid_construction_left_low_();
//    grid_construction_left_high_();
//    grid_construction_high_temperature_();
//    __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,
//                                    &ee,&vv,&flag);
//    Info<< "\n pp,TT,cc" << pp << " " << TT << " " << cc << endl;
//    Info<< "\n xx,a,flag" << x_out<< " " << a_out << " " << flag << endl;
//*************************************************************************
Info<< "\nInitialize conservative fields\n" << endl;
//
double ee=0.0; //-99772;
double vv=0.0; //1.0/128.4;
double pp=0.0;
double TT=0.0;
double cc=0.0;
double x_out;
double a_out;
double dummy;
int    flag;
double vis;
const double e_ref =5.0e5;
scalar initial_time;
scalar M_max = 0.0;
initial_time = runTime.time().value();
//Info <<" time= " <<runTime.timeName()<<endl;
//Info <<" time= " <<initial_time<<endl;
//
if (initial_time < 0.00000001) 
{
Info << "smooth initial field" <<endl;
   const volVectorField&     center = mesh.C();
   scalar xcell;
   scalar ycell;
   scalar zcell;

   scalar e_in1 = e_inlet1 + e_ref;
   scalar e_in2 = e_inlet2 + e_ref;
   scalar e_m   = e_mid    + e_ref;
   scalar e_out = e_outlet + e_ref;

   forAll(e,celli) 
   {
           xcell = center[celli].x();
           ycell = center[celli].y();
           zcell = center[celli].z();

       if (xcell <= ll1) 
       {
	   if (zcell > 3e-3 or ycell > 3e-3 or ycell <-3e-3) {

              rho[celli] = rho_inlet2 + xcell/ll1 * (rho_mid-rho_inlet2);
              e[celli] = e_in2 - xcell/ll1 * (e_in2-e_m);
//
              ee = e[celli] - e_ref;
              vv = 1.0/rho[celli];
              __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,
                                           &ee,&vv,&flag);

              p[celli]  = pp;
              T[celli]  = TT;
              rhoU[celli] = rho[celli]*U[celli];
              rhoE[celli] = rho[celli]*(e[celli] + 0.5*magSqr(U[celli]));
	                     }
	   else{

              rho[celli] = rho_inlet1 - xcell/ll1 * (rho_inlet1-rho_mid);
//              e[celli] = e_in1 - xcell/ll1 * (e_in1-e_m);
              e[celli] = e_in1 + xcell/ll1 * (e_m-e_in1);
//
              ee = e[celli] - e_ref;
              vv = 1.0/rho[celli];
              __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,
                                           &ee,&vv,&flag);

              p[celli]  = pp;
              T[celli]  = TT;
              rhoU[celli] = rho[celli]*U[celli];
              rhoE[celli] = rho[celli]*(e[celli] + 0.5*magSqr(U[celli]));
	        }
       } 
       else 
       {
              rho[celli] = rho_mid + xcell/ll2 * (rho_outlet-rho_mid);
              e[celli]   = e_m - xcell/ll2 * (e_m-e_out);;
//
              ee = e[celli] - e_ref;
              vv = 1.0/rho[celli];
              __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,
                                           &ee,&vv,&flag);

              p[celli]  = pp;
              T[celli]  = TT;
              rhoU[celli] = rho[celli]*U[celli];
              rhoE[celli] = rho[celli]*(e[celli] + 0.5*magSqr(U[celli]));

       }
   }
}
//
//forAll(mesh.boundary(), patchID)
//        {
//         forAll (mesh.boundary()[patchID],facei)
//            {
//               TT = T.boundaryField()[patchID][facei];
//               vv = 1.0/rho.boundaryField()[patchID][facei];
//
//               __properties_MOD_inter_energy(&TT,&vv, &ee);
//
//               e.boundaryFieldRef()[patchID][facei] = ee + e_ref;
//            }
//        }
//
//rhoU = rho * U;
//rhoE = rho*(e + 0.5*magSqr(U));
//
//Info<< "\n e.boundary\n" << e.boundaryFieldRef() << endl;
//Info<< "\n rho\n" << rho << endl;
//Info<< "\n U\n"   << U   << endl;
//Info<< "\n e\n"   << e   << endl; 
