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
if (initial_time < 0.00000001) {
Info << "smooth initial field" <<endl;
   const volVectorField&     center = mesh.C();
   scalar xcell;
   scalar e_min   = e_inlet  + e_ref;
   scalar e_max   = e_outlet;
   forAll(e,celli) {
           xcell = center[celli].x();
           rho[celli] = rho_max - xcell/ll * (rho_max-rho_min);
           e[celli] = e_min + xcell/ll * (e_max-e_min);
//
           ee = e[celli] - e_ref;
           vv = 1.0/rho[celli];
           __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&dummy,
                                           &ee,&vv,&flag);

           p[celli]  = pp;
           T[celli]  = TT;
           rhoU[celli] = rho[celli]*U[celli];
           rhoE[celli] = rho[celli]*(e[celli] + 0.5*magSqr(U[celli]));
//         TT = T[celli];
//
//         __properties_MOD_inter_energy(&TT,&vv, &ee);
         
//         e[celli] = ee + e_ref;
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
