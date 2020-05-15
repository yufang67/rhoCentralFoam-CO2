//
// CO2 table for pressure and speed of sound (on vol)
//
const volVectorField&     center = mesh.C();

//double inv_rhov;
//double  ev ;
        forAll(e,celli)
        {
         scalar xxx = center[celli].x(); //center[celli].x();
//Info<<" ee= "<<e[celli]<<" rho= "<<rho[celli]<<" i= "<<celli<<" x= "<<xxx <<endl;
         ee = e[celli] - e_ref;
         vv = 1.0/rho[celli];
	 
//
//Info<<" ee= "<<ee<<" vv= "<<vv<<" i= "<<celli<<" x= "<<xxx <<endl;
         __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&xxx,&ee,&vv,&flag); 
//
//         press[celli] = pp;
//         temp[celli]  = TT;
         p[celli]  = pp;
         T[celli]  = TT;
         c[celli]  = cc;
	 qual[celli]= x_out;
	 volfrac[celli]= a_out;
//	 muBulk[celli] = 0.125*pp*4e-6;
//       
//         if (flag==5) {
//	 inv_rhov      = (a_out*vv)/x_out;
//	 __properties_MOD_inter_energy(&TT, &inv_rhov, &ev);
//         hv[celli] = (ev + pp*inv_rhov);
//	 hl[celli] = (ee+pp*vv-x_out*(ev + pp*inv_rhov) ) / (1-x_out);
//	              }
//	 else{
//                  hv[celli] = 0;
//		  hl[celli] = 0;
//             }
//
//Info<<"flag"<<flag<<" i "<<celli <<endl;
         __transprop_MOD_cpco2(&cp_out, &vv, &vp_out, &x_out, &TT, &pp, &flag);
         cp[celli] = cp_out;

         __transprop_MOD_cvco2(&cv_out, &vv, &vp_out, &x_out, &TT, &pp, &flag);
         cv[celli] = cv_out;

	 __transprop_MOD_co2visco(&vis, &vv, &TT, &x_out, &pp, &flag);
         mu[celli] = vis;
  
         __transprop_MOD_co2conduc2phase(&lambda_out, &vv, &vp_out, &x_out, &TT, &pp, &flag);
         lambda_l[celli] = lambda_out;

         __transprop_MOD_dedrco2(&dedr_out, &vv, &vp_out, &x_out, &TT, &pp, &flag);
         dedr[celli] = dedr_out;




//         nu[celli] = vis*vv;
//         alpha_lam[celli] = 0.001;
//         Info<< "\n pp" << pp  << endl;
//Info<< "\n TT" << TT  << endl;
	}
        forAll(mesh.boundary(), patchID)
	{
		forAll (mesh.boundary()[patchID],facei)
                {
	                 ee = e.boundaryField()[patchID][facei] - e_ref;
                         vv = 1.0/rho.boundaryField()[patchID][facei];
//Info<< "\n ee bc" << ee  <<" "<<patchID<<" "<<facei <<endl;
//
			 scalar xxx = center[facei].x(); //center[celli].x();
                         __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&xxx,&ee,&vv,&flag);
//
			 p.boundaryFieldRef()[patchID][facei]  = pp;
                	 T.boundaryFieldRef()[patchID][facei]  = TT;
		         c.boundaryFieldRef()[patchID][facei]  = cc;
			 qual.boundaryFieldRef()[patchID][facei]  = x_out;
			 volfrac.boundaryFieldRef()[patchID][facei]  = a_out;
//			 muBulk.boundaryFieldRef()[patchID][facei]   = 0.125*pp*5e-6;
//                e.boundaryFieldRef()[patchID][facei]     = ee + e_ref;										//              
//                rho.boundaryFieldRef()[patchID][facei]   = 1.0/vv;
                         __transprop_MOD_cpco2(&cp_out, &vv, &vp_out, &x_out, &TT, &pp, &flag);
                         cp.boundaryFieldRef()[patchID][facei] = cp_out;

                         __transprop_MOD_cvco2(&cv_out, &vv, &vp_out, &x_out, &TT, &pp, &flag);
                         cv.boundaryFieldRef()[patchID][facei] = cv_out;

			 __transprop_MOD_co2visco(&vis, &vv, &TT, &x_out, &pp, &flag);
                         mu.boundaryFieldRef()[patchID][facei]  = vis;

			 __transprop_MOD_co2conduc2phase(&lambda_out, &vv, &vp_out, &x_out, &TT, &pp, &flag);
			 lambda_l.boundaryFieldRef()[patchID][facei] = lambda_out;

         		__transprop_MOD_dedrco2(&dedr_out, &vv, &vp_out, &x_out, &TT, &pp, &flag);
         		dedr.boundaryFieldRef()[patchID][facei] = dedr_out;



		}
//         nu[celli] = vis*vv;
//         alpha_lam[celli] = 0.001;
//         Info<< "\n pp" << pp  << endl;
//Info<< "\n TT" << TT  << endl;
	}
//        forAll(mesh.boundary(), patchID)
//	{
//		forAll (mesh.boundary()[patchID],facei)
//                {
//	                 ee = e.boundaryField()[patchID][facei] - e_ref;
//                         vv = 1.0/rho.boundaryField()[patchID][facei];
//Info<< "\n ee bc" << ee  <<" "<<patchID<<" "<<facei <<endl;
//
//			 scalar xxx = center[facei].x(); //center[celli].x();
//                         __interp_table_MOD_co2bllt_equi(&pp,&TT,&cc,&x_out,&a_out,&xxx,&ee,&vv,&flag);


//                nu.boundaryFieldRef()[patchID][facei]  = vis*vv;
//                alpha_lam.boundaryFieldRef()[patchID][facei]  = 0.001;
//               }
//       }
Info<<"Calculating properties ... " <<endl;
//
