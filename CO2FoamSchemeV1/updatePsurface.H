//const volVectorField&     center = mesh.C();
double eepos;
double vvpos;
double eeneg;
double vvneg;
double ppos;
double pneg;
double TTpos,ccpos,x_outpos,a_outpos,xxxpos;
int    flagpos;

        forAll(e_pos,celli)
        {
//         scalar xxx = center[celli].x(); //center[celli].x();
//Info<<" ee= "<<e[celli]<<" rho= "<<rho[celli]<<" i= "<<celli<<" x= "<<xxx <<endl;
         eepos = e_pos[celli] - e_ref;
         vvpos = 1.0/rho_pos[celli];
	 
//
//Info<<" ee= "<<ee<<" vv= "<<vv<<" i= "<<celli<<" x= "<<xxx <<endl;
         __interp_table_MOD_co2bllt_equi(&ppos,&TTpos,&ccpos,&x_outpos,&a_outpos,&xxxpos,&eepos,&vvpos,&flagpos); 
//
//         press[celli] = pp;
//         temp[celli]  = TT;
         p_pos[celli]  = ppos;
	}


        forAll(e_neg,celli)
        {
//         scalar xxx = center[celli].x(); //center[celli].x();
//Info<<" ee= "<<e[celli]<<" rho= "<<rho[celli]<<" i= "<<celli<<" x= "<<xxx <<endl;
         eeneg = e_neg[celli] - e_ref;
         vvneg = 1.0/rho_neg[celli];
	 
//
//Info<<" ee= "<<ee<<" vv= "<<vv<<" i= "<<celli<<" x= "<<xxx <<endl;
         __interp_table_MOD_co2bllt_equi(&pneg,&TTpos,&ccpos,&x_outpos,&a_outpos,&xxxpos,&eeneg,&vvneg,&flagpos); 
//
//         press[celli] = pp;
//         temp[celli]  = TT;
         p_neg[celli]  = pneg;
	}
