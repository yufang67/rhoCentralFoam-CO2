//Info<< "\nUpdating inlet/outlet conditions" << endl;
rho.correctBoundaryConditions();
U.correctBoundaryConditions();
e.correctBoundaryConditions();
p.correctBoundaryConditions();
T.correctBoundaryConditions();
//
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
rhoE.boundaryFieldRef() == rho.boundaryField()*(
                           e.boundaryField() + 0.5*magSqr(U.boundaryField())
                           );
//
//
const volVectorField&     C = mesh.C();
scalar M_ii = 0.0;
scalar res2,res3,res4;
scalar dp_dv_u, dp_du_v;
//
label     adj1_x=0;//,adj2_x=0;
scalar dt = runTime.deltaTValue();
scalar p_aver=0, num=0, super=0, num2=0;
//
//
//
forAll (mesh.boundaryMesh(),ipatch) 
{
       word BCTypePhysical = mesh.boundaryMesh().physicalTypes()[ipatch];
       word BCType         = mesh.boundaryMesh().types()[ipatch];
       word BCName         = mesh.boundaryMesh().names()[ipatch];
       const UList<label> &bfaceCells = mesh.boundaryMesh()[ipatch].faceCells();
//Info<< " " << BCTypePhysical <<" "<< BCType <<" "<< BCName << endl;
//
       if ( BCName == "outlet" ) {
        if (BCTypePhysical =="NSCBC"){
// find max mach (defined in initFields.H)
          forAll ( bfaceCells , iface ) {
                label  bcCell = bfaceCells[iface];
                M_ii = mag(U[bcCell])/c[bcCell];
                if (M_max < M_ii) {
                  M_max = M_ii;
                                  }
                                        }
//
// NSCBC for outlet
          forAll ( bfaceCells , iface ) {
                 label  bcCell = bfaceCells[iface];
                 labelList adjI  = mesh.cellCells()[bcCell];
// find adjacent cell
//                 scalar x_bc  = C[bcCell].x();
                 scalar xmin_adj=100.0;
                 forAll (adjI,adjII){
                        scalar x_adj = C[adjI[adjII]].x();
                        if (x_adj<xmin_adj){
                           adj1_x   = adjI[adjII];
                           xmin_adj = x_adj;
                                           }
                                    }
//Info << adj1_x << endl;
//compute gradient at boundary cells
//                 scalar deltax   =  C[bcCell].x()-C[adj1_x].x();
                 scalar drhodx   = (rho[bcCell] - rho[adj1_x]) / (C[bcCell].x()-C[adj1_x].x()); //gradrho[bcCell] & normal;
                 scalar dpdx     = (p[bcCell]   - p[adj1_x])   / (C[bcCell].x()-C[adj1_x].x()); //gradp[bcCell]   & normal;
                 scalar dudx     = (U[bcCell].x() - U[adj1_x].x()) / (C[bcCell].x()-C[adj1_x].x());
                 scalar dvdx     = (U[bcCell].y() - U[adj1_x].y()) / (C[bcCell].x()-C[adj1_x].x());
                 scalar dwdx     = (U[bcCell].z() - U[adj1_x].z()) / (C[bcCell].x()-C[adj1_x].x());
//Info << C[bcCell].x()-C[adj1_x].x() << endl;
//Info << rho[bcCell] - rho[adj1_x] << endl;
                 scalar p_i      = p[bcCell];
                 scalar rho_i    = rho[bcCell];
                 scalar v_i      = 1.0/rho_i;
                 scalar c_i      = c[bcCell];
                 scalar e_i      = e[bcCell] - e_ref;
                 scalar T_i      = T[bcCell];
                 scalar ux       = U[bcCell].x();
                 scalar uy       = U[bcCell].y();
                 scalar uz       = U[bcCell].z();
                 scalar rhoE_i   = rhoE[bcCell];
                 scalar M_i      = mag(U[bcCell])/c[bcCell];
                 p_aver          = p_i + p_aver;
                 num             = 1.0 + num;
//
//                 scalar cbc      = c.boundaryFieldRef()[ipatch][iface];
//                 scalar mbc      = mag(U.boundaryFieldRef()[ipatch][iface])/cbc;
//                 scalar p_i      = p[adj1_x];
//                 scalar rho_i    = rho[adj1_x];
//                 scalar c_i      = c[adj1_x];
//                 scalar ux       = U[adj1_x].x();
//                 scalar uy       = U[adj1_x].y();
//                 scalar uz       = U[adj1_x].x();
//                 scalar M_i      = mag(U[adj1_x])/c[adj1_x];
                 if (M_max<1.0)
//                 if (M_i<1.0)
//                 if (mbc<1.0)
                 {  
//                 scalar K_i      = relax*0.25*c_i*(1.0-M_max*M_max)/lx;
                 scalar K_i      = relax*0.25*c_i*(1.0-M_i*M_i)/lx;
// characteristic waves
                 scalar L1       = ux * (drhodx - dpdx/(c_i*c_i));
                 scalar L2       = ux * dvdx; 
                 scalar L3       = ux * dwdx;
                 scalar L4       = (ux+c_i) * ( dudx+dpdx/(rho_i*c_i) );
                 scalar L5       = K_i * (p_i - p_out);
// characteristic variables
                 scalar d1       =  L1 + rho_i*(L4+L5)/(2.0*c_i);
                 scalar d2       =  0.5*(L4-L5);
                 scalar d3       =  L2;
                 scalar d4       =  L3;
                 scalar d5       =  rho_i*c_i*0.5*(L4+L5);
// compute derivatieves
//
//Info<<bcCell<<" ei= "<<e_i<<" vi= "<<v_i<<"T= "<<T_i<<"p= "<<p_i<<endl;
                 __derivees_MOD_co2der(&dp_dv_u, &dp_du_v , &e_i, &v_i, &T_i, &p_i,&res2,&res3,&res4);
//Info<<bcCell<<" ei= "<<e_i<<" vi= "<<v_i<<" dpdu= "<<dp_du_v<<"dpdv= "<<dp_dv_u<<endl;


                 scalar dvdr    =  -1.0 / (rho_i*rho_i);
                 scalar dedp_r  =   1.0 / dp_du_v;
                 scalar dedr_p  =  -dp_dv_u*dvdr / dp_du_v; 
// varaibles conservative
//
                 scalar Qc1     = rho_i - dt*d1;
                 scalar Qc2     = rho_i*ux - dt*(ux*d1+rho_i*d2);                 
                 scalar Qc3     = rho_i*uy - dt*(uy*d1+rho_i*d3);                 
                 scalar Qc4     = rho_i*uz - dt*(uz*d1+rho_i*d4);                 
                 scalar Qc5     = rhoE_i - dt*(rhoE_i/rho_i*d1 + 
                                              rho_i*ux*d2 + rho_i*uy*d3 + rho_i*uz*d4 + 
                                              rho_i*(dedp_r*d5 + dedr_p*d1));
//
//
                 rho.boundaryFieldRef()[ipatch][iface]       = Qc1;
                 rhoU.boundaryFieldRef()[ipatch][iface].x()  = Qc2;
                 rhoU.boundaryFieldRef()[ipatch][iface].y()  = Qc3;
                 rhoU.boundaryFieldRef()[ipatch][iface].z()  = Qc4;
                 rhoE.boundaryFieldRef()[ipatch][iface]      = Qc5;
//
//
                 U.boundaryFieldRef()[ipatch][iface].x() = Qc2 / Qc1; 
                 U.boundaryFieldRef()[ipatch][iface].y() = Qc3 / Qc1; 
                 U.boundaryFieldRef()[ipatch][iface].z() = Qc4 / Qc1; 
                 e.boundaryFieldRef()[ipatch][iface]     = Qc5 / Qc1 - 0.5*magSqr(U.boundaryFieldRef()[ipatch][iface]); 
                 p.boundaryFieldRef()[ipatch][iface]     = p_i - dt * d5;
                 }
                 else
                 {
                 rho.boundaryFieldRef()[ipatch][iface]   = rho[bcCell];
                 U.boundaryFieldRef()[ipatch][iface]     = U[bcCell];
                 e.boundaryFieldRef()[ipatch][iface]     = e[bcCell];
                 p.boundaryFieldRef()[ipatch][iface]     = p[bcCell];
                 T.boundaryFieldRef()[ipatch][iface]     = T[bcCell];
//
                 rhoU.boundaryFieldRef()[ipatch][iface]  = rho[bcCell] * U[bcCell];
                 rhoE.boundaryFieldRef()[ipatch][iface]  = rho[bcCell] * (e[bcCell]+0.5*magSqr(U[bcCell]));               
//                   rho[bcCell]   =  rho[adj1_x];
//                   U[bcCell]     =  U[adj1_x];
//                   e[bcCell]     =  e[adj1_x];
//                 Info << "supersonic outlet "<<"Mach= " << M_i << endl; 
                 }//endif sub/supersonic
          }//end outlet cells
Info << "outlet NSCBC "<< p_aver<<" "<<num<<" "<<M_max<<endl;
//          rhoU.boundaryFieldRef()[ipatch] = rho.boundaryField()[ipatch]*U.boundaryField()[ipatch];
//          rhoE.boundaryFieldRef()[ipatch] = rho.boundaryField()[ipatch]*
//                                           ( e.boundaryField()[ipatch]+ 
//                                             0.5*magSqr(U.boundaryField()[ipatch]) 
//                                           );
        }
        else if (BCTypePhysical == "zeroExtrapolation"){
          Info << " outlet zero extrapolation not complete "<< endl;
          forAll ( bfaceCells , iface ) {
//Info<< bfaceCells << " iface= "<< iface << endl;             
                 label  bcCell = bfaceCells[iface];
                 labelList adjI  = mesh.cellCells()[bcCell];
// find adjacent cell
                 scalar xmin_adj=100.0;
                 forAll (adjI,adjII){
                        scalar x_adj = C[adjI[adjII]].x();
                        if (x_adj<xmin_adj){
                           adj1_x   = adjI[adjII];
                           xmin_adj = x_adj;
                        }
                 }
//
                 rho[bcCell] = rho[adj1_x];
                 U[bcCell]   = U[adj1_x];
                 p[bcCell]   = p[adj1_x];
                 e[bcCell]   = e[adj1_x];
                 T[bcCell]   = T[adj1_x];

          }//end outlet loop
        }//end type loop
       }//end outlet loop
//
//
       if ( BCName == "inlet1" or BCName == "inlet") 
       {
        if (BCTypePhysical == "NSCBC") 
	{
          forAll ( bfaceCells , iface ) 
	  {
//Info<< bfaceCells << " iface= "<< iface << endl;             
                 label bcCell = bfaceCells[iface];
//
                 rho.boundaryFieldRef()[ipatch][iface]   = rho_inlet1;
                 e.boundaryFieldRef()[ipatch][iface]     = e_inlet1 + e_ref;
                 U.boundaryFieldRef()[ipatch][iface]     = U[bcCell];
                 T.boundaryFieldRef()[ipatch][iface]     = T_inlet1;
                 p.boundaryFieldRef()[ipatch][iface]     = p_inlet1;
//
                 rhoU.boundaryFieldRef()[ipatch][iface]  = rho_inlet1 * U[bcCell];
                 rhoE.boundaryFieldRef()[ipatch][iface]  = rho_inlet1 * ( e_inlet1+e_ref+0.5*magSqr(U[bcCell]) );

                 num2 = 1 + num2;

          }//end inlet cells
          Info << "inlet1  "<<num2<<endl;
//	  num2 = 0;
        }
       }	
       else if (BCName == "inlet2")
       {
        if (BCTypePhysical == "NSCBC") 
	{
          forAll ( bfaceCells , iface ) 
	  {
//Info<< bfaceCells << " iface= "<< iface << endl;             
                 label bcCell = bfaceCells[iface];
//
                 rho.boundaryFieldRef()[ipatch][iface]   = rho_inlet2;
                 e.boundaryFieldRef()[ipatch][iface]     = e_inlet2 + e_ref;
                 U.boundaryFieldRef()[ipatch][iface]     = U[bcCell];
                 T.boundaryFieldRef()[ipatch][iface]     = T_inlet2;
                 p.boundaryFieldRef()[ipatch][iface]     = p_inlet2;
//
                 rhoU.boundaryFieldRef()[ipatch][iface]  = rho_inlet2 * U[bcCell];
                 rhoE.boundaryFieldRef()[ipatch][iface]  = rho_inlet2 * ( e_inlet2+e_ref+0.5*magSqr(U[bcCell]) );

                 num2 = 1 + num2;

          }//end inlet cells
          Info << "inlet2  "<<num2<<endl;

        }
       }
}//end patch loop
