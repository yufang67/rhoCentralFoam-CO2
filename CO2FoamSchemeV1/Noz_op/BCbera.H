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
scalar p_out = 6.0e6;
scalar lx    = 0.04;//100e-3;//80e-3;
scalar relax = 10;
scalar p_in  = 10.0e6;
scalar T_in  = 311.33;//313.15;//310.15;
scalar rho_in= 663.64;//621.475;
scalar e_in  = -218360.5;//-181535.542;//-211417.86;
//vectorField gradrho  = fvc::grad(rho);
//vectorField gradp    = fvc::grad(p);
//volTensorField gradU = fvc::grad(U);
//const surfaceVectorField& Cf  = mesh.Cf();
const volVectorField&     C = mesh.C();
scalar es,out_2,residue,in_2,eguess,vbc;
int Niter,MODE=2; //Newton:3, Brent:2
label     adj1_x=0;//,adj2_x=0;
scalar dt = runTime.deltaTValue();
scalar p_aver=0, num=0, super=0, num2=0;
//const surfaceVectorField& deltas = mesh.delta();

//RectangularMatrix <scalar> L(5,1);
//scalarSquareMatrix X(5, 0.0);
//scalarSquareMatrix S(5, 0.0);
//RectangularMatrix <scalar> dQdx(5,1);
//L(1,1) = 1;
//
forAll (mesh.boundaryMesh(),ipatch) {
       word BCTypePhysical = mesh.boundaryMesh().physicalTypes()[ipatch];
       word BCType         = mesh.boundaryMesh().types()[ipatch];
       word BCName         = mesh.boundaryMesh().names()[ipatch];
       const UList<label> &bfaceCells = mesh.boundaryMesh()[ipatch].faceCells();
//Info<< " " << BCTypePhysical <<" "<< BCType <<" "<< BCName << endl;
//
//       scalarField dpdx   = p.boundaryField()[ipatch].snGrad(); 
//       vectorField dUdx   = U.boundaryField()[ipatch].snGrad(); 
//       vectorField drhodx = gradrho * normal;
//Info<< " "<< normal << endl;             
       if ( BCName == "outlet" ) {
        if (BCTypePhysical =="NSCBC"){
          forAll ( bfaceCells , iface ) {
//Info<< bfaceCells << " iface= "<< iface << endl;             
                 label  bcCell = bfaceCells[iface];
//                 vector normal   = mesh.Sf().boundaryField()[ipatch][iface] / 
//                                   mesh.magSf().boundaryField()[ipatch][iface];
//                 scalar deltax      =  cmptMag(C[bcCell].x() +  Cf[bcCell].x()) * 2.0;
//                 scalar deltax      =  cmptMag(C[bcCell].x() +  Cf[bcCell].x());
//
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
                 scalar c_i      = c[bcCell];
                 scalar ux       = U[bcCell].x();
                 scalar uy       = U[bcCell].y();
                 scalar uz       = U[bcCell].x();
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
//
                 if (M_i<1.0)
//                 if (mbc<1.0)
                 {  
                 scalar K_i      = relax*0.25*c_i*(1.0-M_i*M_i)/lx;
//Info << U[bcCell].x() << endl;
//characteristic waves
                 scalar L1       = ux * (drhodx - dpdx/(c_i*c_i));
                 scalar L2       = ux * dvdx; 
                 scalar L3       = ux * dwdx;
                 scalar L4       = (ux+c_i) * ( dudx+dpdx/(rho_i*c_i) );
                 scalar L5       = K_i * (p_i - p_out);
//characteristic variables
                 scalar d1       =  L1 + rho_i*(L4+L5)/(2.0*c_i);
                 scalar d2       =  0.5*(L4-L5);
                 scalar d3       =  L2;
                 scalar d4       =  L3;
                 scalar d5       =  rho_i*c_i*0.5*(L4+L5);
//
//Info << L1 << " "<< rho_i*(L4+L5)/(2.0*c_i)<<endl;
//Info << ux << " "<< drhodx<<endl;
//primitives variables
                 scalar Qp1         =  rho_i - dt*d1;
                 scalar Qp2         =  ux    - dt*d2;
                 scalar Qp3         =  uy    - dt*d3;
                 scalar Qp4         =  uz    - dt*d4;
                 scalar Qp5         =  p_i   - dt*d5;
//
//primitive variables at outlet interface (average of adjacent cells)
//                 scalar Qs1         = 0.5*(Qp1+rho_i);         // |rho
//                 scalar Qs2         = 0.5*(Qp2+U[bcCell].x()); // |u
//                 scalar Qs3         = 0.5*(Qp3+U[bcCell].y()); // |v
//                 scalar Qs4         = 0.5*(Qp4+U[bcCell].z()); // |w
//                 scalar Qs5         = 0.5*(Qp5+p_i);           // |p
//Info <<Qs1 <<Qs2 <<Qs3 << Qs4 <<Qs5 <<endl; 
//
//primitive variables to conservative variables
                 rho.boundaryFieldRef()[ipatch][iface]   = Qp1;
                 U.boundaryFieldRef()[ipatch][iface].x() = Qp2; 
                 U.boundaryFieldRef()[ipatch][iface].y() = Qp3; 
                 U.boundaryFieldRef()[ipatch][iface].z() = Qp4; 
//
//                 rho[bcCell]        = Qp1;
//                 U[bcCell].x()      = Qp2;
//                 U[bcCell].y()      = Qp3;
//                 U[bcCell].z()      = Qp4;
//                 p[bcCell]          = Qp5;
                 
//
//                 eguess             = e[bcCell]-e_ref;
                 eguess             = e[adj1_x]-e_ref;
                 scalar lower       = eguess - 2e4;//25e3; 
                 scalar upper       = eguess + 2e4;//20e3;
                 vbc                = 1.0/Qp1;
//
//Info << " eguess=" << eguess<< " lower=" <<lower<< " upper="<<upper<<" rho="<<Qp1<<" v="<<vbc<<endl;
//Info << " p_inner="<<p_i<<" p_bc="<< Qp5 <<" x.y="<<C[adj1_x].x()<<" "<<C[adj1_x].y() <<endl;
// 
                 __solver_eos_MOD_brentroots2(&MODE,&es,&out_2,&residue,&Niter,&Qp5,
                                              &lower,&upper,&vbc,&in_2);
//                 __solver_eos_MOD_eos_1d(&MODE,&es, &out_2,&residue,&Niter,&flag,&Qp5, &eguess,&vbc,&in_2);
//Info << "e_correct " << es <<" residu= "<<residue <<endl;
                 
                 e.boundaryFieldRef()[ipatch][iface]     = es + e_ref; 
//                 e[bcCell]                               = es + e_ref;
//                 rhoU[bcCell]                            = rho[bcCell] * U[bcCell];
//                 rhoE[bcCell]                            = rho[bcCell] * ( e[bcCell] + 0.5*magSqr(U[bcCell]) );
//
//Info << dUdx.y() << " " << L5 << endl;
//Info<< gradU[bcCell] << "  "<< p[bcCell] << nx << endl;
//Info<< drhodx << "  "<< dpdx << " "<< dUdx << endl;
//Info << U.boundaryFieldRef()[ipatch][iface].x()<<endl;
                 p.boundaryFieldRef()[ipatch][iface]     = Qp5;
                 rhoU.boundaryFieldRef()[ipatch][iface]  = Qp1 * U.boundaryField()[ipatch][iface];
                 rhoE.boundaryFieldRef()[ipatch][iface]  = Qp1 * (e.boundaryField()[ipatch][iface]
                                                               +0.5*magSqr(U.boundaryField()[ipatch][iface]));
                 }
                 else
                 {
                 super=1;
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
	  Info << "outlet NSCBC "<< p_aver<<" "<<num<<" "<<super<<endl;
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
       if ( BCName == "inlet" ) {
        if (BCTypePhysical == "NSCBC") {
          forAll ( bfaceCells , iface ) {
//Info<< bfaceCells << " iface= "<< iface << endl;             
                 label bcCell = bfaceCells[iface];
//
//                 labelList  adjI  = mesh.cellCells()[bcCell];
//Info << adjI << endl;
// find adjacent cell
//                 scalar x_bc  = C[bcCell].x();
//                 scalar xmax_adj = -100.0;
//                 forAll (adjI,i){
//                       scalar x_adj = C[adjI[i]].x();
//Info << x_adj <<" " <<x_bc << " "<< C[1].x()<<endl;
//                       if (x_adj>xmax_adj){
//                          adj2_x   = adjI[i];
//                          xmax_adj = x_adj; 
//                       }
//                 }
//Info << adj2_x << endl;
//
//                 U[bcCell].x()      = U[adj2_x].x();
//                 U[bcCell].y()      = U[adj2_x].y();
//                 U[bcCell].z()      = U[adj2_x].z();
//                 p[bcCell]          = p_in;//p[adj2_x];//p_in; 
//                 e[bcCell]          = e_in + e_ref;//e[adj2_x];//e_in + e_ref;
//                 rho[bcCell]        = rho_in; //rho[adj2_x];//rho_in;
//
//                 rhoU[bcCell]       = rho[bcCell] * U[bcCell];
//                 rhoE[bcCell]       = rho[bcCell] * ( e[bcCell] + 0.5*magSqr(U[bcCell]) );
//
                 rho.boundaryFieldRef()[ipatch][iface]   = rho_in;
                 e.boundaryFieldRef()[ipatch][iface]     = e_in + e_ref;
                 U.boundaryFieldRef()[ipatch][iface]     = U[bcCell];
                 T.boundaryFieldRef()[ipatch][iface]     = T_in;
                 p.boundaryFieldRef()[ipatch][iface]     = p_in;
//
                 rhoU.boundaryFieldRef()[ipatch][iface]  = rho_in * U[bcCell];
                 rhoE.boundaryFieldRef()[ipatch][iface]  = rho_in * ( e_in+e_ref+0.5*magSqr(U[bcCell]) );
		  
		 num2 = 1 + num2;
                 

          }//end inlet cells
          Info << "inlet NSCBC  " <<num2<< endl;
//
//          rhoU.boundaryFieldRef()[ipatch]         = rho.boundaryField()[ipatch]*U.boundaryField()[ipatch];
//          rhoE.boundaryFieldRef()[ipatch]         = rho.boundaryField()[ipatch]*
//                                                    ( e.boundaryField()[ipatch]+
//                                                      0.5*magSqr(U.boundaryField()[ipatch])
//                                                    );
        }
        else if (BCTypePhysical == "zeroExtrapolation") {
          Info << " inlet zero extrapolation not complete" <<endl;
          forAll ( bfaceCells , iface ) {
//                label bcCell = bfaceCells[iface];
//                labelList  adjI  = mesh.cellCells()[bcCell];
//Info << adjI << endl;
// find adjacent cell
//                scalar x_bc  = C[bcCell].x();
//                scalar xmax_adj = -100.0;
//                forAll (adjI,i){
//                   scalar x_adj = C[adjI[i]].x();
//Info << x_adj <<" " <<x_bc << " "<< C[1].x()<<endl;
//                   if (x_adj>xmax_adj){
//                      adj2_x   = adjI[i];
//                      xmax_adj = x_adj; 
//                   }
//                }
//
//                 U[bcCell].x()      = U[adj2_x].x();
//                 U[bcCell].y()      = U[adj2_x].y();
//                 U[bcCell].z()      = U[adj2_x].z();
//                 p[bcCell]          = p[adj2_x];//p_in; 
//                 e[bcCell]          = e[adj2_x];//e_in + e_ref;
//                 rho[bcCell]        = rho[adj2_x];//rho_in;
          }//end inlet cells
        }//end type loop
       } //end inlet loop
}//end patch loop
//Info<< "-->inlet/outlet conditions finish\n" << endl;
//
