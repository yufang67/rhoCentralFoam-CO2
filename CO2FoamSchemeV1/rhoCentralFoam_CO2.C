/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rhoCentralFoam_CO2

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CO2Thermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
//
#include "CO2table.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"
    #include "createFieldRefs.H"

    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"
//**********************************OP BC************************************//
// nozzle    
//    #include "opCondition.H"
//    #include "BR_opCondition.H"
//    #include "BR85_opCondition.H"
//    #include "JW_opCondition.H"
//    #include "YF_opCondition.H"
//    #include "NK_opCondition.H"
//    #include "NK3_opCondition.H"
//    #include "NK3b_opCondition.H"
//    #include "NK3c_opCondition.H"

// ejector
//    #include "YF10.H"
//    #include "YF10_43.H"
//    #include "YF10_41.H"
    #include "YF_op/YF10_418.H"
//    #include "YF10_422.H"
//    #include "YF10_426.H"
//    #include "YF_op/YF10_434.H"
//    #include "YF10_431.H"
//    #include "YF10_432.H"
//    #include "YF10_38.H"
//    #include "YF10_4418.H"
//    #include "YF10_4457.H"
//    #include "YF10_45.H"
//    #include "ej_opCondition.H"
//    #include "ej_opCondition_sec.H"
//    #include "ej_opCondition_sec_prim.H"
//    #include "ej_ON_opCondition.H"
//    #include "ej_ON2_opCondition.H"
//    #include "ej_78MPa55MPa.H"
//    #include "ej_9MPaON.H"
//    #include "ej_15MPa38.H"
//    #include "ej_15MPaON.H"
//    #include "ej_9MPaONsec4.H"
//    #include "ej_10MPaON.H"
//    #include "ej_10MPa37.H"
//    #include "ej_10MPa3MPa.H"
//    #include "ej_10MPa3MPa33.H"
//    #include "ej_10MPa8MPa.H"
//    #include "ej_10MPa7MPa.H"
//    #include "ej_10MPa64MPa.H"
//    #include "ej_10MPa67MPa.H"
//    #include "ej_10MPa68MPa.H"
//    #include "ej_10MPa69MPa.H"
//
//*******************initial solution***************************//
//
    #include "init_op/initFields_ej.H"
//    #include "initFields_ej3d.H"
//    #include "initFields_out2phase.H"
//    #include "initFields.H"
//
//************************************************************$**//
//
//
    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
	#include "propCalc.H"

        // --- Directed interpolation of primitive fields onto faces

//        volVectorField gradT("gradT", Foam::fvc::grad(T));


        surfaceScalarField rho_pos(interpolate(rho, pos, rho.name() ));
        surfaceScalarField rho_neg(interpolate(rho, neg, rho.name() ));
//Info <<"rho_neg= "<<rho_neg<<endl;
        surfaceScalarField e_pos(interpolate(e, pos, e.name()));
        surfaceScalarField e_neg(interpolate(e, neg, e.name()));

        surfaceVectorField U_pos(interpolate(U, pos, U.name()));
        surfaceVectorField U_neg(interpolate(U, neg, U.name()));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

// interpolation for conservative variables at interfaces ???
//        surfaceScalarField rhoE_pos(interpolate(rhoE, pos, e.name()));
//        surfaceScalarField rhoE_neg(interpolate(rhoE, neg, e.name()));
//        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
//        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);
//        surfaceScalarField e_pos("e_pos", rhoE_pos/rho_pos - 0.5*magSqr(U_pos));
//        surfaceScalarField e_neg("e_neg", rhoE_neg/rho_neg - 0.5*magSqr(U_neg));

//        Info<< "\n rho_pos" << rho_pos  << endl;
//        Info<< "\n e_pos" << e_pos  << endl;
//        Info<< "\n rho" << rho  << endl;
//        Info<< "\n rhoU" << rhoU  << endl;
//        Info<< "\n rhoE" << rhoE  << endl;

//
        surfaceScalarField p_pos(interpolate(p, pos, rho.name()));
        surfaceScalarField p_neg(interpolate(p, neg, rho.name()));
        #include "updatePsurface.H"


//Info<< "\n c " << c  << endl;
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, rho.name())*mesh.magSf()
        );
        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, rho.name())*mesh.magSf()
        );

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());
//Info<< "\n p_pos " << p_pos  << endl;
//Info<< "\n cSf_pos " << cSf_pos  << endl;


        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );
//Info<< "\n ap " << ap  << endl;
//Info<< "\n am " << am  << endl;

        surfaceScalarField a_pos("a_pos", ap/(ap - am));
//Info<< "\n a_pos= " << a_pos << endl;

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
//Info<< "amaxSf= " << amaxSf << endl;
        surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);
//Info<< "a_neg= " << a_neg << endl;
        phiv_pos *= a_pos;
        phiv_neg *= a_neg;
//
        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);
//Info<< "aphiv_pos= " << aphiv_pos << endl;
        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

//Info<< "\n amaxSf " << amaxSf  << endl;
        #include "centralCourantNo.H"
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "setDeltaT.H"
        }
//Info << "timestep "<< runTime.deltaTValue()<<endl;
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
        );

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );

        
// compute heat flux

        volScalarField muEff("muEff", turbulence->muEff());


	lambda_t = cp*(muEff-mu)/0.85;


        lambda_eff = lambda_l + lambda_t;
      
	volScalarField alphaEff("alphaEff", lambda_eff/cv);
	volScalarField beta("beta", alphaEff*dedr);


//        volVectorField gradT("q1", Foam::fvc::grad(alphaEff)*Foam::fvc::grad(e));
//        volVectorField gradT("q2", Foam::fvc::div (Foam::fvc::grad(rho))*beta);
//        volVectorField gradT("q3", Foam::fvc::grad(rho)*Foam::fvc::grad(dedr));
        
//        volScalarField q1("q1", fvc::grad(alphaEff)&fvc::grad(e));
//        volScalarField q2("q2", fvc::div (fvc::grad(rho))*beta);
//        volScalarField q3("q3", fvc::grad(rho)&fvc::grad(beta));

//       volScalarField qsource
//        (
//	 	"qsource",
//		q2
//		fvc::laplacian(beta,rho)
//		 q1+q2+q3
//            -lambda_eff*fvc::grad(T)
//            -lambda_eff*gradT_interp*mesh.Sf()
//              -lambda_eff*
//      
//        );
     
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));
//        volTensorField tauBulk("tauBulk", muBulk*(Foam::fvc::grad(U) - dev(Foam::fvc::grad(U)))*3.0);
//
//###################################
        runTime.write();
        if(runTime.outputTime())
          {
           e.write();
           muEff.write();
	   alphaEff.write();
	   lambda_eff.write();
          }
//###################################

//************************************************
        // --- Solve density(rho)
        solve(fvm::ddt(rho) + fvc::div(phi));
        
//        rho.correctBoundaryConditions();

        // --- Solve momentum(rhoU)
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));


//        U.ref() = rhoU() / rho();
        U       = rhoU   / rho;  
//        U.correctBoundaryConditions();
//        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

        if (!inviscid)
        {
//            Info<<"visco dissip"<<endl;
            U.correctBoundaryConditions();
            solve  //solve U for diffusion
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );
            rhoU = rho*U;
        }

//Info<<"mu_eff "<<muEff<<endl;
//Info<<"tauMC "<<tauMC<<endl;
        // --- Solve energy (rhoE)
        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
                fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
              + fvc::dotInterpolate(mesh.Sf(), tauMC)
            )
          & (a_pos*U_pos + a_neg*U_neg)
        );
//Info <<"sigmaDotU "<<sigmaDotU<<endl;
//
	surfaceScalarField betaDotdRho
        (
	  "betaDotdRho",
          fvc::interpolate(beta)*fvc::snGrad(rho)*mesh.magSf()
        );
        
        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
//	  + fvc::div(betaDotdRho) 
        );

        e = rhoE/rho - 0.5*magSqr(U);
//Info<< rhoE.boundaryFieldRef()<< endl;

        if (!inviscid)
        {
//              volScalarField alphaEff("alphaEff", turbulence->nu()/Pr)
            e.correctBoundaryConditions();
            solve //solve e
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvm::laplacian(alphaEff, e)
	      + fvc::div(betaDotdRho) 
            );
//Info<<" alphaEff= "<<alphaEff<<endl;
//" alpha= "<<alpha<<" alphat= "<<alphat<<" alphah= "<<alphah<<endl;
            rhoE = rho*(e + 0.5*magSqr(U));
        }
//        p.correctBoundaryConditions();
//        rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
//*****************************
        turbulence->correct();
//******************************
//        #include "BCbera.H"
//        #include "BCupdate.H"
        #include "BCupdate2.H"
//******************************

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
