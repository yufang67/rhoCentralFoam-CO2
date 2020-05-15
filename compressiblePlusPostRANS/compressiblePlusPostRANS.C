/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    compressiblePlusPostRANS

Description

    calculates y+ and u+ fields for wall-bounded compressible flows 
    computed with one of the available low-Re RANS (no wall function!)
    turbulence models. More specifically it

    :: 	calculates and outputs y+ (avg., min., max.) based on the 
	velocity gradient at the wall  

    ::	calculates and outputs the wall averaged friction velocity 

    ::  writes fields of y+ and U+ to the corresponding time directory

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "RASModel.H"
#include "wallFvPatch.H" 
//aali
#include "compressibleTurbulenceModel.H"
#include "basicThermo.H"
//aali end
#include "nearWallDist.H"
#include "wallDist.H"
#include "volFields.H"

#include "CO2Thermo.H"
#include "turbulentFluidThermoModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

//        fvMesh::readUpdateState state = mesh.readUpdate();

//        wallDist y(mesh, true);

//        if (timeI == 0 || state != fvMesh::UNCHANGED)
//        {
//            Info<< "Calculating wall distance\n" <<endl;
//            Info<< "Writing wall distance to field " << y.name() << nl << endl;
//            y.write();
//        }

	#include "createFields.H"

//        volScalarField::Boundary& yPlusBf = yPlus.boundaryFieldRef();

        const fvPatchList& patches = mesh.boundary();
        
//        dimensionedScalar uTauAvg("uTauAvg", dimVelocity, 0);

//        scalar nPatch = 0;
                
	Info<< "Summary: " << nl << endl;

        forAll(patches, patchi)
        {
            const fvPatch& currPatch = patches[patchi];

            if (typeid(currPatch) == typeid(wallFvPatch))
            {
//Info<<" dist: "<<d[patchi]<<endl;         
                yPlusTemp.boundaryFieldRef()[patchi] =
                    (d[patchi]*0.5) /( mu.boundaryField()[patchi])
                    *(rho.boundaryField()[patchi])
                    *sqrt
                    (
                       mag(U.boundaryField()[patchi].snGrad())
                      *(muEff.boundaryField()[patchi])
                      /(rho.boundaryField()[patchi])
                    );
                   
                
		const scalarField& YpTemp = yPlusTemp.boundaryField()[patchi];
		const scalarField& dist   = d[patchi];


                Info<< "  y+ for Patch " << patchi
                    << " named " << currPatch.name() << ":" 
                    << " min: " << min(YpTemp) << " max: " << max(YpTemp)
                    << " average: " << average(YpTemp) 
                    << " average distance: " << average(dist) 
                    << " avgUGradWall: " <<  average(mag(U.boundaryField()[patchi].snGrad())) << nl << endl;
            }
       }

        yPlus.boundaryFieldRef() = yPlusTemp.boundaryField();
//        uTauAvg /= nPatch; 
        
//        Info << "  avg. friction velocity uTau is: "
//             << uTauAvg.value() << " (averaged over " << nPatch << " wall(s))" << nl <<endl;

//        yPlus = y.y() * uTauAvg / (mu) *(rho);
//        uPlus = UMean / uTauAvg;
        
//        Info << "Writing yPlus and uPlus to corresponding fields." << nl <<endl;
        yPlus.write();
//        uPlus.write();

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
