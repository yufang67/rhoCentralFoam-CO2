/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "CO2Thermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CO2Thermo, 0);
    defineRunTimeSelectionTable(CO2Thermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CO2Thermo::CO2Thermo(const fvMesh& mesh, const word& phaseName)
:
    fluidThermo(mesh, phaseName),
    rho_
    (
        IOobject
        (
//            phasePropertyName("thermo:rho"),
            phasePropertyName("rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh//,
//        dimDensity
    ),

//    psi_
//    (
//        IOobject
//        (
//            phasePropertyName("thermo:psi"),
//            mesh.time().timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh,
//        dimensionSet(0, -2, 2, 0, 0)
//    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


Foam::CO2Thermo::CO2Thermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    fluidThermo(mesh, dict, phaseName),
    rho_
    (
        IOobject
        (
//            phasePropertyName("thermo:rho"),
            phasePropertyName("rho"),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh//,
//        dimDensity
    ),

//    psi_
//    (
//        IOobject
//        (
//            phasePropertyName("thermo:psi"),
//            mesh.time().timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh,
//        dimensionSet(0, -2, 2, 0, 0)
//    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::CO2Thermo> Foam::CO2Thermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<CO2Thermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CO2Thermo::~CO2Thermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::CO2Thermo::rho() const
{
   return rho_;
}


Foam::tmp<Foam::scalarField> Foam::CO2Thermo::rho(const label patchi) const
{
   return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::CO2Thermo::rho()
{
    return rho_;
}


const Foam::volScalarField& Foam::CO2Thermo::psi() const
{
//    return psi_;
    return mu_;
}


Foam::tmp<Foam::volScalarField> Foam::CO2Thermo::mu() const
{
//    Info<<"CO2mu"<<endl;
    return mu_;
}
Foam::volScalarField&  Foam::CO2Thermo::mu()
{
    Info<<"CO2mu_vary"<<endl;
    return mu_;
}

//Foam::volScalarField&  Foam::CO2Thermo::nu()
Foam::tmp<Foam::volScalarField> Foam::CO2Thermo::nu() const
{
    return mu_/rho_;
}

Foam::tmp<Foam::scalarField> Foam::CO2Thermo::nu(const label patchi) const
{
//    Info<<"CO2mu_BC"<<endl;
    return mu_.boundaryField()[patchi]/rho_.boundaryField()[patchi];
}

Foam::tmp<Foam::scalarField> Foam::CO2Thermo::mu(const label patchi) const
{
//    Info<<"CO2mu_BC"<<endl;
    return mu_.boundaryField()[patchi];
}


// ************************************************************************* //
