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
#include "makeThermo.H"

#include "specie.H"
#include "SpanWagner.H"

#include "CO2SW.H"
#include "CO2InternalEnergy.H"
#include "thermo_CO2.H"
//#include "thermo.H"

#include "CO2Transport.H"


#include "CO2heRhoThermo.H"
#include "pureMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

//makeThermo
//(
//    rhoThermo,
//    heRhoThermo,
//    pureMixture,
//    constTransport,
//    sensibleEnthalpy,
//    hConstThermo,
//    perfectGas,
//    specie
//);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// *** CO2 ****
CO2makeThermo
(   
    CO2Thermo,
    CO2heRhoThermo,
    pureMixture,
    CO2Transport,
    CO2InternalEnergy,
    CO2SW,
    SpanWagner,
    specie
);

// *************************************************************************//

//makeThermo
//(
//    rhoThermo,
//    heRhoThermo,
//    pureMixture,
//    constTransport,
//    sensibleInternalEnergy,
//    hConstThermo,
//    perfectGas,
//    specie
//);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
