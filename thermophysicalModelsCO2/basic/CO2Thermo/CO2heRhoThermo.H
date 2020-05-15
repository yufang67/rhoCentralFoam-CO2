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

Class
    Foam::CO2heRhoThermo

Description
    Energy for a mixture based on density

SourceFiles
    CO2heRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef CO2heRhoThermo_H
#define CO2heRhoThermo_H

#include "CO2Thermo.H"
#include "heThermo_CO2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class CO2heRhoThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class CO2heRhoThermo
:
    public heThermo_CO2<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        CO2heRhoThermo(const CO2heRhoThermo<BasicPsiThermo, MixtureType>&);


public:

    //- Runtime type information
    TypeName("CO2heRhoThermo");


    // Constructors

        //- Construct from mesh and phase name
        CO2heRhoThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Destructor
    virtual ~CO2heRhoThermo();


    // Member functions

        //- Update properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CO2heRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
