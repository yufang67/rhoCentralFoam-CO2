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
    Foam::CO2SW

Description
    Constant properties thermodynamics package
    templated into the EquationOfState.

SourceFiles
    CO2SWI.H
    CO2SW.C

\*---------------------------------------------------------------------------*/

#ifndef CO2SW_H
#define CO2SW_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class CO2SW;

template<class EquationOfState>
inline CO2SW<EquationOfState> operator+
(
    const CO2SW<EquationOfState>&,
    const CO2SW<EquationOfState>&
);

template<class EquationOfState>
inline CO2SW<EquationOfState> operator-
(
    const CO2SW<EquationOfState>&,
    const CO2SW<EquationOfState>&
);

template<class EquationOfState>
inline CO2SW<EquationOfState> operator*
(
    const scalar,
    const CO2SW<EquationOfState>&
);

template<class EquationOfState>
inline CO2SW<EquationOfState> operator==
(
    const CO2SW<EquationOfState>&,
    const CO2SW<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const CO2SW<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class CO2SW Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class CO2SW
:
    public EquationOfState
{
    // Private data

//        scalar Cp_;
//        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline CO2SW
        (
            const EquationOfState& st
//            const scalar cp,
//            const scalar hf
        );


public:

    // Constructors

        //- Construct from Istream
        CO2SW(Istream&);

        //- Construct from dictionary
        CO2SW(const dictionary& dict);

        //- Construct as named copy
        inline CO2SW(const word&, const CO2SW&);

        //- Construct and return a clone
        inline autoPtr<CO2SW> clone() const;

        //- Selector from Istream
        inline static autoPtr<CO2SW> New(Istream& is);

        //- Selector from dictionary
        inline static autoPtr<CO2SW> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "CO2SW<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
//        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kmol K)]
//            inline scalar cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kmol]
//            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kmol]
//            inline scalar hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kmol]
//            inline scalar hc() const;

            //- Entropy [J/(kmol K)]
//            inline scalar s(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const CO2SW&);
        inline void operator-=(const CO2SW&);


    // Friend operators

//        friend CO2SW operator+ <EquationOfState>
//        (
//            const CO2SW&,
//            const CO2SW&
//        );

//        friend CO2SW operator- <EquationOfState>
//        (
//            const CO2SW&,
//            const CO2SW&
//        );

//        friend CO2SW operator* <EquationOfState>
//        (
//            const scalar,
//            const CO2SW&
//        );

//        friend CO2SW operator== <EquationOfState>
//        (
//            const CO2SW&,
//            const CO2SW&
//        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const CO2SW&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CO2SWI.H"

#ifdef NoRepository
    #include "CO2SW.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
