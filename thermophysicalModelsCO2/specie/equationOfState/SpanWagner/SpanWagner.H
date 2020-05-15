/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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
    Foam::SpanWagner

Description
    SpanWagner gas equation of state.

SourceFiles
    SpanWagnerI.H
    SpanWagner.C

\*---------------------------------------------------------------------------*/

#ifndef SpanWagner_H
#define SpanWagner_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class SpanWagner;

template<class Specie>
inline SpanWagner<Specie> operator+
(
    const SpanWagner<Specie>&,
    const SpanWagner<Specie>&
);

template<class Specie>
inline SpanWagner<Specie> operator-
(
    const SpanWagner<Specie>&,
    const SpanWagner<Specie>&
);

template<class Specie>
inline SpanWagner<Specie> operator*
(
    const scalar,
    const SpanWagner<Specie>&
);

template<class Specie>
inline SpanWagner<Specie> operator==
(
    const SpanWagner<Specie>&,
    const SpanWagner<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const SpanWagner<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class SpanWagner Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class SpanWagner
:
    public Specie
{
    // Private data
        //- Critical Temperature [K]
//        scalar Tc_;

public:

    // Constructors

        //- Construct from components
        inline SpanWagner( const Specie& sp );

        //- Construct from Istream
        SpanWagner(Istream&);

        //- Construct from dictionary
        SpanWagner(const dictionary& dict);

        //- Construct as named copy
        inline SpanWagner(const word& name, const SpanWagner&);

        //- Construct and return a clone
        inline autoPtr<SpanWagner> clone() const;

        // Selector from Istream
        inline static autoPtr<SpanWagner> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<SpanWagner> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "SpanWagner<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
//            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kmol]
//            inline scalar h(const scalar p, const scalar T) const;

            //- Return cp departure [J/(kmol K]
//            inline scalar cp(scalar p, scalar T) const;

            //- Return entropy [J/(kmol K)]
//            inline scalar s(const scalar p, const scalar T) const;

// IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator+=(const SpanWagner&);
        inline void operator-=(const SpanWagner&);

        inline void operator*=(const scalar);


    // Friend operators

//        friend SpanWagner operator+ <Specie>
//        (
//            const SpanWagner&,
//            const SpanWagner&
//        );

//        friend SpanWagner operator- <Specie>
//        (
//            const SpanWagner&,
//            const SpanWagner&
//        );

//        friend SpanWagner operator* <Specie>
//        (
//            const scalar s,
//            const SpanWagner&
//        );

//        friend SpanWagner operator== <Specie>
//        (
//            const SpanWagner&,
//            const SpanWagner&
//        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const SpanWagner&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SpanWagnerI.H"

#ifdef NoRepository
    #include "SpanWagner.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
