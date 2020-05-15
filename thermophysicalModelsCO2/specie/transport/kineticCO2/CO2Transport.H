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
    Foam::CO2Transport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    CO2TransportI.H
    CO2Transport.C

\*---------------------------------------------------------------------------*/

#ifndef CO2Transport_H
#define CO2Transport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class CO2Transport;

template<class Thermo>
inline CO2Transport<Thermo> operator+
(
    const CO2Transport<Thermo>&,
    const CO2Transport<Thermo>&
);

template<class Thermo>
inline CO2Transport<Thermo> operator-
(
    const CO2Transport<Thermo>&,
    const CO2Transport<Thermo>&
);

template<class Thermo>
inline CO2Transport<Thermo> operator*
(
    const scalar,
    const CO2Transport<Thermo>&
);

template<class Thermo>
inline CO2Transport<Thermo> operator==
(
    const CO2Transport<Thermo>&,
    const CO2Transport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const CO2Transport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class CO2Transport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class CO2Transport
:
    public Thermo
{
    // Private data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;


    // Private Member Functions

        //- Construct from components
        inline CO2Transport
        (
            const Thermo& t,
            const scalar mu,
            const scalar Pr
        );


public:

    // Constructors

        //- Construct as named copy
        inline CO2Transport(const word&, const CO2Transport&);

        //- Construct from Istream
        CO2Transport(Istream&);

        //- Construct from dictionary
        CO2Transport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<CO2Transport> clone() const;

        // Selector from Istream
        inline static autoPtr<CO2Transport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<CO2Transport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "CO2_kinetic<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar rho, const scalar T);

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar rho, const scalar T);

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar rho, const scalar T);

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const CO2Transport&);

        inline void operator+=(const CO2Transport&);

        inline void operator-=(const CO2Transport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend CO2Transport operator+ <Thermo>
        (
            const CO2Transport&,
            const CO2Transport&
        );

        friend CO2Transport operator- <Thermo>
        (
            const CO2Transport&,
            const CO2Transport&
        );

        friend CO2Transport operator* <Thermo>
        (
            const scalar,
            const CO2Transport&
        );

        friend CO2Transport operator== <Thermo>
        (
            const CO2Transport&,
            const CO2Transport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const CO2Transport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CO2TransportI.H"

#ifdef NoRepository
    #include "CO2Transport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
