/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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

#include "SpanWagner.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::SpanWagner<Specie>::SpanWagner(Istream& is)
:
    Specie(is)
//    Tc_(readScalar(is)),
{
    is.check("SpanWagner<Specie>::SpanWagner(Istream& is)");
}


template<class Specie>
Foam::SpanWagner<Specie>::SpanWagner
(
    const dictionary& dict
)
:
    Specie(dict)
//    Tc_(readScalar(dict.subDict("equationOfState").lookup("Tc"))),
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::SpanWagner<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const SpanWagner<Specie>& pg
)
{
    os  << static_cast<const Specie&>(pg);
//        << token::SPACE << pg.omega_;

    os.check
    (
        "Ostream& operator<<(Ostream& os, const SpanWagner<Specie>& st)"
    );
    return os;
}


// ************************************************************************* //
