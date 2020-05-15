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

\*---------------------------------------------------------------------------*/

#include "SpanWagner.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::SpanWagner<Specie>::SpanWagner
(
    const Specie& sp
)
:
    Specie(sp)
//    Tc_(Tc),
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::SpanWagner<Specie>::SpanWagner
(
    const word& name,
    const SpanWagner& pg
)
:
//    Specie(name, pg)
    Specie(name)
//    omega_(pg.omega_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::SpanWagner <Specie>>
Foam::SpanWagner<Specie>::clone() const
{
    return autoPtr<SpanWagner<Specie>>
    (
        new SpanWagner<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::SpanWagner<Specie>>
Foam::SpanWagner<Specie>::New
(
    Istream& is
)
{
    return autoPtr<SpanWagner<Specie>>(new SpanWagner<Specie>(is));
}


template<class Specie>
inline Foam::autoPtr<Foam::SpanWagner<Specie>>
Foam::SpanWagner<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<SpanWagner<Specie>>
    (
        new SpanWagner<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//template<class Specie>
//inline Foam::scalar Foam::SpanWagner<Specie>::rho
//(
//    scalar p,
//    scalar T
//) const
//{
//    scalar Z = this->Z(p, T);
//    return p/(Z*this->R()*T);
//}


//template<class Specie>
//inline Foam::scalar Foam::SpanWagner<Specie>::h(scalar p, scalar T) const
//{
//    scalar Pr = p/Pc_;
//    scalar Tr = T/Tc_;
//    scalar B = 0.07780*Pr/Tr;
//    scalar kappa = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);
//    scalar alpha = sqr(1 + kappa*(1 - sqrt(Tr)));

//    scalar Z = this->Z(p, T);

//    return 0;
//        RR*Tc_
//       *(
//           Tr*(Z - 1)
//         - 2.078*(1 + kappa)*sqrt(alpha)
//          *log((Z + 2.414*B)/(Z - 0.414*B))
//        );
//}


//template<class Specie>
//inline Foam::scalar Foam::SpanWagner<Specie>::cp(scalar p, scalar T) const
//{
//    scalar Tr = T/Tc_;
//    scalar a = 0.45724*sqr(RR*Tc_)/Pc_;
//    scalar b = 0.07780*RR*Tc_/Pc_;
//    scalar kappa = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);
//    scalar alpha = sqr(1 + kappa*(1 - sqrt(Tr)));
//
//    scalar A = a*alpha*p/sqr(RR*T);
//   scalar B = b*p/(RR*T);

//   scalar Z = this->Z(p, T);

//    scalar ap = kappa*a*(kappa/Tc_ - (1 + kappa)/sqrt(T*Tc_));
//    scalar app = kappa*a*(1 + kappa)/(2*sqrt(pow3(T)*Tc_));

//    scalar M = (sqr(Z) + 2*B*Z - sqr(B))/(Z - B);
//    scalar N = ap*B/(b*RR);

//    const scalar root2 = sqrt(2.0);

//    return 0;
//        app*(T/(2*root2*b))*log((Z + (root2 + 1)*B)/(Z - (root2 - 1)*B))
//      + RR*sqr(M - N)/(sqr(M) - 2*A*(Z + B))
//      - RR;
//}


//template<class Specie>
//inline Foam::scalar Foam::SpanWagner<Specie>::s
//(
//    scalar p,
//    scalar T
//) const
//{
 //   scalar Pr = p/Pc_;
 //   scalar Tr = T/Tc_;
 //   scalar B = 0.07780*Pr/Tr;
 //   scalar kappa = 0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_);

 //   scalar Z = this->Z(p, T);

//    return 0;
 //     - RR*log(p/Pstd)
 //     + RR
 //      *(
 //          log(Z - B)
//         - 2.078*kappa*((1 + kappa)/sqrt(Tr) - kappa)
//          *log((Z + 2.414*B)/(Z - 0.414*B))
 //       );
//}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::SpanWagner<Specie>::operator+=
(
    const SpanWagner<Specie>& pg
)
{
}


template<class Specie>
inline void Foam::SpanWagner<Specie>::operator-=
(
    const SpanWagner<Specie>& pg
)
{

}


template<class Specie>
inline void Foam::SpanWagner<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


//template<class Specie>
//Foam::SpanWagner<Specie> Foam::operator+
//(
//    const SpanWagner<Specie>& pg1,
//    const SpanWagner<Specie>& pg2
//)
//{

//}


//template<class Specie>
//Foam::SpanWagner<Specie> Foam::operator-
//(
//    const SpanWagner<Specie>& pg1,
//    const SpanWagner<Specie>& pg2
//)
//{
//}


//template<class Specie>
//Foam::SpanWagner<Specie> Foam::operator*
//(
//    const scalar s,
//    const SpanWagner<Specie>& pg
//)
//{
//}


//template<class Specie>
//Foam::SpanWagner<Specie> Foam::operator==
//(
//    const SpanWagner<Specie>& pg1,
//    const SpanWagner<Specie>& pg2
//)
//{
//    return 0;
//}


// ************************************************************************* //
