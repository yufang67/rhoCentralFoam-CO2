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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::CO2SW<EquationOfState>::CO2SW
(
    const EquationOfState& st
//    const scalar cp,
//    const scalar hf
)
:
    EquationOfState(st)
//    Cp_(cp),
//    Hf_(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::CO2SW<EquationOfState>::CO2SW
(
    const word& name,
    const CO2SW& ct
)
:
    EquationOfState(name, ct)
//    Cp_(ct.Cp_),
//    Hf_(ct.Hf_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::CO2SW<EquationOfState>>
Foam::CO2SW<EquationOfState>::clone() const
{
    return autoPtr<CO2SW<EquationOfState>>
    (
        new CO2SW<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::CO2SW<EquationOfState>>
Foam::CO2SW<EquationOfState>::New(Istream& is)
{
    return autoPtr<CO2SW<EquationOfState>>
    (
        new CO2SW<EquationOfState>(is)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::CO2SW<EquationOfState>>
Foam::CO2SW<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<CO2SW<EquationOfState>>
    (
        new CO2SW<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//template<class EquationOfState>
//inline Foam::scalar Foam::CO2SW<EquationOfState>::limit
//(
//    const scalar T
//) const
//{
//    return T;
//}


//template<class EquationOfState>
//inline Foam::scalar Foam::CO2SW<EquationOfState>::cp
//(
//    const scalar p,
//    const scalar T
//) const
//{
//    return Cp_ + EquationOfState::cp(p, T);
//}


//template<class EquationOfState>
//inline Foam::scalar Foam::CO2SW<EquationOfState>::ha
//(
//    const scalar p, const scalar T
//) const
//{
//    return Cp_*T + Hf_ + EquationOfState::h(p, T);
//}


//template<class EquationOfState>
//inline Foam::scalar Foam::CO2SW<EquationOfState>::hs
//(
//    const scalar p, const scalar T
//) const
//{
//    return Cp_*T + EquationOfState::h(p, T);
//}


//template<class EquationOfState>
//inline Foam::scalar Foam::CO2SW<EquationOfState>::hc() const
//{
//    return Hf_;
//}


//template<class EquationOfState>
//inline Foam::scalar Foam::CO2SW<EquationOfState>::s
//(
//    const scalar p, const scalar T
//) const
//{
//    return Cp_*log(T/Tstd) + EquationOfState::s(p, T);
//}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//template<class EquationOfState>
//inline void Foam::CO2SW<EquationOfState>::operator+=
//(
//    const CO2SW<EquationOfState>& ct
//)
//{
//    scalar molr1 = this->nMoles();

//    EquationOfState::operator+=(ct);

//    molr1 /= this->nMoles();
//    scalar molr2 = ct.nMoles()/this->nMoles();

//    Cp_ = molr1*Cp_ + molr2*ct.Cp_;
//    Hf_ = molr1*Hf_ + molr2*ct.Hf_;
//}


//template<class EquationOfState>
//inline void Foam::CO2SW<EquationOfState>::operator-=
//(
//    const CO2SW<EquationOfState>& ct
//)
//{
//    scalar molr1 = this->nMoles();

//    EquationOfState::operator-=(ct);

//    molr1 /= this->nMoles();
//    scalar molr2 = ct.nMoles()/this->nMoles();

//    Cp_ = molr1*Cp_ - molr2*ct.Cp_;
//    Hf_ = molr1*Hf_ - molr2*ct.Hf_;
//}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

//template<class EquationOfState>
//inline Foam::CO2SW<EquationOfState> Foam::operator+
//(
//    const CO2SW<EquationOfState>& ct1,
//    const CO2SW<EquationOfState>& ct2
//)
//{
//    EquationOfState eofs
//    (
//        static_cast<const EquationOfState&>(ct1)
//      + static_cast<const EquationOfState&>(ct2)
//    );

//    return CO2SW<EquationOfState>
//    (
//        eofs,
//        ct1.nMoles()/eofs.nMoles()*ct1.Cp_
//      + ct2.nMoles()/eofs.nMoles()*ct2.Cp_,
//        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
//      + ct2.nMoles()/eofs.nMoles()*ct2.Hf_
//    );
//}

// ************************************************************************* //
