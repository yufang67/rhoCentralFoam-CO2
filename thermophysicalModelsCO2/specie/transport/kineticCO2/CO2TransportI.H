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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::CO2Transport<Thermo>::CO2Transport
(
    const Thermo& t,
    const scalar mu,
    const scalar Pr
)
:
    Thermo(t),
    mu_(mu),
    rPr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::CO2Transport<Thermo>::CO2Transport
(
    const word& name,
    const CO2Transport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    rPr_(ct.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::CO2Transport<Thermo>>
Foam::CO2Transport<Thermo>::clone() const
{
    return autoPtr<CO2Transport<Thermo>>
    (
        new CO2Transport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::CO2Transport<Thermo>>
Foam::CO2Transport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<CO2Transport<Thermo>>
    (
        new CO2Transport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::CO2Transport<Thermo>>
Foam::CO2Transport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<CO2Transport<Thermo>>
    (
        new CO2Transport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::CO2Transport<Thermo>::mu
(
    const scalar rho,
    const scalar T
)
{
    Info<<"CO2transport_mu"<<endl;
    return mu_;
}


template<class Thermo>
inline Foam::scalar Foam::CO2Transport<Thermo>::kappa
(
    const scalar p,
    const scalar T
)
{
//    return this->Cp(p, T)*mu(p, T)*rPr_;
    return 1.0;
}


template<class Thermo>
inline Foam::scalar Foam::CO2Transport<Thermo>::alphah
(
    const scalar rho,
    const scalar T
) 
{
    return mu_*rPr_;
//    return 1.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::CO2Transport<Thermo>::operator=
(
    const CO2Transport<Thermo>& ct
)
{
    Thermo::operator=(ct);

    mu_ = ct.mu_;
    rPr_ = ct.rPr_;
}


template<class Thermo>
inline void Foam::CO2Transport<Thermo>::operator+=
(
    const CO2Transport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

        mu_ = molr1*mu_ + molr2*st.mu_;
        rPr_ = 1.0/(molr1/rPr_ + molr2/st.rPr_);
    }
}


template<class Thermo>
inline void Foam::CO2Transport<Thermo>::operator-=
(
    const CO2Transport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

        mu_ = molr1*mu_ - molr2*st.mu_;
        rPr_ = 1.0/(molr1/rPr_ - molr2/st.rPr_);
    }
}


template<class Thermo>
inline void Foam::CO2Transport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::CO2Transport<Thermo> Foam::operator+
(
    const CO2Transport<Thermo>& ct1,
    const CO2Transport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(ct1.nMoles()) + mag(ct2.nMoles()) < SMALL)
    {
        return CO2Transport<Thermo>
        (
            t,
            0,
            ct1.rPr_
        );
    }
    else
    {
        scalar molr1 = ct1.nMoles()/t.nMoles();
        scalar molr2 = ct2.nMoles()/t.nMoles();

        return CO2Transport<Thermo>
        (
            t,
            molr1*ct1.mu_ + molr2*ct2.mu_,
            1.0/(molr1/ct1.rPr_ + molr2/ct2.rPr_)
        );
    }
}


template<class Thermo>
inline Foam::CO2Transport<Thermo> Foam::operator-
(
    const CO2Transport<Thermo>& ct1,
    const CO2Transport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) - static_cast<const Thermo&>(ct2)
    );

    if (mag(ct1.nMoles()) + mag(ct2.nMoles()) < SMALL)
    {
        return CO2Transport<Thermo>
        (
            t,
            0,
            ct1.rPr_
        );
    }
    else
    {
        scalar molr1 = ct1.nMoles()/t.nMoles();
        scalar molr2 = ct2.nMoles()/t.nMoles();

        return CO2Transport<Thermo>
        (
            t,
            molr1*ct1.mu_ - molr2*ct2.mu_,
            1.0/(molr1/ct1.rPr_ - molr2/ct2.rPr_)
        );
    }
}


template<class Thermo>
inline Foam::CO2Transport<Thermo> Foam::operator*
(
    const scalar s,
    const CO2Transport<Thermo>& ct
)
{
    return CO2Transport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        1.0/ct.rPr_
    );
}


template<class Thermo>
inline Foam::CO2Transport<Thermo> Foam::operator==
(
    const CO2Transport<Thermo>& ct1,
    const CO2Transport<Thermo>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
