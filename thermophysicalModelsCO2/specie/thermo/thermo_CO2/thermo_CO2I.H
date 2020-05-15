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

#include "thermo_CO2.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::thermo_CO2<Thermo, Type>::thermo_CO2
(
    const Thermo& sp
)
:
    Thermo(sp)
{}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo_CO2<Thermo, Type>::T
(
    scalar f,
    scalar p,
    scalar T0,
    scalar (thermo_CO2<Thermo, Type>::*F)(const scalar, const scalar) const,
    scalar (thermo_CO2<Thermo, Type>::*dFdT)(const scalar, const scalar)
        const,
    scalar (thermo_CO2<Thermo, Type>::*limit)(const scalar) const
) const
{
    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    int    iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            (this->*limit)
            (Test - ((this->*F)(p, Test) - f)/(this->*dFdT)(p, Test));

        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded"
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::thermo_CO2<Thermo, Type>::thermo_CO2
(
    const word& name,
    const thermo_CO2& st
)
:
    Thermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::word
Foam::species::thermo_CO2<Thermo, Type>::heName()
{
    return Type<thermo_CO2<Thermo, Type>>::name();
}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar
//Foam::species::thermo_CO2<Thermo, Type>::he(const scalar p, const scalar T) const
//{
//    return Type<thermo_CO2<Thermo, Type>>::he(*this, p, T);
//}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar
//Foam::species::thermo_CO2<Thermo, Type>::cv(const scalar p, const scalar T) const
//{
//    return this->cp(p, T) - this->cpMcv(p, T);
//}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar
//Foam::species::thermo_CO2<Thermo, Type>::cpv(const scalar p, const scalar T) const
//{
//    return Type<thermo_CO2<Thermo, Type>>::cpv(*this, p, T);
//}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar
//Foam::species::thermo_CO2<Thermo, Type>::gamma(const scalar p, const scalar T) const
//{
//    scalar cp = this->cp(p, T);
//    return cp/(cp - this->cpMcv(p, T));
//}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar
//Foam::species::thermo_CO2<Thermo, Type>::cpBycpv
//(
//    const scalar p,
//    const scalar T
//) const
//{
//    return Type<thermo_CO2<Thermo, Type>>::cpBycpv(*this, p, T);
//}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::es(const scalar p, const scalar T) const
{
    return this->hs(p, T) - p*this->W()/this->rho(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::ea(const scalar p, const scalar T) const
{
    return this->ha(p, T) - p*this->W()/this->rho(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::g(const scalar p, const scalar T) const
{
    return this->ha(p, T) - T*this->s(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::a(const scalar p, const scalar T) const
{
    return this->ea(p, T) - T*this->s(p, T);
}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar
//Foam::species::thermo_CO2<Thermo, Type>::Cpv(const scalar p, const scalar T) const
//{
//    return this->cpv(p, T)/this->W();
//}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar
//Foam::species::thermo_CO2<Thermo, Type>::Cp(const scalar e, const scalar v) const
//{
//    return 1.0;
//}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar
//Foam::species::thermo_CO2<Thermo, Type>::Cv(const scalar e, const scalar v) const
//{
//    return 1.0;
//}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar
//Foam::species::thermo_CO2<Thermo, Type>::HE(const scalar e, const scalar v) const
//{
//    return 1.0;
//}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::Hs(const scalar p, const scalar T) const
{
    return this->hs(p, T)/this->W();
}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar
//Foam::species::thermo_CO2<Thermo, Type>::Hc() const
//{
//    return 1.0;
//}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::Ha(const scalar p, const scalar T) const
{
    return this->ha(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::S(const scalar p, const scalar T) const
{
    return this->s(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::E(const scalar p, const scalar T) const
{
    return this->e(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::Es(const scalar p, const scalar T) const
{
    return this->es(p, T)/this->W();
}

template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::Ea(const scalar p, const scalar T) const
{
    return this->ea(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::G(const scalar p, const scalar T) const
{
    return this->g(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::A(const scalar p, const scalar T) const
{
    return this->a(p, T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::K(const scalar p, const scalar T) const
{
    scalar arg = -this->nMoles()*this->g(Pstd, T)/(RR*T);

    if (arg < 600.0)
    {
        return exp(arg);
    }
    else
    {
        return VGREAT;
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::Kp(const scalar p, const scalar T) const
{
    return K(p, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::species::thermo_CO2<Thermo, Type>::Kc(const scalar p, const scalar T) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(Pstd/(RR*T), this->nMoles());
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo_CO2<Thermo, Type>::Kx
(
    const scalar p,
    const scalar T
) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(Pstd/p, this->nMoles());
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::species::thermo_CO2<Thermo, Type>::Kn
(
    const scalar p,
    const scalar T,
    const scalar n
) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(n*Pstd/p, this->nMoles());
    }
}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar Foam::species::thermo_CO2<Thermo, Type>::THE
//(
//    const scalar he,
//    const scalar p,
//    const scalar T0
//) const
//{
//    return 1.0;
//}


//template<class Thermo, template<class> class Type>
//inline Foam::scalar Foam::species::thermo_CO2<Thermo, Type>::THs
//(
//    const scalar hs,
//    const scalar p,
//    const scalar T0
//) const
//{
//    return T
//    (
//        hs,
//        p,
//        T0,
//        &thermo_CO2<Thermo, Type>::Hs,
//        &thermo_CO2<Thermo, Type>::Cp,
//        &thermo_CO2<Thermo, Type>::limit
//    );
//}

template<class Thermo, template<class> class Type>
inline void Foam::species::thermo_CO2<Thermo, Type>::operator+=
(
    const thermo_CO2<Thermo, Type>& st
)
{
    Thermo::operator+=(st);
}


template<class Thermo, template<class> class Type>
inline void Foam::species::thermo_CO2<Thermo, Type>::operator-=
(
    const thermo_CO2<Thermo, Type>& st
)
{
    Thermo::operator-=(st);
}


template<class Thermo, template<class> class Type>
inline void Foam::species::thermo_CO2<Thermo, Type>::operator*=(const scalar s)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::species::thermo_CO2<Thermo, Type> Foam::species::operator+
(
    const thermo_CO2<Thermo, Type>& st1,
    const thermo_CO2<Thermo, Type>& st2
)
{
    return thermo_CO2<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::species::thermo_CO2<Thermo, Type> Foam::species::operator-
(
    const thermo_CO2<Thermo, Type>& st1,
    const thermo_CO2<Thermo, Type>& st2
)
{
    return thermo_CO2<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::species::thermo_CO2<Thermo, Type> Foam::species::operator*
(
    const scalar s,
    const thermo_CO2<Thermo, Type>& st
)
{
    return thermo_CO2<Thermo, Type>
    (
        s*static_cast<const Thermo&>(st)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::species::thermo_CO2<Thermo, Type> Foam::species::operator==
(
    const thermo_CO2<Thermo, Type>& st1,
    const thermo_CO2<Thermo, Type>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
