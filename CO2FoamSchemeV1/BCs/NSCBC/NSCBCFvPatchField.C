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

#include "NSCBCFvPatchField.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::NSCBCFvPatchField<Type>::NSCBCFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::NSCBCFvPatchField<Type>::NSCBCFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict)
{
    fvPatchField<Type>::operator=(this->patchInternalField());
}


template<class Type>
Foam::NSCBCFvPatchField<Type>::NSCBCFvPatchField
(
    const NSCBCFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::NSCBCFvPatchField<Type>::NSCBCFvPatchField
(
    const NSCBCFvPatchField& zgpf
)
:
    fvPatchField<Type>(zgpf)
{}


template<class Type>
Foam::NSCBCFvPatchField<Type>::NSCBCFvPatchField
(
    const NSCBCFvPatchField& zgpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(zgpf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::NSCBCFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{

    Info<< "Inside evaluate before updateCoeffs" << endl;

    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<Type>::operator==(this->patchInternalField());
//    fvPatchField<Type>::operator==patch().boundaryFieldRef();
    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::NSCBCFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), pTraits<Type>::one)
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::NSCBCFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::NSCBCFvPatchField<Type>::gradientInternalCoeffs() const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::NSCBCFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}


// ************************************************************************* //
