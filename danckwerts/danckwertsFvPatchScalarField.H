/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::danckwertsFvPatchScalarField

Description
    This boundary condition provides danckwerts condition for scalar field
    calculated as:

        (u*n)(1-F) + D*n*grad(F) = 0

    where
    \vartable
        F       | field
        u       | velocity field
        D       | diffusion coefficient
        n       | normals
    \endtable

    \heading Patch usage

    \table
        Property     | Description             | Required
        u            | field U                 | yes
        D            | diffusion coef D        | yes
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            danckwerts;
        value           uniform 0; // optional initial value
    }
    \endverbatim

SeeAlso
    Foam::mixedFvPatchScalarField

SourceFiles
    danckwertsFvPatchScalarField.C

Contributors/Copyright:
    Vitaliy Starchenko (2015) <vitaliy.starchenko@gmail.com>

\*---------------------------------------------------------------------------*/

#ifndef danckwertsFvPatchScalarField_H
#define danckwertsFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class danckwertsFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class danckwertsFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

public:

    //- Runtime type information
    TypeName("danckwerts");


    // Constructors

        //- Construct from patch and internal field
        danckwertsFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        danckwertsFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  danckwertsFvPatchScalarField onto a new patch
        
        danckwertsFvPatchScalarField
        (
            const danckwertsFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        danckwertsFvPatchScalarField
        (
            const danckwertsFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new danckwertsFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        danckwertsFvPatchScalarField
        (
            const danckwertsFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new danckwertsFvPatchScalarField(*this, iF)
            );
        }

    // Member functions
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
