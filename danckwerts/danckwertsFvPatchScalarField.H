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

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides static pressure condition for p_rgh,
    calculated as:

        \f[
            p_rgh = p - \rho g h
        \f]

    where
    \vartable
        p_rgh   | Pseudo hydrostatic pressure [Pa]
        p       | Static pressure [Pa]
        h       | Height in the opposite direction to gravity
        \rho    | density
        g       | acceleration due to gravity [m/s2]
    \endtable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        rhoName      | rho field name          | no          | rho
        p            | static pressure         | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            danckwerts;
        p               uniform 0;
        value           uniform 0; // optional initial value
    }
    \endverbatim

SeeAlso
    Foam::mixedFvPatchScalarField

SourceFiles
    danckwertsFvPatchScalarField.C

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

    // Protected data
  
  
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
        // Evaluation functions
        virtual void evaluate( const Pstream::commsTypes commsType=Pstream::blocking );

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
