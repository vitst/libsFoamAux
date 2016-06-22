/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
  
Class
    Foam::danckwertsFvPatchScalarField

Description
    This is the patch for nonlinear boundary concentration field.
    Where reaction kinetics is proportional to:

        ~ c^n_1  if c<c_th
        ~ c^n_2  if c>c_th

    where
    \vartable
        c       | concentration field
        c_th    | threshold
        n_1     | exponent 1 for nonlinear kinetics
        n_2     | exponent 2 for nonlinear kinetics
    \endtable

    \heading Patch usage

    \table
        Property     | Description                                 | Required
        l_T          | thickness of concentration boundary layer   | yes
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
      type            nonLinear;
      value           uniform 0.0;
      Cth   0.2;
      n1    1.0;
      n2    4.0;
    }
    \endverbatim

SeeAlso
    Foam::mixedFvPatchScalarField

SourceFiles
    danckwertsFvPatchScalarField.C

  
Contributors/Copyright:
    Vitaliy Starchenko (2015) <vitaliy.starchenko@gmail.com>

\*---------------------------------------------------------------------------*/

#ifndef nonLinearFvPatchField_H
#define nonLinearFvPatchField_H

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class nonLinearFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class nonLinearFvPatchField
:
    public mixedFvPatchField<Type>
{

protected:
    // Protected data
    scalar l_T;
    
    scalar Cth; // threshold concentration
    
    scalar n1, n2; // power low for the dissolution kinetics
public:

    //- Runtime type information
    TypeName("nonLinear");

    // Constructors

    //- Construct from patch and internal field
    nonLinearFvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    nonLinearFvPatchField
    (
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given nonLinearFvPatchField onto a new patch
    nonLinearFvPatchField
    (
        const nonLinearFvPatchField<Type>&,
        const fvPatch&,
        const DimensionedField<Type, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    nonLinearFvPatchField
    (
        const nonLinearFvPatchField<Type>&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
                (
                    new nonLinearFvPatchField<Type>(*this)
                );
        }

    //- Construct as copy setting internal field reference
    nonLinearFvPatchField
    (
        const nonLinearFvPatchField<Type>&,
        const DimensionedField<Type, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchField<Type> > clone
    (
        const DimensionedField<Type, volMesh>& iF
    ) const
        {
            return tmp<fvPatchField<Type> >
                (
                    new nonLinearFvPatchField<Type>(*this, iF)
                );
        }

    // Member functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #ifdef NoRepository
    #   include "nonLinearFvPatchField.C"
    #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #endif

// ************************************************************************* //
