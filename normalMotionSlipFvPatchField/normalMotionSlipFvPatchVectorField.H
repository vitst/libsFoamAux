/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-201X OpenFOAM Foundation
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
    Foam::normalMotionSlipFvPatchVectorField

Description

    This boundary condition provides a velocity condition for cases with
    moving walls.  In addition, it should also be applied to 'moving' walls
    for moving reference frame (MRF) calculations.


    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    patchName
    {
        type            normalMotionSlip;
        value           uniform (0 0 0);   // initial value
    }
    \endverbatim

SourceFiles
    normalMotionSlipFvPatchVectorField.C
    

\*---------------------------------------------------------------------------*/

#ifndef normalMotionSlipFvPatchVectorField_H
#define normalMotionSlipFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class normalMotionSlipFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class normalMotionSlipFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
  
  
public:

    //- Runtime type information
    TypeName("normalMotionSlip");

    
    // Constructors

        //- Construct from patch and internal field
        normalMotionSlipFvPatchVectorField
        (
          const fvPatch&,
          const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        normalMotionSlipFvPatchVectorField
        (
          const fvPatch&,
          const DimensionedField<vector, volMesh>&,
          const dictionary&
        );

        //- Construct by mapping given normalMotionSlipFvPatchVectorField
        //  onto a new patch
        normalMotionSlipFvPatchVectorField
        (
          const normalMotionSlipFvPatchVectorField&,
          const fvPatch&,
          const DimensionedField<vector, volMesh>&,
          const fvPatchFieldMapper&
        );

        //- Construct as copy
        normalMotionSlipFvPatchVectorField
        (
            const normalMotionSlipFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new normalMotionSlipFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        normalMotionSlipFvPatchVectorField
        (
            const normalMotionSlipFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new normalMotionSlipFvPatchVectorField(*this, iF)
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

#endif

// ************************************************************************* //
