/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::normalMotionSlipPointPatchVectorField

Description
    This patch can be used for the mesh motion using velocityLaplacian 
    dynamic mesh solver. It should be used along with 
    normalMotionSlipPointPatchVectorField for cellMotionU volume filed. 
    It saves face displacement at the begin of Laplace solver. At the end it
    interpolates displacement to the point field.

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    patchName
    {
        type            normalMotionSlip;
        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    normalMotionSlipPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef normalMotionSlipPointPatchVectorField_H
#define normalMotionSlipPointPatchVectorField_H

#include "valuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class normalMotionSlipPointPatchVectorField;


Ostream& operator<<
(
    Ostream&,
    const normalMotionSlipPointPatchVectorField&
);


/*---------------------------------------------------------------------------*\
                        Class normalMotionSlipPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class normalMotionSlipPointPatchVectorField
:
    public valuePointPatchField<vector>
{
    // Private data
    
    //- saves displacement field on the boundary for cellMotionU
    vectorField faceDispl;
    
    word fieldName;
    
    word scalarName;

public:

    //- Runtime type information
    TypeName("normalMotionSlip");


    // Constructors

        //- Construct from patch and internal field
        normalMotionSlipPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        normalMotionSlipPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        normalMotionSlipPointPatchVectorField
        (
            const normalMotionSlipPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new normalMotionSlipPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        normalMotionSlipPointPatchVectorField
        (
            const normalMotionSlipPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new normalMotionSlipPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions
        const Field<vector>& getDisp() const 
        { 
          return faceDispl; 
        }
        
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
        
        //- Update the patch field
        virtual void evaluate
        (
          const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );
        
        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "normalMotionSlipPointPatchVectorField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
