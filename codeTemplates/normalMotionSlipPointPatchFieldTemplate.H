/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR AUTHOR,AFFILIATION
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

Description
    Template for use with dynamic code generation of a
    normalMotionSlip pointPatchField.

    - without state

SourceFiles
    normalMotionSlipPointPatchFieldTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef fixedValuePointPatchTemplate${FieldType}_H
#define normalMotionSlipPointPatchTemplate${FieldType}_H

#include "normalMotionSlipBasePointPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                A templated NormalMotionSlipPointPatchField
\*---------------------------------------------------------------------------*/

class ${typeName}NormalMotionSlipPointPatch${FieldType}
:
    public normalMotionSlipBasePointPatchVectorField
{
public:

    //- Information about the SHA1 of the code itself
    static const char* const SHA1sum;

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from patch and internal field
        ${typeName}NormalMotionSlipPointPatch${FieldType}
        (
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ${typeName}NormalMotionSlipPointPatch${FieldType}
        (
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping a copy onto a new patch
        ${typeName}NormalMotionSlipPointPatch${FieldType}
        (
            const ${typeName}NormalMotionSlipPointPatch${FieldType}&,
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        ${typeName}NormalMotionSlipPointPatch${FieldType}
        (
            const ${typeName}NormalMotionSlipPointPatch${FieldType}&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<${TemplateType}>> clone() const
        {
            return autoPtr<pointPatchField<${TemplateType}>>
            (
                new ${typeName}NormalMotionSlipPointPatch${FieldType}(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ${typeName}NormalMotionSlipPointPatch${FieldType}
        (
            const ${typeName}NormalMotionSlipPointPatch${FieldType}&,
            const DimensionedField<${TemplateType}, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<${TemplateType}>> clone
        (
            const DimensionedField<${TemplateType}, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<${TemplateType}>>
            (
                new ${typeName}NormalMotionSlipPointPatch${FieldType}(*this, iF)
            );
        }


    //- Destructor
    virtual ~${typeName}NormalMotionSlipPointPatch${FieldType}();


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
