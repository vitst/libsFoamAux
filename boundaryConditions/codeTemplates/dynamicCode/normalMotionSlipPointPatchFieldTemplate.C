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

\*---------------------------------------------------------------------------*/

#include "normalMotionSlipPointPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "pointPatchFieldMapper.H"
#include "pointFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void ${typeName}_${SHA1sum}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatch${FieldType},
    ${typeName}NormalMotionSlipPointPatch${FieldType}
);


const char* const ${typeName}NormalMotionSlipPointPatch${FieldType}::SHA1sum =
    "${SHA1sum}";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}NormalMotionSlipPointPatch${FieldType}::
${typeName}NormalMotionSlipPointPatch${FieldType}
(
    const pointPatch& p,
    const DimensionedField<${TemplateType}, pointMesh>& iF
)
:
//    normalMotionSlipPointPatchField<${TemplateType}>(p, iF)
    normalMotionSlipBasePointPatchVectorField(p, iF)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/DimensionedField\n";
    }
}


${typeName}NormalMotionSlipPointPatch${FieldType}::
${typeName}NormalMotionSlipPointPatch${FieldType}
(
    const ${typeName}NormalMotionSlipPointPatch${FieldType}& ptf,
    const pointPatch& p,
    const DimensionedField<${TemplateType}, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    normalMotionSlipBasePointPatchVectorField(ptf, p, iF, mapper)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/DimensionedField/mapper\n";
    }
}


${typeName}NormalMotionSlipPointPatch${FieldType}::
${typeName}NormalMotionSlipPointPatch${FieldType}
(
    const pointPatch& p,
    const DimensionedField<${TemplateType}, pointMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    //normalMotionSlipBasePointPatchVectorField(p, iF, dict, valueRequired)
    normalMotionSlipBasePointPatchVectorField(p, iF, dict)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " from patch/dictionary\n";
    }
}


${typeName}NormalMotionSlipPointPatch${FieldType}::
${typeName}NormalMotionSlipPointPatch${FieldType}
(
    const ${typeName}NormalMotionSlipPointPatch${FieldType}& ptf
)
:
    normalMotionSlipBasePointPatchVectorField(ptf)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum}"
            " as copy\n";
    }
}


${typeName}NormalMotionSlipPointPatch${FieldType}::
${typeName}NormalMotionSlipPointPatch${FieldType}
(
    const ${typeName}NormalMotionSlipPointPatch${FieldType}& ptf,
    const DimensionedField<${TemplateType}, pointMesh>& iF
)
:
    normalMotionSlipBasePointPatchVectorField(ptf, iF)
{
    if (${verbose:-false})
    {
        Info<<"construct ${typeName} sha1: ${SHA1sum} "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}NormalMotionSlipPointPatch${FieldType}::
~${typeName}NormalMotionSlipPointPatch${FieldType}()
{
    if (${verbose:-false})
    {
        Info<<"destroy ${typeName} sha1: ${SHA1sum}\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ${typeName}NormalMotionSlipPointPatch${FieldType}::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (${verbose:-false})
    {
        Info<<"updateCoeffs ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${code}
//}}} end code

    this->normalMotionSlipBasePointPatchVectorField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
