/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "codedNormalMotionSlipPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "pointPatchFieldMapper.H"
#include "pointFields.H"

#include "volFields.H"
#include "pointPatchFields.H"

#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//template<Type>
const Foam::word Foam::codedNormalMotionSlipPointPatchVectorField::codeTemplateC
    = "normalMotionSlipPointPatchFieldTemplate.C";

//template<Type>
const Foam::word Foam::codedNormalMotionSlipPointPatchVectorField::codeTemplateH
    = "normalMotionSlipPointPatchFieldTemplate.H";

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::codedNormalMotionSlipPointPatchVectorField::setFieldTemplates
(
    dynamicCode& dynCode
)
{
    word fieldType("vector");

    // Template type for pointPatchField
    dynCode.setFilterVariable("TemplateType", fieldType);

    // Name for pointPatchField - eg, ScalarField, VectorField, ...
    fieldType[0] = toupper(fieldType[0]);
    dynCode.setFilterVariable("FieldType", fieldType + "Field");
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::IOdictionary& Foam::codedNormalMotionSlipPointPatchVectorField::dict()
const
{
    const objectRegistry& obr = this->db();

    if (obr.foundObject<IOdictionary>("codeDict"))
    {
        return obr.lookupObject<IOdictionary>("codeDict");
    }
    else
    {
        return obr.store
        (
            new IOdictionary
            (
                IOobject
                (
                    "codeDict",
                    this->db().time().system(),
                    this->db(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


Foam::dlLibraryTable& Foam::codedNormalMotionSlipPointPatchVectorField::libs() const
{
    return const_cast<dlLibraryTable&>(this->db().time().libs());
}


void Foam::codedNormalMotionSlipPointPatchVectorField::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Take no chances - typeName must be identical to name_
    dynCode.setFilterVariable("typeName", name_);

    // Set TemplateType and FieldType filter variables
    // (for pointPatchField)
    setFieldTemplates(dynCode);

    // Compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH);


    // Debugging: make BC verbose
    //   dynCode.setFilterVariable("verbose", "true");
    //   Info<<"compile " << name_ << " sha1: "
    //       << context.sha1() << endl;

    // Define Make/options
    dynCode.setMakeOptions
        (
            "EXE_INC = -g \\\n"
            "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
            + context.options()
            + "\n\nLIB_LIBS = \\\n"
            + "    -lOpenFOAM \\\n"
            + "    -lfiniteVolume \\\n"
            + context.libs()
        );
}


const Foam::dictionary& Foam::codedNormalMotionSlipPointPatchVectorField::codeDict()
const
{
    // Use system/codeDict or in-line
    return
    (
        dict_.found("code")
      ? dict_
      : this->dict().subDict(name_)
    );
}


Foam::string Foam::codedNormalMotionSlipPointPatchVectorField::description() const
{
    return
        "patch "
      + this->patch().name()
      + " on field "
      + this->internalField().name();
}


void Foam::codedNormalMotionSlipPointPatchVectorField::clearRedirect() const
{
    // Remove instantiation of pointPatchField provided by library
    redirectPatchFieldPtr_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedNormalMotionSlipPointPatchVectorField::codedNormalMotionSlipPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    normalMotionSlipBasePointPatchVectorField(p, iF),
    codedBase(),
    redirectPatchFieldPtr_()
{}


Foam::codedNormalMotionSlipPointPatchVectorField::codedNormalMotionSlipPointPatchVectorField
(
    const codedNormalMotionSlipPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    normalMotionSlipBasePointPatchVectorField(ptf, p, iF, mapper),
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
    redirectPatchFieldPtr_()
{}


Foam::codedNormalMotionSlipPointPatchVectorField::codedNormalMotionSlipPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    normalMotionSlipBasePointPatchVectorField(p, iF, dict),
    codedBase(),
    dict_(dict),
    name_
    (
        dict.found("redirectType")
      ? dict.lookup("redirectType")
      : dict.lookup("name")
    ),
    redirectPatchFieldPtr_()
{
    updateLibrary(name_);
}


Foam::codedNormalMotionSlipPointPatchVectorField::codedNormalMotionSlipPointPatchVectorField
(
    const codedNormalMotionSlipPointPatchVectorField& ptf
)
:
    normalMotionSlipBasePointPatchVectorField(ptf),
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
    redirectPatchFieldPtr_()
{}


Foam::codedNormalMotionSlipPointPatchVectorField::codedNormalMotionSlipPointPatchVectorField
(
    const codedNormalMotionSlipPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    normalMotionSlipBasePointPatchVectorField(ptf, iF),
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
    redirectPatchFieldPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//const Foam::pointPatchField<Foam::vector>&
const Foam::normalMotionSlipBasePointPatchVectorField&
Foam::codedNormalMotionSlipPointPatchVectorField::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        // Construct a patch
        // Make sure to construct the patchfield with up-to-date value

        OStringStream os;
        //os.writeEntry("type", name_);
        os.writeKeyword("type") << name_ << token::END_STATEMENT << nl;
        static_cast<const Field<vector>&>(*this).writeEntry("value", os);
        IStringStream is(os.str());
        dictionary dict(is);

        redirectPatchFieldPtr_.reset
        (
            dynamic_cast<normalMotionSlipBasePointPatchVectorField*>
            (
                pointPatchField<vector>::New
                (
                    this->patch(),
                    this->internalField(),
                    dict
                ).ptr()
            )
        );
    }
    return redirectPatchFieldPtr_();
}


void Foam::codedNormalMotionSlipPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Make sure library containing user-defined pointPatchField is up-to-date
    updateLibrary(name_);

    //const pointPatchField<vector>& fvp = redirectPatchField();
    const normalMotionSlipBasePointPatchVectorField& fvp = redirectPatchField();

    //const_cast<pointPatchField<vector>&>(fvp).updateCoeffs();
    const_cast<normalMotionSlipBasePointPatchVectorField&>(fvp).updateCoeffs();

    // Copy through value
    //this->operator==(fvp);
    this->setDisp(fvp.getDisp());

    normalMotionSlipBasePointPatchVectorField::updateCoeffs();
}


void Foam::codedNormalMotionSlipPointPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Make sure library containing user-defined pointPatchField is up-to-date
    updateLibrary(name_);

    //const pointPatchField<vector>& fvp = redirectPatchField();
    const normalMotionSlipBasePointPatchVectorField& fvp = redirectPatchField();

    //const_cast<pointPatchField<vector>&>(fvp).evaluate(commsType);
    const_cast<normalMotionSlipBasePointPatchVectorField&>(fvp).evaluate(commsType);

    normalMotionSlipBasePointPatchVectorField::evaluate(commsType);
}


void Foam::codedNormalMotionSlipPointPatchVectorField::write(Ostream& os) const
{
    normalMotionSlipBasePointPatchVectorField::write(os);
    //os.writeEntry("name", name_);
    os.writeKeyword("type") << name_ << token::END_STATEMENT << nl;

    //codedBase::writeCodeDict(os, dict_);
    if (dict_.found("codeInclude"))
    {
        os.writeKeyword("codeInclude")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["codeInclude"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("localCode"))
    {
        os.writeKeyword("localCode")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["localCode"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("code"))
    {
        os.writeKeyword("code")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["code"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("codeOptions"))
    {
        os.writeKeyword("codeOptions")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["codeOptions"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("codeLibs"))
    {
        os.writeKeyword("codeLibs")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["codeLibs"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

}


namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        codedNormalMotionSlipPointPatchVectorField
    );
}

// ************************************************************************* //
