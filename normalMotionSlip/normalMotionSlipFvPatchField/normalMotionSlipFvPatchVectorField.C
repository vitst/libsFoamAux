/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "normalMotionSlipFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "symmTransformField.H"

#include "normalMotionSlipPointPatchVectorField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::normalMotionSlipFvPatchVectorField::
normalMotionSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::normalMotionSlipFvPatchVectorField::
normalMotionSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::normalMotionSlipFvPatchVectorField::
normalMotionSlipFvPatchVectorField
(
    const normalMotionSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::normalMotionSlipFvPatchVectorField::
normalMotionSlipFvPatchVectorField
(
    const normalMotionSlipFvPatchVectorField& mwvpvf
)
:
    fixedValueFvPatchVectorField(mwvpvf)
{}


Foam::normalMotionSlipFvPatchVectorField::
normalMotionSlipFvPatchVectorField
(
    const normalMotionSlipFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  
void Foam::normalMotionSlipFvPatchVectorField::updateCoeffs()
{
    if(debug)
    {
        Info<<"   normalMotionSlipFvPatchVectorField::updateCoeffs()"<<nl;
    }

    if (updated())
    {
        return;
    }
    
    const vectorField n(this->patch().nf());
    label patchID = this->patch().index();

    const pointVectorField& pmU = 
            this->db().objectRegistry::lookupObject<pointVectorField>
            (
                "pointMotionU"
            );

    const Foam::normalMotionSlipPointPatchVectorField& pmuBC =
            refCast<const Foam::normalMotionSlipPointPatchVectorField>
            (
                pmU.boundaryField()[patchID]
            );
    vectorField disp( pmuBC.getDisp() );

    vectorField::operator=
    ( 
        pmuBC.getDisp() + transform(I - sqr(n), this->patchInternalField())
    );

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::normalMotionSlipFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        normalMotionSlipFvPatchVectorField
    );
}

// ************************************************************************* //
