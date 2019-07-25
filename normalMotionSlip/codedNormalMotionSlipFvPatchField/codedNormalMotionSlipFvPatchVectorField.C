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

#include "codedNormalMotionSlipFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "symmTransformField.H"

#include "codedNormalMotionSlipPointPatchVectorField.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedNormalMotionSlipFvPatchVectorField::
codedNormalMotionSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    normalMotionSlipFvPatchVectorField(p, iF)
{}


Foam::codedNormalMotionSlipFvPatchVectorField::
codedNormalMotionSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    normalMotionSlipFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::codedNormalMotionSlipFvPatchVectorField::
codedNormalMotionSlipFvPatchVectorField
(
    const codedNormalMotionSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    normalMotionSlipFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::codedNormalMotionSlipFvPatchVectorField::
codedNormalMotionSlipFvPatchVectorField
(
    const codedNormalMotionSlipFvPatchVectorField& mwvpvf
)
:
    normalMotionSlipFvPatchVectorField(mwvpvf)
{}


Foam::codedNormalMotionSlipFvPatchVectorField::
codedNormalMotionSlipFvPatchVectorField
(
    const codedNormalMotionSlipFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    normalMotionSlipFvPatchVectorField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::codedNormalMotionSlipFvPatchVectorField::updateCoeffs()
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
    
    const Foam::codedNormalMotionSlipPointPatchVectorField& pmuBC =
            refCast<const Foam::codedNormalMotionSlipPointPatchVectorField>
            (
                pmU.boundaryField()[patchID]
            );
    
    //vectorField disp( pmuBC.getDisp() );
    bool rlxON = pmuBC.getRlxON();
    if(rlxON)
    {
        vectorField::operator=
        ( 
            pmuBC.getDisp() + transform(I - sqr(n), this->patchInternalField())
        );
    }
    else
    {
        vectorField::operator=
        ( 
            pmuBC.getDisp()
        );
    }

    fixedValueFvPatchVectorField::updateCoeffs();
  
}
 
 
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        codedNormalMotionSlipFvPatchVectorField
    );
}

// ************************************************************************* //
