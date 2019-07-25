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

\*---------------------------------------------------------------------------*/

#include "normalMotionSlipBasePointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointPatchFields.H"
#include "IFstream.H"
#include "symmTransformField.H"

#include "coupledPatchInterpolation.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::normalMotionSlipBasePointPatchVectorField::normalMotionSlipBasePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    valuePointPatchField<vector>(p, iF),
    rlxON(true)
{}


Foam::normalMotionSlipBasePointPatchVectorField::normalMotionSlipBasePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    valuePointPatchField<vector>(p, iF, dict),
    rlxON(true)
{
}


Foam::normalMotionSlipBasePointPatchVectorField::normalMotionSlipBasePointPatchVectorField
(
    const normalMotionSlipBasePointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    valuePointPatchField<vector>(ptf, p, iF, mapper),
    rlxON(ptf.rlxON),
    faceDispl( ptf.faceDispl )
{}

Foam::normalMotionSlipBasePointPatchVectorField::normalMotionSlipBasePointPatchVectorField
(
    const normalMotionSlipBasePointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    valuePointPatchField<vector>(ptf, iF),
    rlxON(ptf.rlxON),
    faceDispl( ptf.faceDispl )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::normalMotionSlipBasePointPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if(debug) 
    {
        Info<<"   normalMotionSlipBasePointPatchVectorField::evaluate()"<<nl;
    }
  
    const polyMesh& mesh = this->internalField().mesh()();
    label patchID = this->patch().index();
    const fvMesh& fvmesh_ = refCast<const fvMesh>(mesh);
    coupledPatchInterpolation patchInterpolator
    ( 
        mesh.boundaryMesh()[patchID], fvmesh_
    );
  
    const volVectorField& cmu =
            this->db().objectRegistry::lookupObject<volVectorField>
            (
                "cellMotionU"
            );
    /*
    //const fvMesh& mesh = this->internalField().mesh();
    const polyPatch& curPatch = mesh.boundaryMesh()[patchID];
    const List<face>& llf = curPatch.localFaces();
    const pointField& curPP  = curPatch.localPoints();

    vectorField pMotion = patchInterpolator.faceToPointInterpolate(this->getDisp());
    
    pointField movedPoints(curPP + pMotion);

    vectorField fn( llf.size() );
    forAll(fn, facei)
    {
      fn[facei]  = llf[facei].normal(movedPoints);
      fn[facei] /= mag(fn[facei]) + VSMALL;
    }

    vectorField proj
    (
    //cmu.boundaryField()[patchID].patchInternalField()
    transform(I - sqr(fn), cmu.boundaryField()[patchID].patchInternalField())
    );
    */
    
    
    this->operator==
    ( 
        patchInterpolator.faceToPointInterpolate
        (
            cmu.boundaryField()[patchID]
            //+
            //proj
        )
    );
    
    valuePointPatchField<vector>::evaluate();
}

void Foam::normalMotionSlipBasePointPatchVectorField::updateCoeffs()
{
    if(debug) 
    {
        Info<<"   normalMotionSlipBasePointPatchVectorField::updateCoeffs()"<<nl;
    }

    valuePointPatchField<vector>::updateCoeffs();
}


void Foam::normalMotionSlipBasePointPatchVectorField::write
(
  Ostream& os
) const
{
    valuePointPatchField<vector>::write(os);
}


namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        normalMotionSlipBasePointPatchVectorField
    );
}


// ************************************************************************* //
