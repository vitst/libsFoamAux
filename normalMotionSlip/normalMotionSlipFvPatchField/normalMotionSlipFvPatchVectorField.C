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
#include "coupledPatchInterpolation.H"

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
    // calculate updated face normals and face centers from new point positions
    // without the relaxation
    //const polyMesh& mesh = this->internalField().mesh();
    //const fvMesh& mesh = this->internalField().mesh();
    //const polyPatch& curPatch = mesh.boundaryMesh()[patchID];
    //const List<face>& llf = curPatch.localFaces();
    //const pointField& curPP  = curPatch.localPoints();
    
    //const fvMesh& fvmesh_ = refCast<const fvMesh>(mesh);
    //const fvMesh& fvmesh_ = this->internalField().mesh();// patchField_.patch().boundaryMesh().mesh();
    
    //coupledPatchInterpolation patchInterpolator
    //(
    //  curPatch, fvmesh_
    //);
    //vectorField pMotion = patchInterpolator.faceToPointInterpolate(pmuBC.getDisp());
    
    //pointField movedPoints(curPP + pMotion);

    // normals for the displaced mesh without relaxation
    //const vectorField n(nextFaceNormals(movedPoints, llf));
    const vectorField n(this->patch().nf());
    // cell centers for the displaced mesh without relaxation
    //const pointField fc(faceCentres(movedPoints, llf));
    
    //vectorField aux(pmuBC.getDisp() + transform(I - sqr(n), this->patchInternalField()));

    vectorField::operator=
    ( 
        pmuBC.getDisp() + transform(I - sqr(n), this->patchInternalField())
        //aux + transform(I - sqr(n), this->patchInternalField())
    );
    //pmuBC.setDisp(aux + transform(I - sqr(n), this->patchInternalField()));

    fixedValueFvPatchVectorField::updateCoeffs();
}

//Foam::tmp<Foam::vectorField>
Foam::vectorField
Foam::normalMotionSlipFvPatchVectorField::
nextFaceNormals(const pointField& points, const List<face>& flist) const
{
  vectorField fn( flist.size() );
  forAll(fn, facei)
  {
    fn[facei]  = flist[facei].normal(points);
    fn[facei] /= mag(fn[facei]) + VSMALL;
  }
  return fn;
}

Foam::pointField
Foam::normalMotionSlipFvPatchVectorField::
nextFaceCentres(const pointField& points, const List<face>& flist) const
{
  pointField fc( flist.size() );
  forAll(fc, facei)
  {
    fc[facei] = flist[facei].centre(points);
  }
  return fc;
}


void Foam::normalMotionSlipFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if(debug)
    {
        Info<<"   normalMotionSlipFvPatchVectorField::evaluate()"<<nl;
    }
    fixedValueFvPatchVectorField::evaluate();
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
