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

#include "normalMotionSlipPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointPatchFields.H"

#include "coupledPatchInterpolation.H"
#include "fixedNormalSlipPointPatchField.H"
#include "cyclicSlipPointPatchField.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::normalMotionSlipPointPatchVectorField::normalMotionSlipPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    valuePointPatchField<vector>(p, iF)
{}


Foam::normalMotionSlipPointPatchVectorField::normalMotionSlipPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    valuePointPatchField<vector>(p, iF, dict)
{}


Foam::normalMotionSlipPointPatchVectorField::normalMotionSlipPointPatchVectorField
(
    const normalMotionSlipPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    valuePointPatchField<vector>(ptf, p, iF, mapper)
{}

Foam::normalMotionSlipPointPatchVectorField::normalMotionSlipPointPatchVectorField
(
    const normalMotionSlipPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    valuePointPatchField<vector>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::normalMotionSlipPointPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    if(debug) 
    {
        Info<<"   normalMotionSlipPointPatchVectorField::evaluate()"<<nl;
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
    
    this->operator==
    ( 
        patchInterpolator.faceToPointInterpolate(cmu.boundaryField()[patchID]) 
    );
    
    valuePointPatchField<vector>::evaluate();
}

void Foam::normalMotionSlipPointPatchVectorField::updateCoeffs()
{
    if(debug) 
    {
        Info<<"   normalMotionSlipPointPatchVectorField::updateCoeffs()"<<nl;
    }

    const polyMesh& mesh = this->internalField().mesh()();

    const volScalarField& C =
            this->db().objectRegistry::lookupObject<volScalarField>("C");
    const label& patchID = this->patch().index();
    const IOdictionary& IOd
          = this->db().lookupObject<IOdictionary>("transportProperties");
    scalar lR =  (new dimensionedScalar(IOd.lookup("lR")))->value();
    scalarField grC = -C.boundaryField()[patchID].snGrad();
    vectorField pNf = mesh.boundaryMesh()[patchID].faceNormals();
    
    const scalar dt = this->db().time().deltaTValue();
    faceDispl = dt * lR * grC * pNf;
    
    valuePointPatchField<vector>::updateCoeffs();
}
 
  namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        normalMotionSlipPointPatchVectorField
    );
}


// ************************************************************************* //
