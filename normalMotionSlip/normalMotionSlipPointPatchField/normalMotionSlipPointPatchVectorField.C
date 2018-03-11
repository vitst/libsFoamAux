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
#include "IFstream.H"

#include "coupledPatchInterpolation.H"
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
{
    if (!dict.readIfPresent<word>("field", fieldName))
    {
        SeriousErrorIn("normalMotionSlipPointPatchVectorField")
            << "No field name parameter 'field' in "<< dict.name()
            << exit(FatalError);
    }
    
    if( debug )
    {
      Info << "field name: " << fieldName << nl;
    }
    
    if (!dict.readIfPresent<word>("scalarName", scalarName))
    {
        SeriousErrorIn("normalMotionSlipPointPatchVectorField")
            << "No parameter 'scalarName' in "<< dict.name()
            << exit(FatalError);
    }
    
    if( debug )
    {
      Info << "scalar name: " << scalarName << nl;
    }
}


Foam::normalMotionSlipPointPatchVectorField::normalMotionSlipPointPatchVectorField
(
    const normalMotionSlipPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    valuePointPatchField<vector>(ptf, p, iF, mapper),
    faceDispl( ptf.faceDispl ),
    fieldName(ptf.fieldName),
    scalarName(ptf.scalarName)
{}

Foam::normalMotionSlipPointPatchVectorField::normalMotionSlipPointPatchVectorField
(
    const normalMotionSlipPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    valuePointPatchField<vector>(ptf, iF),
    faceDispl( ptf.faceDispl ),
    fieldName(ptf.fieldName),
    scalarName(ptf.scalarName)
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
    
    if( debug )
    {
      Info << "Reading field: " << fieldName 
              << "  and scalar: "<< scalarName << nl;
    }

    const volScalarField& field =
            this->db().objectRegistry::lookupObject<volScalarField>(fieldName);
    const label& patchID = this->patch().index();
    const IOdictionary& IOd
          = this->db().lookupObject<IOdictionary>("transportProperties");
    scalar scalarVal =  (new dimensionedScalar(IOd.lookup(scalarName)))->value();
    scalarField gradField = -field.boundaryField()[patchID].snGrad();
    vectorField faceNorm = mesh.boundaryMesh()[patchID].faceNormals();
    
    const scalar dt = this->db().time().deltaTValue();
    faceDispl = dt * scalarVal * gradField * faceNorm;

    valuePointPatchField<vector>::updateCoeffs();
}


void Foam::normalMotionSlipPointPatchVectorField::write
(
  Ostream& os
) const
{
    valuePointPatchField<vector>::write(os);
    
    if( debug )
    {
      Info << "Writing to BC field name: " << fieldName 
              << "  scalar: "<< scalarName << nl;
    }

    os.writeKeyword("field")<< fieldName << token::END_STATEMENT << nl;
    os.writeKeyword("scalarName")<< scalarName << token::END_STATEMENT << nl;
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
