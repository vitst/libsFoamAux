/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "danckwertsFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::danckwertsFvPatchScalarField::
danckwertsFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
  fixedValueFvPatchScalarField(p, iF),
  D_(1.0)
{
}

Foam::danckwertsFvPatchScalarField::
danckwertsFvPatchScalarField
(
    const danckwertsFvPatchScalarField& ptf, 
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
  fixedValueFvPatchScalarField(ptf, p, iF, mapper),
  D_(ptf.D_)
{
}

Foam::danckwertsFvPatchScalarField::
danckwertsFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
  fixedValueFvPatchScalarField(p, iF)
{
  D_ = 1.0; //default value (later we read it from the file)

  if (dict.found("value"))
  {
      fvPatchScalarField::operator=
      (
          scalarField("value", dict, p.size())
      );
  }
}

Foam::danckwertsFvPatchScalarField::
danckwertsFvPatchScalarField
(
    const danckwertsFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    D_(ptf.D_)
{
}


Foam::danckwertsFvPatchScalarField::
danckwertsFvPatchScalarField
(
    const danckwertsFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    D_(ptf.D_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::danckwertsFvPatchScalarField::updateCoeffs()
{
  if (updated()) return;
  
  const IOdictionary& iod = this->db().lookupObject<IOdictionary>("transportProperties");
  D_ = (new dimensionedScalar(iod.lookup("D")))->value();
  
  scalarField newValues(this->patchInternalField());

  vectorField del = this->patch().delta();

  const vectorField& boundaryU =
          this->patch().template lookupPatchField<volVectorField, vector>("U");
  
  scalar Dinv_ = 1.0 / D_;

  forAll(newValues, ii){
    scalar aa = mag( boundaryU[ii].z()*del[ii].z() ) * Dinv_;
    newValues[ii] = ( aa+newValues[ii] )/(aa+1.0);
  }
  
  operator==(newValues);

  fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::danckwertsFvPatchScalarField::write(Ostream& os) const
{
  fvPatchScalarField::write(os);
  writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  makePatchTypeField
  (
    fvPatchScalarField,
    danckwertsFvPatchScalarField
  );
}

// ************************************************************************* //
