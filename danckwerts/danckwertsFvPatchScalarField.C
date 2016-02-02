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


#include <typeinfo>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::danckwertsFvPatchScalarField::
danckwertsFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
  fixedValueFvPatchScalarField(p, iF)
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
  fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{
}

Foam::danckwertsFvPatchScalarField::danckwertsFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
  fixedValueFvPatchScalarField(p, iF)
{
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
    fixedValueFvPatchScalarField(ptf)
{
}


Foam::danckwertsFvPatchScalarField::
danckwertsFvPatchScalarField
(
    const danckwertsFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::danckwertsFvPatchScalarField::updateCoeffs()
{
  if (updated()) return;
  
  vectorField n   = (-1)*this->patch().nf();
  scalarField iF  = this->patchInternalField();
  vectorField del = (-1)*this->patch().delta();
  vectorField boundaryU =
          this->patch().template lookupPatchField<volVectorField, vector>("U");
  
  scalarField AA;
  if(this->db().find("D") != this->db().end()){
    if( this->db().find("D")()->type() == "volScalarField" ){
      const scalarField& patchD =
            this->patch().template lookupPatchField<volScalarField, scalar>("D");
      AA = (boundaryU & n) * (del & n) / patchD;
    }
    else if( this->db().find("D")()->type() == "volSphericalTensorField" ){
      const sphericalTensorField& patchD =
            this->patch().template lookupPatchField<volSphericalTensorField, sphericalTensor>("D");
      AA = (boundaryU & n) * (del & n) / (n & patchD & n);
    }
    else if( this->db().find("D")()->type() == "volSymmTensorField" ){
      const symmTensorField& patchD =
            this->patch().template lookupPatchField<volSymmTensorField, symmTensor>("D");
      AA = (boundaryU & n) * (del & n) / (n & patchD & n);
    }
    else if( this->db().find("D")()->type() == "volTensorField" ){
      const tensorField& patchD =
            this->patch().template lookupPatchField<volTensorField, tensor>("D");
      AA = (boundaryU & n) * (del & n) / (n & patchD & n);
    }
    else{
      SeriousErrorIn("danckwertsFvPatchScalarField::updateCoeffs()")
              <<"D type is not implemented. D type is "<< this->db().find("D")()->type()
              <<exit(FatalError);
    }
  }
  else{
    const IOdictionary& iod = this->db().lookupObject<IOdictionary>("transportProperties");
    scalar patchD = (new dimensionedScalar(iod.lookup("D")))->value();
    AA = (boundaryU & n) * (del & n) / patchD;
  }
  
  operator==( (AA+iF) / (AA+1.0) );
  
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
