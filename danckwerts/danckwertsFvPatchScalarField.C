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
#include "surfaceFields.H"
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
  mixedFvPatchScalarField(p, iF)
{
  if(debug) {
      Info << "danckwertsFvPatchField<Type>::danckwertsFvPatchField 1" << endl;
  }

  this->refValue() = pTraits<scalar>::zero;
  this->refGrad() = pTraits<scalar>::zero;
  this->valueFraction() = 1;
  
  vf_ = 1;
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
  mixedFvPatchScalarField(ptf, p, iF, mapper),
        vf_(1.0)
{
  if(debug) {
      Info << "danckwertsFvPatchField<Type>::danckwertsFvPatchField 1.1" << endl;
  }
}

Foam::danckwertsFvPatchScalarField::
danckwertsFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
  mixedFvPatchScalarField(p, iF)
{
  if(debug) {
      Info << "danckwertsFvPatchField<Type>::danckwertsFvPatchField 2" << endl;
  }
  
  
  if (dict.found("value"))
  {
    fvPatchScalarField::operator=
    (
      scalarField("value", dict, p.size())
    );
  }
  else
  {
    fvPatchScalarField::operator=(this->refValue());
    WarningIn(
        "!!!nonLinearFvPatchField<Type>::nonLinearFvPatchField"
        "("
        "const fvPatch& p,"
        "const DimensionedField<Type, volMesh>& iF,"
        "const dictionary& dict"
        ")"
    ) << "No value defined for " << this->dimensionedInternalField().name()
        << " on " << this->patch().name() << " therefore using "
        << this->refValue()
        << endl;
  }

  if(debug) {
      Info << "danckwertsFvPatchField<Type>::danckwertsFvPatchField 2  "
              "value set" << endl;
  }

  if (dict.found("refValue"))
    this->refValue() = scalarField("refValue", dict, p.size());
  else
  {
    if (dict.found("value")) {
      // make sure that refValue has a sensible value for the "update" below
      this->refValue() = scalarField(this->patchInternalField());
    }
    else{
      this->refValue() = pTraits<scalar>::zero;
    }
  }

  if(debug) {
      Info << "danckwertsFvPatchField<Type>::danckwertsFvPatchField 2  "
              "refValue set" << endl;
  }

  if (dict.found("refGradient")) {
      this->refGrad() = scalarField("refGradient", dict, p.size());
  } else {
      this->refGrad() = pTraits<scalar>::zero;
  }

  if(debug) {
      Info << "danckwertsFvPatchField<Type>::danckwertsFvPatchField 2  "
              "refGrad set" << endl;
  }

  if (dict.found("valueFraction")) {
      this->valueFraction() = Field<scalar>("valueFraction", dict, p.size());
  } else {
      this->valueFraction() = 1;
  }

  if(debug) {
      Info << "danckwertsFvPatchField<Type>::danckwertsFvPatchField 2  "
              "valueFraction set" << endl;
  }

  if (!this->updated())
  {
    this->mixedFvPatchScalarField::updateCoeffs();
  }
  
  vf_ = 1.0;
  
  /*
  scalarField aux = 
          this->valueFraction()*this->refValue()
          +
          (1.0 - this->valueFraction())*
          (
              this->patchInternalField()
              + this->refGrad()/this->patch().deltaCoeffs()
          );

  scalarField::operator=
      (
          aux
      );
   */

  if(debug) {
      Info << "danckwertsFvPatchField<Type>::danckwertsFvPatchField 2  "
              "operator" << endl;
  }
  
  fvPatchScalarField::evaluate();
}

Foam::danckwertsFvPatchScalarField::
danckwertsFvPatchScalarField
(
    const danckwertsFvPatchScalarField& ptf
)
:
  mixedFvPatchScalarField(ptf)
{
  if(debug) {
    Info << "danckwertsFvPatchField<Type>::danckwertsFvPatchField 3" << endl;
  }
}

Foam::danckwertsFvPatchScalarField::
danckwertsFvPatchScalarField
(
    const danckwertsFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF)
{
  if(debug) {
    Info << "danckwertsFvPatchField<Type>::danckwertsFvPatchField 4" << endl;
  }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::danckwertsFvPatchScalarField::updateCoeffs()
{
  if(debug) {
    Info<<"danckwerts::updateCoeffs:  updated: "<<updated()<<endl;
  }
  
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
    else if( this->db().find("D")()->type() == "surfaceScalarField" ){
      const scalarField& patchD =
            this->patch().template lookupPatchField<surfaceScalarField, scalar>("D");
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
  
  this->refValue() = pTraits<scalar>::one;
  this->refGrad() = pTraits<scalar>::zero;
  this->valueFraction() = AA / (AA+1.0);
  //this->valueFraction() = 1.0;
  
  vf_ = AA / (AA+1.0);
  
  mixedFvPatchScalarField::updateCoeffs();
}


void Foam::danckwertsFvPatchScalarField::evaluate(const Pstream::commsTypes)
{
  if(debug) {
      Info << "danckwertsFvPatchField<Type>::evaluate" << endl;
  }
  
  if (!this->updated()) this->updateCoeffs();

  scalarField iF = this->patchInternalField();

  scalarField aux = //(1-vf_) * iF;
          
    this->valueFraction()*this->refValue()
    +
    (1.0 - this->valueFraction())*
    (
      iF
      +
      this->refGrad()/this->patch().deltaCoeffs()
    );
           
  
  //Info<<nl<< "this->patchInternalField()  "<<this->patchInternalField()<<nl<<endl;
  //Info<<nl<< "*this  "<<*this<<nl<<endl;

  scalarField::operator=
  (
    aux
  );
  
  //Info<<nl<< "!!!!!AAAAA this->patchInternalField()  "<<this->patchInternalField()<<nl<<endl;
  //Info<<nl<< "!!!!!AAAAA *this  "<<*this<<nl<<endl;

  fvPatchScalarField::evaluate();
}




Foam::tmp<Foam::Field<Foam::scalar> > 
Foam::danckwertsFvPatchScalarField::snGrad() const
{
  if(debug) {
      Info << "danckwerts::snGrad()" << endl;
  }
    return
        this->valueFraction()
       *(this->refValue() - this->patchInternalField())
       *this->patch().deltaCoeffs()
      + (1.0 - this->valueFraction())*this->refGrad();
}

Foam::tmp<Foam::Field<Foam::scalar> > 
Foam::danckwertsFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
  if(debug) {
      Info << "danckwerts::valueInternalCoeffs" << endl;
  }
  
    return scalar(pTraits<scalar>::one)*(1.0 - this->valueFraction());
}

Foam::tmp<Foam::Field<Foam::scalar> > 
Foam::danckwertsFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
  if(debug) {
      Info << "danckwerts::valueBoundaryCoeffs" << endl;
  }
    return
         this->valueFraction() * this->refValue()
       + 
         (1.0 - this->valueFraction())
            *
         this->refGrad()/this->patch().deltaCoeffs();
}

Foam::tmp<Foam::Field<Foam::scalar> > 
Foam::danckwertsFvPatchScalarField::gradientInternalCoeffs() const
{
  if(debug) {
      Info << "danckwerts::gradientInternalCoeffs()" << endl;
  }
    return 
          - scalar(pTraits<scalar>::one)
            *
            this->valueFraction()
            *
            this->patch().deltaCoeffs();
}

Foam::tmp<Foam::Field<Foam::scalar> > 
Foam::danckwertsFvPatchScalarField::gradientBoundaryCoeffs() const
{
  if(debug) {
      Info << "danckwerts::gradientBoundaryCoeffs()" << endl;
  }
    return
        this->valueFraction() * this->patch().deltaCoeffs() * this->refValue()
      + (1.0 - this->valueFraction()) * this->refGrad();
}





void Foam::danckwertsFvPatchScalarField::write(Ostream& os) const
{
  //fvPatchScalarField::write(os);
  //writeEntry("value", os);
  mixedFvPatchScalarField::write(os);
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
