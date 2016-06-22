/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
  
\*---------------------------------------------------------------------------*/

#include "nonLinearFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nonLinearFvPatchScalarField::nonLinearFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    l_T(1.0),
    Cth(1.0),
    n1(1.0),
    n2(1.0)
{
  if(debug) {
      Info << "nonLinearFvPatchScalarField::nonLinearFvPatchScalarField 1" << endl;
  }

  this->refValue() = pTraits<scalar>::zero;
  this->refGrad() = pTraits<scalar>::zero;
  this->valueFraction() = 1;
}

nonLinearFvPatchScalarField::nonLinearFvPatchScalarField
(
    const nonLinearFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    l_T(ptf.l_T),
    Cth(ptf.Cth),
    n1(ptf.n1),
    n2(ptf.n2)
{
    if(debug) {
        Info << "nonLinearFvPatchScalarField::nonLinearFvPatchScalarField 2" << endl;
    }
}

nonLinearFvPatchScalarField::nonLinearFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF)
{
  if(debug) {
    Info << "nonLinearFvPatchScalarField::nonLinearFvPatchScalarField 3" << endl;
  }

  l_T = 1.0;
  
  if (!dict.readIfPresent<scalar>("Cth", Cth))
  {
    Cth = 1.0;
    WarningIn(
        "nonLinearFvPatchScalarField::nonLinearFvPatchScalarField"
        "("
        "const fvPatch& p,"
        "const DimensionedField<Type, volMesh>& iF,"
        "const dictionary& dict"
        ")"
    ) << "No value defined for Cth"
        << " on " << this->patch().name() << " therefore using "
        << Cth
        << endl;
  }
  
  if (!dict.readIfPresent<scalar>("n1", n1))
  {
    n1 = 1.0;
    WarningIn(
        "nonLinearFvPatchScalarField::nonLinearFvPatchScalarField"
        "("
        "const fvPatch& p,"
        "const DimensionedField<Type, volMesh>& iF,"
        "const dictionary& dict"
        ")"
    ) << "No value defined for n1"
        << " on " << this->patch().name() << " therefore using "
        << n1
        << endl;
  }
  
  if (!dict.readIfPresent<scalar>("n2", n2))
  {
    n2 = 1.0;
    WarningIn(
        "nonLinearFvPatchScalarField::nonLinearFvPatchScalarField"
        "("
        "const fvPatch& p,"
        "const DimensionedField<Type, volMesh>& iF,"
        "const dictionary& dict"
        ")"
    ) << "No value defined for n2"
        << " on " << this->patch().name() << " therefore using "
        << n2
        << endl;
  }

  if(debug) {
    Info << "Cth: "<< Cth << endl;
    Info << "n1: "<< n1 << endl;
    Info << "n2: "<< n2 << endl;
  }
  
  if (dict.found("value"))
  {
    fvPatchField<scalar>::operator=
    (
      scalarField("value", dict, p.size())
    );
  }
  else
  {
    fvPatchField<scalar>::operator=(this->refValue());
    WarningIn(
        "nonLinearFvPatchScalarField<Type>::nonLinearFvPatchScalarField"
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
      Info << "nonLinearFvPatchScalarField<Type>::nonLinearFvPatchScalarField 2  "
              "refValue set" << endl;
  }
  
  if (dict.found("refGradient")) {
      this->refGrad() = Field<scalar>("refGradient", dict, p.size());
  } else {
      this->refGrad() = pTraits<scalar>::zero;
  }

  if (dict.found("valueFraction")) {
      this->valueFraction() = Field<scalar>("valueFraction", dict, p.size());
  } else {
      this->valueFraction() = 1;
  }
  
  if (!this->updated())
  {
    this->mixedFvPatchScalarField::updateCoeffs();
  }

  Field<scalar>::operator=
                (
                  this->valueFraction()*this->refValue()
                  +
                  (1.0 - this->valueFraction())*
                  (
                    this->patchInternalField()
                    + this->refGrad()/this->patch().deltaCoeffs()
                  )
                );
  
  fvPatchField<scalar>::evaluate();
}


nonLinearFvPatchScalarField::nonLinearFvPatchScalarField
(
    const nonLinearFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    l_T(ptf.l_T),
    Cth(ptf.Cth),
    n1(ptf.n1),
    n2(ptf.n2)
{
  if(debug) {
      Info << "nonLinearFvPatchScalarField<Type>::nonLinearFvPatchScalarField 4" << endl;
  }
}


nonLinearFvPatchScalarField::nonLinearFvPatchScalarField
(
    const nonLinearFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    l_T(ptf.l_T),
    Cth(ptf.Cth),
    n1(ptf.n1),
    n2(ptf.n2)
{
    if(debug) {
        Info << "nonLinearFvPatchScalarField<Type>::nonLinearFvPatchScalarField 5" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void nonLinearFvPatchScalarField::updateCoeffs()
{
  if(debug) {
      Info << "nonLinearFvPatchScalarField<Type>::updateCoeffs" << endl;
  }
  if (this->updated())
  {
      return;
  }

  if(debug) {
      Info << "nonLinearFvPatchScalarField<Type>::updateCoeffs - updating" << endl;
  }

  const IOdictionary& iod = this->db().objectRegistry::template 
                            lookupObject<IOdictionary>("transportProperties");
  if( !iod.readIfPresent<scalar>("l_T", l_T) ){
    SeriousErrorIn("nonLinearFvPatchScalarField<Type>::nonLinearFvPatchScalarField")
            <<"There is no l_T parameter in transportProperties dictionary"
            <<exit(FatalError);
  }

  if(debug) {
      Info << "l_T: "<< l_T << endl;
  }

  Field<scalar>& val = *this;

  scalarField del = mag( 1. / this->patch().deltaCoeffs() );

  scalarField fb = mag(val);

  scalarField f( fb.size(), 0.0 );
  scalarField refVal( fb.size(), 0.0 );
  scalarField refGr( fb.size(), 0.0 );

  scalar n1_ = n1-1;
  scalar n2_ = n2-1;

  scalar A1         = 1.0;
  scalar A2         = A1 * pow(Cth, n1_-n2_);
  scalar gamma      = 0.01;

  forAll(fb, i){
    scalar ff   = fb[i];
    scalar c1   = A1 * pow(ff, n1);
    scalar c2   = A2 * pow(ff, n2);
    scalar dc1  = A1 * n1 * pow(ff, n1_);
    scalar dc2  = A2 * n2 * pow(ff, n2_);
    scalar w    = 0.5 + 0.5 * tanh( (ff-Cth)/gamma );
    scalar dw   = 0.5 / gamma / pow( cosh( (ff-Cth)/gamma ), 2 );

    scalar R    = w * c1 + (1-w) * c2;
    scalar dR   = dw * (c1-c2) + w * (dc1-dc2) + dc2;

    scalar alpha  = del[i] / l_T;
    scalar adR    = alpha * dR;

    f[i]      = adR / (1 + adR );
    refVal[i] = ff;
    refGr[i]  = - R / l_T;
  }

  this->refValue() = pTraits<scalar>::one * refVal;
  this->refGrad() = pTraits<scalar>::one * refGr;
  this->valueFraction() = f;

  if(debug) {
    Field<scalar> iF = this->patchInternalField();
    Field<scalar> newF = (1.0 - this->valueFraction()) * iF;
    scalar residual = ( max( mag( newF ) )<SMALL ) ? (1.0) : (max( mag( newF-val ) ) / max( mag( newF ) )) ;
    Info << "  Nonlinear boundary field residual " << residual << nl << nl;
  }

  mixedFvPatchScalarField::updateCoeffs();
}

void nonLinearFvPatchScalarField::write(Ostream& os) const
{
  if(debug) {
    Info << "nonLinearFvPatchScalarField<Type>::write" << endl;
  }
  mixedFvPatchScalarField::write(os);
  os.writeKeyword("Cth")<< Cth << token::END_STATEMENT << nl;
  os.writeKeyword("n1")<< n1 << token::END_STATEMENT << nl;
  os.writeKeyword("n2")<< n2 << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
