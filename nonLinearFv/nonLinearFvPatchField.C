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
  
Description
    This is the patch for nonlinear boundary concentration field.
  
Contributors/Copyright:
    Vitaliy Starchenko (2015) <vitaliy.starchenko@gmail.com>

\*---------------------------------------------------------------------------*/

#include "nonLinearFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
nonLinearFvPatchField<Type>::nonLinearFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    l_T(1.0),
    Cth(1.0),
    n1(1.0),
    n2(1.0)
{
    if(debug) {
        Info << "nonLinearFvPatchField<Type>::nonLinearFvPatchField 1" << endl;
    }

    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 1;
}


template<class Type>
nonLinearFvPatchField<Type>::nonLinearFvPatchField
(
    const nonLinearFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    l_T(ptf.l_T),
    Cth(ptf.Cth),
    n1(ptf.n1),
    n2(ptf.n2)
{
    if(debug) {
        Info << "nonLinearFvPatchField<Type>::nonLinearFvPatchField 2" << endl;
    }
}


template<class Type>
nonLinearFvPatchField<Type>::nonLinearFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF)
{
    if(debug) {
        Info << "nonLinearFvPatchField<Type>::nonLinearFvPatchField 3" << endl;
    }

    const IOdictionary& iod = this->db().objectRegistry::template 
                              lookupObject<IOdictionary>("transportProperties");
    
    if( !iod.readIfPresent<scalar>("l_T", l_T) ){
      SeriousErrorIn("nonLinearFvPatchField<Type>::nonLinearFvPatchField")
              <<"There is no l_T parameter in transportProperties dictionary"
              <<exit(FatalError);
    }
    Cth = dict.lookupOrDefault<scalar>("Cth", 1.0);
    n1 = dict.lookupOrDefault<scalar>("n1", 1.0);
    n2 = dict.lookupOrDefault<scalar>("n2", 1.0);
    
    if(debug) {
        Info << "l_T: "<< l_T << endl;
        Info << "Cth: "<< Cth << endl;
        Info << "n1: "<< n1 << endl;
        Info << "n2: "<< n2 << endl;
    }
    
    if (dict.found("refValue")) {
        this->refValue() = Field<Type>("refValue", dict, p.size());
    } else {
        this->refValue() = pTraits<Type>::zero;
    }

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
	if (!dict.found("refValue")) {
 	    // make sure that refValue has a sensible value for the "update" below
	    this->refValue() = Field<Type>("value", dict, p.size());
	}
    }
    else
    {
        fvPatchField<Type>::operator=(this->refValue());
        WarningIn(
            "nonLinearFvPatchField<Type>::nonLinearFvPatchField"
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

    if (dict.found("refGradient")) {
        this->refGrad() = Field<Type>("refGradient", dict, p.size());
    } else {
        this->refGrad() = pTraits<Type>::zero;
    }

    if (dict.found("valueFraction")) {
        this->valueFraction() = Field<scalar>("valueFraction", dict, p.size());
    } else {
        this->valueFraction() = 1;
    }

    // emulate mixedFvPatchField<Type>::evaluate, but avoid calling "our" updateCoeffs
    if (!this->updated())
    {
        this->mixedFvPatchField<Type>::updateCoeffs();
    }

    Field<Type>::operator=
        (
            this->valueFraction()*this->refValue()
            +
            (1.0 - this->valueFraction())*
            (
                this->patchInternalField()
                + this->refGrad()/this->patch().deltaCoeffs()
            )
        );

    fvPatchField<Type>::evaluate();
}


template<class Type>
nonLinearFvPatchField<Type>::nonLinearFvPatchField
(
    const nonLinearFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    l_T(ptf.l_T),
    Cth(ptf.Cth),
    n1(ptf.n1),
    n2(ptf.n2)
{
    if(debug) {
        Info << "nonLinearFvPatchField<Type>::nonLinearFvPatchField 4" << endl;
    }
}


template<class Type>
nonLinearFvPatchField<Type>::nonLinearFvPatchField
(
    const nonLinearFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    l_T(ptf.l_T),
    Cth(ptf.Cth),
    n1(ptf.n1),
    n2(ptf.n2)
{
    if(debug) {
        Info << "nonLinearFvPatchField<Type>::nonLinearFvPatchField 5" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void nonLinearFvPatchField<Type>::updateCoeffs()
{
    if(debug) {
        Info << "nonLinearFvPatchField<Type>::updateCoeffs" << endl;
        //Info << "Value: " << this->valueExpression_ << endl;
        //Info << "Gradient: " << this->gradientExpression_ << endl;
        //Info << "Fraction: " << this->fractionExpression_ << endl;
    }
    if (this->updated())
    {
        return;
    }

    if(debug) {
        Info << "nonLinearFvPatchField<Type>::updateCoeffs - updating" << endl;
    }
    
    Field<Type>& val = *this;
    
    scalarField del = mag( 1. / this->patch().deltaCoeffs() );
    
    scalarField fb = mag(val);
    
    scalar n1_ = n1-1;
    scalar n2_ = n2-1;
    
    scalar A1 = 1.0;
    scalar A2 = A1 * std::pow(Cth, n1_-n2_);
    
    forAll(fb, i){
      scalar ff = fb[i];
      if( ff < Cth ){
        fb[i] = A1 * std::pow(ff, n1_);
      }
      else{
        fb[i] = A2 * std::pow(ff, n2_);
      }
    }
    
    scalarField delval = mag( del * fb );
    scalarField f   = delval / (l_T + delval);
    
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = f;
    
    if(debug) {
      Field<Type> iF = this->patchInternalField();
      Field<Type> newF = (1.0 - this->valueFraction()) * iF;
      scalar residual = ( max( mag( newF ) )<SMALL ) ? (1.0) : (max( mag( newF-val ) ) / max( mag( newF ) )) ;
      Info << "  Nonlinear boundary field residual " << residual << nl << nl;
    }

    mixedFvPatchField<Type>::updateCoeffs();
}

template<class Type>
void nonLinearFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    
    Field<Type> iF = this->patchInternalField();

    Field<Type>::operator=
    (
        (1.0 - this->valueFraction()) * iF
    );
    fvPatchField<Type>::evaluate();
}


template<class Type>
void nonLinearFvPatchField<Type>::write(Ostream& os) const
{
    if(debug) {
        Info << "nonLinearFvPatchField<Type>::write" << endl;
    }
    mixedFvPatchField<Type>::write(os);
    os.writeKeyword("Cth")<< Cth << token::END_STATEMENT << nl;
    os.writeKeyword("n1")<< n1 << token::END_STATEMENT << nl;
    os.writeKeyword("n2")<< n2 << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
