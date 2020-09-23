//  Square nonlinear kinetics  //

Field<scalar>& val = *this;
scalarField    del = mag( 1. / this->patch().deltaCoeffs() );
scalarField    field  = mag(val);
scalarField    frac( field.size(), 0.0 );
scalarField    refVal( field.size(), 0.0 );
scalarField    refGr( field.size(), 0.0 );

const IOdictionary& IOd = this->db().lookupObject<IOdictionary>("transportProperties");
scalar lR =  (new dimensionedScalar("lR", dimLength, IOd))->value();

forAll(field, i)
{
  //scalar R    = field[i];
  scalar R    = field[i] * field[i];
  //scalar dR   = 1.0;
  scalar dR   = 2 * field[i];

  scalar alpha  = del[i] / lR;
  scalar adR    = alpha * dR;
  frac[i]   = adR / ( 1 + adR );
  refVal[i] = field[i]; // reference concentration is a boundary value from previous iteration
  refGr[i]  = -R / lR;
}

this->refValue() = pTraits<scalar>::one * refVal;
this->refGrad()  = pTraits<scalar>::one * refGr;
this->valueFraction() = frac;
