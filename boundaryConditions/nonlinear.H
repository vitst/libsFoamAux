//  Power-law kinetics  //

const IOdictionary& IOd                 //  Get auxiliary input
        = this->db().lookupObject<IOdictionary>("transportProperties");
scalar lR =  (new dimensionedScalar(IOd.lookup("lR")))->value();

scalar n1     = 1.0;                    //  Reaction parameters
scalar n2     = 1.0;
scalar Cth    = 0.0;
scalar gamma  = 0.01;

// Update coefficients

Field<scalar>& val = *this;
scalarField    del = mag( 1. / this->patch().deltaCoeffs() );
scalarField    field  = mag(val);
scalarField    frac( field.size(), 0.0 );
scalarField    refVal( field.size(), 0.0 );
scalarField    refGr( field.size(), 0.0 );

scalar Rth  = pow(Cth, n1-n2);
forAll(field, i)
{
  scalar fi   = field[i];
  scalar c1   = pow(fi, n1);
  scalar c2   = pow(fi, n2) * Rth;
  scalar dc1  = n1 * pow(fi, n1-1);
  scalar dc2  = n2 * pow(fi, n2-1) * Rth;

  scalar w    = 0.5 + 0.5 * tanh( (fi-Cth)/gamma );
  scalar dw   = 0.5 / gamma / pow( cosh( (fi-Cth)/gamma ), 2 );
  scalar R    = w * c1 + (1-w) * c2;
  scalar dR   = dw * (c1-c2) + w * (dc1-dc2) + dc2;

  scalar alpha  = del[i] / lR;
  scalar adR    = alpha * dR;
  frac[i]   = adR / ( 1 + adR );
  refVal[i] = fi;
  refGr[i]  = -R / lR;
}

this->refValue() = pTraits<scalar>::one * refVal;
this->refGrad()  = pTraits<scalar>::one * refGr;
this->valueFraction() = frac;
