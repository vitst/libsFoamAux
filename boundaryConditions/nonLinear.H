//  Power-law kinetics with limiting plane  //

scalar n1     = 1.0;                // Reaction parameters
scalar n2     = 1.0;
scalar Cth    = 0.0;
scalar gamma  = 0.01;

bool limitingPlane = false;         // Limiting plane
vector refNorm(0,0,1);
point  refPoint (0,0,-4.7);
plane limitPlane(refPoint, refNorm);

//  Get auxiliary input data

const scalar dt = this->db().time().deltaTValue();
const pointField& localFcs = this->patch().patch().faceCentres();
const pointField& localFns = this->patch().patch().faceNormals();
const IOdictionary& IOd = this->db().objectRegistry::template
                 lookupObject<IOdictionary>("transportProperties");
scalar l_T = 1.0;
if( !IOd.readIfPresent<scalar>("l_T", l_T) )
{
  SeriousErrorIn("nonLinearFvPatchScalarField::updateCoeffs()")
      << "There is no l_T parameter in transportProperties dictionary"
      << exit(FatalError);
}

//  Update coefficients

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

  scalar theta = 1.0;
  if( limitingPlane )               // Apply limiting plane
  {
    scalar displ = R * dt;
    scalar distToInters = limitPlane.normalIntersect(localFcs[i], localFns[i]);
    if (displ > distToInters) theta = distToInters / displ;
  }
  if (theta < SMALL) theta = 0.0;
  R  *= theta;
  dR *= theta;

  scalar alpha  = del[i] / l_T;
  scalar adR    = alpha * dR;
  frac[i]   = adR / ( 1 + adR );
  refVal[i] = fi;
  refGr[i]  = -R / l_T;
}

this->refValue() = pTraits<scalar>::one * refVal;
this->refGrad()  = pTraits<scalar>::one * refGr;
this->valueFraction() = frac;
