scalar Cth =  0.0;
scalar n1  =  1.0;
scalar n2  =  1.0;

bool limitingPlane = false;
vector refNorm(0,0,1);
point  refPoint (0,0,-4.7);

const IOdictionary& iod = this->db().objectRegistry::template
                 lookupObject<IOdictionary>("transportProperties");
scalar l_T = 1.0;
if( !iod.readIfPresent<scalar>("l_T", l_T) )
{
  SeriousErrorIn("nonLinearFvPatchScalarField::updateCoeffs()")
      <<"There is no l_T parameter in transportProperties dictionary"
      <<exit(FatalError);
}

// ***************************************************************************
plane limitPlane(refPoint, refNorm);
const scalar dt = this->db().time().deltaTValue();
const pointField& localFcs = this->patch().patch().faceCentres();
const pointField& localFns = this->patch().patch().faceNormals();
// ***************************************************************************

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

forAll(fb, i)
{
  scalar ff   = fb[i];
  scalar c1   = A1 * pow(ff, n1);
  scalar c2   = A2 * pow(ff, n2);
  scalar dc1  = A1 * n1 * pow(ff, n1_);
  scalar dc2  = A2 * n2 * pow(ff, n2_);

  scalar w    = 0.5 + 0.5 * tanh( (ff-Cth)/gamma );
  scalar dw   = 0.5 / gamma / pow( cosh( (ff-Cth)/gamma ), 2 );

  scalar R    = w * c1 + (1-w) * c2;
  scalar dR   = dw * (c1-c2) + w * (dc1-dc2) + dc2;

  // ***************************************************************************
  scalar theta = 1.0;
  if( limitingPlane )
  {
    scalar displ = R * dt;
    scalar distToInters = limitPlane.normalIntersect(localFcs[i], localFns[i]);
    if((displ > distToInters))
    {
      theta = distToInters / displ;
    }
  }
  if(theta<SMALL) theta = 0.0;
  R  *= theta;
  dR *= theta;
  // ***************************************************************************

  scalar alpha  = del[i] / l_T;
  scalar adR    = alpha * dR;

  f[i]      = adR / (1 + adR );
  refVal[i] = ff;
  refGr[i]  = - R / l_T;
}

this->refValue() = pTraits<scalar>::one * refVal;
this->refGrad() = pTraits<scalar>::one * refGr;
this->valueFraction() = f;
