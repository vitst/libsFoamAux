//  Linear kinetics with limiting plane  //
//  Add "#include plane.H" in boundary specification
//  Edit plane: refNorm and refPoint (below)

const IOdictionary& IOd                     //  Get auxiliary data
        = this->db().lookupObject<IOdictionary>("transportProperties");
scalar lR =  (new dimensionedScalar(IOd.lookup("lR")))->value();

bool limitingPlane = true;                  //  Limiting plane
vector refNorm(0,0,1);
point  refPoint (0,0,-4.5);
plane limitPlane(refPoint, refNorm);
const scalar dt = this->db().time().deltaTValue();
const pointField& localFcs = this->patch().patch().faceCentres();
const pointField& localFns = this->patch().patch().faceNormals();

//  Update coefficients

Field<scalar>& val = *this;
scalarField    del = mag( 1. / this->patch().deltaCoeffs() );
scalarField    field  = mag(val);
scalarField    frac( field.size(), 0.0 );
scalarField    refVal( field.size(), 0.0 );
scalarField    refGr( field.size(), 0.0 );

forAll(field, i)
{
  scalar R    = field[i];
  scalar dR   = 1.0;

  scalar theta = 1.0;
  if( limitingPlane )                   // Apply limiting plane
  {
    scalar displ = R * dt;
    scalar distToInters = limitPlane.normalIntersect(localFcs[i], localFns[i]);
    if (displ > distToInters) theta = distToInters / displ;
  }
  if (theta < SMALL) theta = 0.0;
  R  *= theta;
  dR *= theta;

  scalar alpha  = del[i] / lR;
  scalar adR    = alpha * dR;
  frac[i]   = adR / ( 1 + adR );
  refVal[i] = R;
  refGr[i]  = -R / lR;
}

this->refValue() = pTraits<scalar>::one * refVal;
this->refGrad()  = pTraits<scalar>::one * refGr;
this->valueFraction() = frac;
