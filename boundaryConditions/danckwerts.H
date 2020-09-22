//  Danckwerts boundary condition - read D from transport properties  //

const scalarField delInv = this->patch().deltaCoeffs();
const scalarField Un     = this->patch().template
                           lookupPatchField<volVectorField, vector>("U")
                         & (-1)*this->patch().nf();

const IOdictionary& iod = this->db().lookupObject<IOdictionary>("transportProperties");
scalar D = (new dimensionedScalar(iod.lookup("D")))->value();

this->refValue() = pTraits<scalar>::one;
this->refGrad()  = pTraits<scalar>::zero;
this->valueFraction() = Un / (D*delInv + Un);
