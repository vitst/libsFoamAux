//  Danckwerts boundary condition - read D from transport properties  //

const vectorField boundaryU = this->patch().template lookupPatchField<volVectorField, vector>("U");
const IOdictionary& iod = this->db().lookupObject<IOdictionary>("transportProperties");
vectorField n   = (-1)*this->patch().nf();
vectorField del = (-1)*this->patch().delta();

scalar D = (new dimensionedScalar(iod.lookup("D")))->value();

scalarField alpha = (boundaryU & n) * (del & n) / D;

this->refValue() = pTraits<scalar>::one;
this->refGrad() = pTraits<scalar>::zero;
this->valueFraction() = alpha / (1 + alpha);
