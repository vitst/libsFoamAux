//  Danckwerts boundary condition - read D from transport properties  //

const IOdictionary& iod =               // Get auxiliary data 
    this->db().lookupObject<IOdictionary>("transportProperties");
scalar D = (new dimensionedScalar(iod.lookup("D")))->value();
const vectorField boundaryU = \
    this->patch().template lookupPatchField<volVectorField, vector>("U");
vectorField n   = (-1)*this->patch().nf();
vectorField del = (-1)*this->patch().delta();

scalarField alpha = (boundaryU & n) * (del & n) / D;

this->refValue() = pTraits<scalar>::one;
this->refGrad() = pTraits<scalar>::zero;
this->valueFraction() = alpha / (1 + alpha);
