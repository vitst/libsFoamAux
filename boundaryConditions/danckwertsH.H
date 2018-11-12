//  Danckwerts boundary condition - read D from transport properties  //

const scalarField delInv    = this->patch().deltaCoeffs();
const volScalarField& field = this->db().objectRegistry::lookupObject<volScalarField>("H");
scalarField Un              = field.boundaryField()[patchID].snGrad();

const IOdictionary& iod = this->db().lookupObject<IOdictionary>("transportProperties");
scalar D = (new dimensionedScalar(iod.lookup("D")))->value();

this->refValue() = pTraits<scalar>::one;
this->refGrad()  = pTraits<scalar>::zero;
this->valueFraction() = Un / (D * delInv + Un);
