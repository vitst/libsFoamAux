//  Coded Dirichlet condition //

const label&    patchID     = this->patch().index();
const polyMesh& mesh        = this->internalField().mesh();
vectorField     faceCenters = mesh.boundaryMesh()[patchID].faceCentres();
scalarField     refVal      = faceCenters.component(2);

this->refValue()      = refVal;
this->refGrad()       = pTraits<scalar>::zero;
this->valueFraction() = pTraits<scalar>::one;

