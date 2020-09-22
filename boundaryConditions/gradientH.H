//  Coded Dirichlet condition //

const label&    patchID  = this->patch().index();
const polyMesh& mesh     = this->internalField().mesh();
vectorField     faceNorm = mesh.boundaryMesh()[patchID].faceNormals();
scalarField     refGrad  = faceNorm.component(2);

const IOdictionary& iod  = this->db().lookupObject<IOdictionary>("transportProperties");
scalar R = (new dimensionedScalar(iod.lookup("R")))->value();

this->refValue()         = pTraits<scalar>::one;
this->refGrad()          = R * refGrad;
this->valueFraction()    = pTraits<scalar>::zero;

