scalarField delInv      = this->patch().deltaCoeffs();

const polyMesh& mesh    = this->internalField().mesh();
const label& patchID    = this->patch().index();
vectorField faceCent    = mesh.boundaryMesh()[patchID].faceCentres();
vectorField faceNorm    = mesh.boundaryMesh()[patchID].faceNormals();

const IOdictionary& iod = this->db().lookupObject<IOdictionary>("transportProperties");
scalar G                = (new dimensionedScalar(iod.lookup("G")))->value();
scalarField alpha       = -G * faceNorm.component(2) / delInv;

this->refValue()        = faceCent.component(2) - 1;
this->refGrad()         = pTraits<scalar>::zero;
this->valueFraction()   = alpha / (1 + alpha);

