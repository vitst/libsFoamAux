const label& patchID        = this->patch().index();
const volScalarField& field = this->db().objectRegistry::lookupObject<volScalarField>("C");
scalarField gradField       = field.boundaryField()[patchID].snGrad();

const scalar dt      = this->db().time().deltaTValue();
const polyMesh& mesh = this->internalField().mesh()();
vectorField faceNorm = mesh.boundaryMesh()[patchID].faceNormals();
vectorField fD       = -dt * gradField * faceNorm;

const IOdictionary& IOd
        = this->db().lookupObject<IOdictionary>("perchedFoamDict");
bool moveBoundary (Foam::readBool( IOd.lookup("moveBoundary")));
if (!moveBoundary) fD = fD * 0;

this->setDisp(fD);
