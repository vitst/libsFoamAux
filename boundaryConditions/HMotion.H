const label& patchID        = this->patch().index();
const volScalarField& field = this->db().objectRegistry::lookupObject<volScalarField>("H");
scalarField boundaryField   = field.boundaryField()[patchID];
scalarField gradField       = field.boundaryField()[patchID].snGrad();

const polyMesh& mesh = this->internalField().mesh()();
vectorField faceCent = mesh.boundaryMesh()[patchID].faceCentres();
vectorField faceNorm = mesh.boundaryMesh()[patchID].faceNormals();

scalarField delta    = -(faceCent.component(2) - boundaryField) / (faceNorm.component(2) - gradField);
vectorField fD       = faceNorm * delta;

this->setDisp(fD);
