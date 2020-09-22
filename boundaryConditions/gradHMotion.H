const label& patchID        = this->patch().index();
const volScalarField& field = this->db().objectRegistry::lookupObject<volScalarField>("H");
scalarField gradField       = field.boundaryField()[patchID].snGrad();

const polyMesh& mesh = this->internalField().mesh()();
vectorField faceNorm = mesh.boundaryMesh()[patchID].faceNormals();

const IOdictionary& iod = this->db().lookupObject<IOdictionary>("transportProperties");
scalar R = (new dimensionedScalar(iod.lookup("R")))->value();
    
vectorField fD = 0.1 * faceNorm * (-gradField / R + faceNorm.component(2));

this->setDisp(fD);
