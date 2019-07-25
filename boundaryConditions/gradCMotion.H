const polyMesh& mesh = this->internalField().mesh()();

const volScalarField& field =
            this->db().objectRegistry::lookupObject<volScalarField>("C");
const label& patchID = this->patch().index();
const IOdictionary& IOd
          = this->db().lookupObject<IOdictionary>("transportProperties");
scalar scalarVal =  (new dimensionedScalar(IOd.lookup("lR")))->value();
scalarField gradField = -field.boundaryField()[patchID].snGrad();
vectorField faceNorm = mesh.boundaryMesh()[patchID].faceNormals();
    
const scalar dt = this->db().time().deltaTValue();
vectorField fD = dt * scalarVal * gradField * faceNorm;

this->setDisp(fD);
this->setRlxON(true);
