/*
  Here we read a scalar volume field C to get its gradiend
  and a scalar lR from transportProperties dictionary.
  Also, we calculate surface curvature.
*/
const polyMesh& mesh = this->internalField().mesh()();

const volScalarField& field =
            this->db().objectRegistry::lookupObject<volScalarField>("C");
const label& patchID = this->patch().index();
const IOdictionary& IOd
          = this->db().lookupObject<IOdictionary>("transportProperties");
scalar scalarVal =  (new dimensionedScalar("lR", dimLength, IOd))->value();

scalarField gradField = -field.boundaryField()[patchID].snGrad();
vectorField faceNorm = mesh.boundaryMesh()[patchID].faceNormals();

labelHashSet includePatches(1);
// insert only current patch
includePatches.insert(patchID);
labelList faceMap;
triSurface wallTriSurface
(
    triSurfaceTools::triangulate( mesh.boundaryMesh(), includePatches, faceMap)
);

scalarField curv = triSurfaceTools::curvatures(wallTriSurface);

Info<<"Curvature: "<<nl<<curv<<nl;
std::exit(0);

const scalar dt = this->db().time().deltaTValue();
vectorField fD = dt * scalarVal * gradField * faceNorm;

this->setDisp(fD);
this->setRlxON(true);
