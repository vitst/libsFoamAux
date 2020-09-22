/////////////////////////////////////////////
// parameters for individual Fourier modes //
/////////////////////////////////////////////

int majDir =              0;
int minDir =              1;
int dispDr =              2;
int majK   =              4;
int minK   =              4;
scalar majLen =         200;
scalar minLen =         200;

scalar roughness =        1;

/////////////////////////////////////////////

const polyMesh& mesh = this->internalField().mesh()();
const label&    patchID = this->patch().index();
pointField      pointFace = mesh.boundaryMesh()[patchID].faceCentres();
vectorField     faceVec(pointFace.size(), vector::zero);

scalar k1 = majK/majLen;
scalar k2 = minK/minLen;

forAll(pointFace, i)
{
    scalar kx = k1*pointFace[i].component(majDir);
    scalar ky = k2*pointFace[i].component(minDir);
    faceVec[i].component(dispDr) = roughness*std::exp(-kx*kx-ky*ky);
}

this->setDisp(faceVec);
