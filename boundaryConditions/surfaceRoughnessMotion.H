/////////////////////////////////////////////
// parameters for RoughnessGenerator class //
/////////////////////////////////////////////

int majDir =              0;
int minDir =              1;
scalar majLen =         380;
scalar minLen =         330;
int majNum =            380;
int minNum =            330;

int seed =           513287;
scalar roughness =     0.05;
scalar dHurst =           0;
scalar cutLen =           2;
scalar maxDisp =        0.1;

/////////////////////////////////////////////

RoughnessGenerator rg( seed, majNum, minNum, majLen, minLen, roughness, dHurst,
        cutLen, maxDisp);

const polyMesh& mesh = this->internalField().mesh()();
// Get patch ID for moving boundary
const label& patchID = this->patch().index();

vectorField faceNorm = mesh.boundaryMesh()[patchID].faceNormals();

scalarField fftDisplacement(faceNorm.size(), 0.0);

scalarField faceDisp(faceNorm.size(), 0.0);

rg.getFFTdisp(fftDisplacement);

pointField pointFace = mesh.boundaryMesh()[patchID].faceCentres();
scalar maxMaj = max( pointFace.component(majDir) );
scalar maxLat = max( pointFace.component(minDir) );
scalar minMaj = min( pointFace.component(majDir) );
scalar minLat = min( pointFace.component(minDir) );
double Llat = maxLat - minLat;
double Lmaj = maxMaj - minMaj;

forAll(pointFace, i)
{
    scalar curMaj = pointFace[i].component(majDir) - minMaj;
    scalar curLat = pointFace[i].component(minDir) - minLat;

    int curm = std::floor(curMaj / Lmaj * (majNum-1));
    int curn = std::floor(curLat / Llat * (minNum-1));
    int ind = curn + minNum * curm;

    faceDisp[i] = fftDisplacement[ind];
}

Info <<  "Maximum and minimum face displacements  " 
       <<  max(fftDisplacement)  <<  "  " << min(fftDisplacement) << endl;

forAll( faceDisp, i )
{
    faceDisp[i] = std::min(faceDisp[i],  maxDisp);
    faceDisp[i] = std::max(faceDisp[i], -maxDisp);
}

Info <<  "Maximum and minimum face displacements  " 
       <<  max(faceDisp)  <<  "  " << min(faceDisp) << endl;


vectorField faceVec(faceNorm.size(), vector::zero);

forAll( faceVec, i )
    faceVec[i] = faceDisp[i] * faceNorm[i];

this->setDisp(faceVec);
