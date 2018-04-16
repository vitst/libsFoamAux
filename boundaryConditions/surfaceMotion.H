/////////////////////////////////////////////
// parameters for RoughnessGenerator class //
/////////////////////////////////////////////

int majDir =              0;
int minDir =              1;
scalar majLen =        1000;
scalar minLen =         500;
int majNum =           1000;
int minNum =            500;

int seed =           513287;
scalar roughness =      0.1;
scalar dHurst =         0.5;
scalar cutLen =           4;
scalar maxDisp =        0.2;

/////////////////////////////////////////////

RoughnessGenerator rg( seed, majNum, minNum, majLen, minLen,
        roughness, dHurst, cutLen, maxDisp);

const polyMesh& mesh = this->internalField().mesh()();
// Get patch ID for moving boundary
const label& patchID = this->patch().index();

vectorField faceNorm = mesh.boundaryMesh()[patchID].faceNormals();

scalarField fftDisplacement(faceNorm.size(), 0.0);

scalarField faceDisp(faceNorm.size(), 0.0);

rg.getFFTdisp(fftDisplacement);

pointField pointFace = mesh.boundaryMesh()[patchID].faceCentres();
scalar majMin = min( pointFace.component(majDir) );
scalar minMin = min( pointFace.component(minDir) );

forAll(pointFace, i)
{
    scalar majPf = pointFace[i].component(majDir) - majMin;
    scalar minPf = pointFace[i].component(minDir) - minMin;

    int mPf = std::floor(majPf * majNum / majLen);
    int nPf = std::floor(minPf * minNum / minLen);
    int ind = mPf * minNum + nPf;

    faceDisp[i] = fftDisplacement[ind];
}
Info << "Maximum and minimum face displacements from FFT " 
     <<  max(faceDisp)  <<  "  " << min(faceDisp) << endl;

forAll( faceDisp, i )
{
    faceDisp[i] = std::min(faceDisp[i],  maxDisp);
    faceDisp[i] = std::max(faceDisp[i], -maxDisp);
}
Info << "Displacements trucated at +/- " << maxDisp << endl;

scalar avD  = average(faceDisp);
scalar avD2 = average(faceDisp*faceDisp);

Info << "Mean displacement " << abs(avD) << nl;
Info << "Roughness         " << sqrt(avD2 - avD*avD) << nl << endl;

vectorField faceVec(faceNorm.size(), vector::zero);

forAll( faceVec, i )
    faceVec[i] = faceDisp[i] * faceNorm[i];

this->setDisp(faceVec);
