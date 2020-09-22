/////////////////////////////////////////////
// parameters for RoughnessGenerator class //
/////////////////////////////////////////////

int majDir =              0;            // Directions: x,y,z
int minDir =              1;
int majNum =            200;            // Number of cells
int minNum =            200;
scalar majLen =         200;            // Dimensions
scalar minLen =         200;

int seed   =         199881;
scalar roughness =      1.0;
scalar dHurst =         1.0;
scalar cutoff =          10;
scalar maxDisp =       10.0;

/////////////////////////////////////////////

RoughnessGenerator rg(seed, majNum, minNum, roughness, dHurst, cutoff);

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
scalar majDel = majLen/majNum;
scalar minDel = minLen/minNum;
Info << "Dels: " << majDel << "  " << minDel << endl;

forAll(pointFace, i)
{
    scalar majPf = (pointFace[i].component(majDir) - majMin)/majDel;
    scalar minPf = (pointFace[i].component(minDir) - minMin)/minDel;

    int mPf = std::floor(majPf);
    int nPf = std::floor(minPf);
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
