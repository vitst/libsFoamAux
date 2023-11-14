//  parabolic boundary conditions for velocity   //

/*
const scalarField delInv = this->patch().deltaCoeffs();
const scalarField Un     = this->patch().template
                           lookupPatchField<volVectorField, vector>("U")
                         & (-1)*this->patch().nf();
*/
//const IOdictionary& iod = this->db().lookupObject<IOdictionary>("transportProperties");
//scalar D = (new dimensionedScalar("D", (dimArea/dimTime), iod))->value();

const polyMesh& mesh = this->internalField().mesh()();

const label& patchID = this->patch().index();
vectorField& Ub = *this;

vector n(1, 0, 0);
scalar averU = 2.0;
scalar minY = 1.0
scalar maxY = 150.0;

const pointField& fc = mesh.boundaryMesh()[patchID].faceCentres();

const scalarField yc = Foam::mag(fc.component(vector::Y));

scalarField newValues(fv.size(), 0.0);

newValues = 2.0 * averU * (1 - Foam::pow(2.0*(yc-minY)/(maxY-minY)-1.0,2))

this->refValue() = n * newValues;
this->refGrad()  = pTraits<vector>::zero;
this->valueFraction() = pTraits<scalar>::one;
