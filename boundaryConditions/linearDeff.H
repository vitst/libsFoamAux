//  Linear kinetics  //

const IOdictionary& IOd                     //  Get auxiliary data
        = this->db().lookupObject<IOdictionary>("transportProperties");
scalar lR =  (new dimensionedScalar(IOd.lookup("lR")))->value();
scalar D  =  (new dimensionedScalar(IOd.lookup("D")))->value();
const volScalarField& auxField =
        this->db().objectRegistry::lookupObject<volScalarField>("Deff");
const label& patchID = this->patch().index();
scalarField patchD = auxField.boundaryField()[patchID] / D;

//  Update coefficients

Field<scalar>& val = *this;
scalarField    del = mag( 1. / this->patch().deltaCoeffs() );
scalarField    field  = mag(val);
scalarField    frac( field.size(), 0.0 );
scalarField    refVal( field.size(), 0.0 );
scalarField    refGr( field.size(), 0.0 );

forAll(field, i)
{
  scalar alpha = del[i] / lR / patchD[i];
  frac[i]   =  alpha / ( 1 + alpha );
}

this->refValue() = pTraits<scalar>::one * refVal;
this->refGrad()  = pTraits<scalar>::one * refGr;
this->valueFraction() = frac;
