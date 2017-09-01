/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "normalMotionSlipPointPatchVectorField.H"
#include "volFields.H"
#include "pointPatchFields.H"

#include "coupledPatchInterpolation.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::normalMotionSlipPointPatchVectorField::normalMotionSlipPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
//    slipPointPatchField<vector>(p, iF)
//    fixedValuePointPatchField<vector>(p, iF)
    valuePointPatchField<vector>(p, iF)
//    calculatedPointPatchField<vector>(p, iF)
{}


Foam::normalMotionSlipPointPatchVectorField::normalMotionSlipPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
//    slipPointPatchField<vector>(p, iF, dict)
//    fixedValuePointPatchField<vector>(p, iF, dict)
    valuePointPatchField<vector>(p, iF, dict)
//    calculatedPointPatchField<vector>(p, iF, dict)
{}


Foam::normalMotionSlipPointPatchVectorField::normalMotionSlipPointPatchVectorField
(
    const normalMotionSlipPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
//    slipPointPatchField<vector>(ptf, p, iF, mapper)
//    fixedValuePointPatchField<vector>(ptf, p, iF, mapper)
    valuePointPatchField<vector>(ptf, p, iF, mapper)
//    calculatedPointPatchField<vector>(ptf, p, iF, mapper)
{}


Foam::normalMotionSlipPointPatchVectorField::normalMotionSlipPointPatchVectorField
(
    const normalMotionSlipPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
//    slipPointPatchField<vector>(ptf, iF)
//    fixedValuePointPatchField<vector>(ptf, iF)
    valuePointPatchField<vector>(ptf, iF)
//    calculatedPointPatchField<vector>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::normalMotionSlipPointPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
  //Info<<"normalMotionSlipPointPatchVectorField::evaluate"<<nl;
  
  const polyMesh& mesh = this->internalField().mesh()();
        label patchID = this->patch().index();
  const fvMesh& fvmesh_ = refCast<const fvMesh>(mesh);
  coupledPatchInterpolation patchInterpolator
  ( 
    mesh.boundaryMesh()[patchID], fvmesh_
  );
  
    Field<vector>& iF = const_cast<Field<vector>&>(this->primitiveField());

        const volVectorField& cmu =
          this->db().objectRegistry::lookupObject<volVectorField>("cellMotionU");
        vectorField cmuB = cmu.boundaryField()[patchID];
		vectorField pm = patchInterpolator.faceToPointInterpolate( cmuB );
        
        vectorField pif = this->patchInternalField();
  vectorField pNf = mesh.boundaryMesh()[patchID].faceNormals();
  
    vectorField tvalues =
        ( 
            pm
            /*
						pif
            + 
            transform(I - pNf*pNf, pm) 
             */
        );
    
    //Info<< "Set iniDisp in normalMotionSlipFvPatchVectorField"<<nl;
    //cmu.boundaryFieldRef()[patchID].setIniDisp(true);
    
	//	Info<<"EV poMo: "<<max(pif)<<"  f0 "<<max(pm)<<"    tval: "<<max(tvalues)<<nl;
    //this->operator==( this->primitiveField() );
    //this->operator==( pointMotion );
    //this->operator==( pm );
    //this->operator==( tvalues() );
    this->operator==( tvalues );
    //this->setInInternalField(iF, tvalues);
    
	//	Info<<"IIInternField: "<<max(this->patchInternalField())<<nl;
  
  valuePointPatchField<vector>::evaluate();
  //fixedValuePointPatchField<vector>::evaluate();
  //calculatedPointPatchField<vector>::evaluate();
}

void Foam::normalMotionSlipPointPatchVectorField::updateCoeffs()
{
  //Info<<"normalMotionSlipPointPatchVectorField::updateCoeffs"<<nl;
  //iniDisp = true;
  
  const polyMesh& mesh = this->internalField().mesh()();
  
  const volScalarField& C =
          this->db().objectRegistry::lookupObject<volScalarField>("C");
  label patchID = this->patch().index();
  const IOdictionary& IOd
        = this->db().lookupObject<IOdictionary>("transportProperties");
  scalar lR =  (new dimensionedScalar(IOd.lookup("lR")))->value();
  scalarField grC = -C.boundaryField()[patchID].snGrad();
  vectorField pNf = mesh.boundaryMesh()[patchID].faceNormals();

  faceDispl = lR * grC * pNf;
  valuePointPatchField<vector>::updateCoeffs();
}
/*
void Foam::normalMotionSlipPointPatchVectorField::updateCoeffs()
{
  Info<<"normalMotionSlip  updateCoeffsp()"<<nl;
 
  / *
    tmp<Field<vector>> tvalues =
        transform(I - n_*n_, this->patchInternalField());

    // Get internal field to insert values into
    Field<vector>& iF = const_cast<Field<vector>&>(this->primitiveField());

    this->setInInternalField(iF, tvalues());
   * /
  
  
        const volScalarField& C =
          this->db().objectRegistry::lookupObject<volScalarField>("C");
  
        label patchID = this->patch().index();
        
        //const vectorField nHat(this->patch().nf());
//        const vectorField nHat(mesh.boundaryMesh()[patchID].faceNormals());

        / *
        const IOdictionary& IOd
              = this->db().lookupObject<IOdictionary>("transportProperties");
        
        scalar lR =  (new dimensionedScalar(IOd.lookup("lR")))->value();
         * /
        scalar lR = 0.5;
//        scalarField grC = -C.boundaryField()[patchID].snGrad();

  const fvMesh& fvmesh_ = refCast<const fvMesh>(mesh);
  coupledPatchInterpolation patchInterpolator
  ( 
    mesh.boundaryMesh()[patchID], fvmesh_
  );

  //surfaceScalarField ssc = fvc::snGrad(C);
  
  scalarField pCf = -C.boundaryField()[patchID].snGrad();
//  scalarField pCf = -ssc.boundaryField()[patchID];
  vectorField pNf = mesh.boundaryMesh()[patchID].faceNormals();
  
  scalarField motionC    = patchInterpolator.faceToPointInterpolate(pCf);
  vectorField motionN    = patchInterpolator.faceToPointInterpolate(pNf);
  
  // normalize point normals to 1
  forAll(motionN, ii) motionN[ii]/=mag(motionN[ii]);
  
  // set the velocity for the point motion
  vectorField pointMotion = lR * motionC * motionN;
  //vectorField pointMotion = lR * pCf * pNf;

        
//        Info<<"AAAAAAAAA  maxvpointMo= "<<max(pointMotion)<<"   GGG  "<<max(C)<<nl;
//        Info<<"CCCCCCCCC  maxPrimitiv= "<<max(this->primitiveField())<<nl;

    //tmp<Field<vector>> tvalues =
    //    ( 
  //          lR * grC * nHat
	//					pointMotion
    //        + 
  //          transform(I - nHat*nHat, this->patchInternalField()) 
    //        transform(I - pNf*pNf, this->primitiveField()) 
    //    );
            
    //    transform(I - n_*n_, this->patchInternalField());

    // Get internal field to insert values into
    Field<vector>& iF = const_cast<Field<vector>&>(this->primitiveField());

    //this->setInInternalField(iF, tvalues());
   
		//vectorField pm = patchInterpolator.faceToPointInterpolate( this->primitiveField() );
		//vectorField pm = patchInterpolator.faceToPointInterpolate( this->patchInternalField() );
        const volVectorField& cmu =
          this->db().objectRegistry::lookupObject<volVectorField>("cellMotionU");
        vectorField cmuB = cmu.boundaryField()[patchID];
		vectorField pm = patchInterpolator.faceToPointInterpolate( cmuB );
		//vectorField pm = this->patchInternalField();
  
    //tmp<Field<vector>> tvalues =
    vectorField tvalues =
        ( 
  //          lR * grC * nHat
						pointMotion
            -
            transform(I - pNf*pNf, pointMotion) 
//            + 
  //          transform(I - nHat*nHat, this->patchInternalField()) 
//            transform(I - pNf*pNf, pm) 
        );
    
		Info<<"poMo: "<<max(pointMotion)<<"  f0 "<<max(pm)<<"    tval: "<<max(tvalues)<<nl;
    //this->operator==( this->primitiveField() );
    //this->operator==( pointMotion );
    //this->operator==( pm );
    //this->operator==( tvalues() );
    this->operator==( tvalues );
    this->setInInternalField(iF, tvalues);
    
		Info<<"IIInternField: "<<max(this->patchInternalField())<<nl;
        
valuePointPatchField<vector>::updateCoeffs();
//fixedValuePointPatchField<vector>::updateCoeffs();
//calculatedPointPatchField<vector>::updateCoeffs();
}
*/
 

/*
void Foam::normalMotionSlipPointPatchVectorField::updateCoeffs()
{
  if (this->updated())
  {
    return;
  }
  
  Info<<"SSSSSSSSSSSSSSSSSSSSS "<<nl;
  if( this->db().foundObject<IOdictionary>("transportProperties") )
  {
    const polyMesh& mesh = this->internalField().mesh()();
  
        const volScalarField& C =
          this->db().objectRegistry::lookupObject<volScalarField>("C");
  
        label patchID = this->patch().index();
        
        //const vectorField nHat(this->patch().nf());
        const vectorField nHat(mesh.boundaryMesh()[patchID].faceNormals());

        / *
        const IOdictionary& IOd
              = this->db().lookupObject<IOdictionary>("transportProperties");
        
        scalar lR =  (new dimensionedScalar(IOd.lookup("lR")))->value();
         * /
        scalar lR = 0.5;
        scalarField grC = -C.boundaryField()[patchID].snGrad();
        
        Info<<"AAAAAAAAAAAAAAAAAAAAAA  grC[0]= "<<grC[0]<<nl;

    tmp<Field<vector>> tvalues =
        ( 
            lR * grC * nHat
            + 
            transform(I - nHat*nHat, this->patchInternalField()) 
        );
            
    //    transform(I - n_*n_, this->patchInternalField());

    // Get internal field to insert values into
    Field<vector>& iF = const_cast<Field<vector>&>(this->primitiveField());

    this->setInInternalField(iF, tvalues());
  
    this->operator==( this->primitiveField() );
  }
  fixedValuePointPatchField<vector>::updateCoeffs();
  
}
*/


/*
void Foam::normalMotionSlipPointPatchVectorField::write(Ostream& os) const
{
    //slipPointPatchField<vector>::write(os);
    //calculatedPointPatchField<vector>::write(os);
    valuePointPatchField<vector>::write(os);
    this->writeEntry("value", os);
    //os.writeKeyword("value") << this->patchInternalField() << Foam::token::END_STATEMENT << '\n';
    //os.writeKeyword("n") << n_ << token::END_STATEMENT << nl;
}
 */
 

bool Foam::normalMotionSlipPointPatchVectorField::getIniDisp()
{
  return iniDisp;
}





namespace Foam
{
  makePointPatchTypeField
  (
    pointPatchVectorField,
    normalMotionSlipPointPatchVectorField
  );
}


// ************************************************************************* //
