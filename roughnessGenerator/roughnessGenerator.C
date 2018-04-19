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

#include <complex>
#include <vector>
#include <fftw3.h>

#include "roughnessGenerator.H"
#include "Random.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label
Foam::RoughnessGenerator::index
(
    const label i,
    const label j
)
{
    return j + minNum * i;
}

Foam::scalar
Foam::RoughnessGenerator::power(double ksq)
{
    if (ksq == 0)    return 0;        // <rad^2> ~ 1/ksq^(1+H)
    if (ksq >  0.5)  return 0;        // cutoff wavelength = cutLen
    scalar p = Foam::pow(ksq, -(dHurst+1) );
    return std::sqrt(p);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RoughnessGenerator::RoughnessGenerator
(
    int seed_,
    int majNum_,
    int minNum_,
    double majLen_,
    double minLen_,
    double rgh_,
    double dHurst_,
    double cutLen_,
    double maxDisp_
)
:
    seed(seed_),
    majNum(majNum_),
    minNum(minNum_),
    majLen(majLen_),
    minLen(minLen_),
    rgh(rgh_),
    dHurst(dHurst_),
    cutLen(cutLen_),
    maxDisp(maxDisp_)
{
}


// ************************************************************************* //

void Foam::RoughnessGenerator::getFFTdisp(scalarField& disp)
{
    unsigned int MN = majNum*minNum;
    std::vector<std::complex<double> > f, F;
    f.resize(MN);
    F.resize(MN);

    Random rnd( seed );
    scalar TwoPi = constant::mathematical::twoPi;
    
    Info <<  "majLen:        "  <<  majLen      <<  endl;
    Info <<  "minLen:        "  <<  minLen      <<  endl;
    Info <<  "majNum:        "  <<  majNum      <<  endl;
    Info <<  "minNum:        "  <<  minNum      <<  endl;
    Info <<  "seed:          "  <<  seed        <<  endl;
    Info <<  "roughness:     "  <<  rgh         <<  endl;
    Info <<  "dHurst:        "  <<  dHurst      <<  endl;
    Info <<  "cutLen:        "  <<  cutLen      <<  endl;
    Info <<  "maxDisp:       "  <<  maxDisp     <<  endl;

    Info <<  "Displacement calc starts...." << nl;
    /*
     *   --- ---
     *  | 1 | 2 |
     *   --- ---
     *  | 3 | 4 |
     *   --- ---
     */
    // calculating 1 and 4
    for(int m=0; m<majNum/2+1; ++m)
    {
      for(int n=0; n<minNum/2+1; ++n)
      {
        scalar p = TwoPi * rnd.sample01<scalar>();
        scalar rad;
        double majk, mink, ksq;
        majk = m*cutLen/majLen;
        mink = n*cutLen/minLen;
        ksq   = majk*majk + mink*mink;
        rad = power(ksq) * rnd.GaussNormal<scalar>();

        f[ index(m,n) ] = 
                rad * std::complex<double>(Foam::cos(p),  Foam::sin(p));
        f[ index(((majNum-m)%majNum),(minNum-n)%minNum) ] = 
                rad * std::complex<double>(Foam::cos(p), -Foam::sin(p));
      }
    }
    f[ index(majNum/2,0)   ].imag(0.0);
    f[ index(0,  minNum/2) ].imag(0.0);
    f[ index(majNum/2,minNum/2) ].imag(0.0);

    // calculating 2 and 3
    for(int m=1; m<majNum/2; ++m)
    {
      for(int n=1; n<minNum/2; ++n)
      {
        scalar p = TwoPi * rnd.sample01<scalar>();
        scalar rad;
        double majk, mink, ksq;
        majk = TwoPi*m/majLen;
        mink = TwoPi*n/minLen;
        ksq  = majk*majk + mink*mink;
        rad  = power(ksq) * rnd.GaussNormal<scalar>();

        f[ index(  m, minNum-n) ] = 
                rad * std::complex<double>(Foam::cos(p),  Foam::sin(p));
        f[ index(majNum-m,   n) ] = 
                rad * std::complex<double>(Foam::cos(p), -Foam::sin(p));
      }
    }
    
    fftw_plan plan;
    plan = fftw_plan_dft_2d(majNum, minNum,
                             reinterpret_cast<fftw_complex*>(&f[0]),
                             reinterpret_cast<fftw_complex*>(&F[0]),
                             FFTW_BACKWARD,
                             FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);


    scalarField sF(MN);
    forAll(sF, ii)
    {
      sF[ii] = F[ii].real();
    }

    scalarField sF2 = sqr(sF);
    scalar avSF     = average(sF);
    scalar avSF2    = average(sF2);

    scalar factor = rgh / Foam::sqrt( mag(avSF2 - sqr(avSF)) );

    sF *= factor;

    disp = sF;
}
