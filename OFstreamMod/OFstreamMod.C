/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "OFstreamMod.H"
#include "OSspecific.H"
#include "gzstream.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(OFstreamMod, 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

OFstreamModAllocator::OFstreamModAllocator
(
    const fileName& pathname,
    ios_base::openmode mode,
    IOstream::compressionType compression
)
:
    ofPtr_(NULL)
{
    if (!pathname.size())
    {
        if (OFstreamMod::debug)
        {
            Info
                << "OFstreamModAllocator::OFstreamModAllocator"
                   "(const fileName& pathname) : "
                   "can't open null file "
                << endl;
        }
    }

    if (compression == IOstream::COMPRESSED)
    {
        if (isFile(pathname,false))
        {
            rm(pathname);
        }

        ogzstream* gzofPtr = new ogzstream((pathname + ".gz").c_str(), mode);
        ofPtr_ = gzofPtr;
    }
    else
    {
        if (isFile(pathname + ".gz", false))
        {
            rm(pathname + ".gz");
        }

        ofPtr_ = new ofstream(pathname.c_str(), mode);
    }
}


OFstreamModAllocator::~OFstreamModAllocator()
{
    delete ofPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

OFstreamMod::OFstreamMod
(
    const fileName& pathname,
    ios_base::openmode mode,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OFstreamModAllocator(pathname, mode, compression),
    OSstream(*ofPtr_, "OFstreamMod.sinkFile_", format, version, compression),
    pathname_(pathname)
{
    setClosed();

    setState(ofPtr_->rdstate());
                
    if (!good())
    {
        if (debug)
        {
            Info<< "IFstream::IFstream(const fileName& pathname,"
                   "streamFormat format=ASCII,"
                   "versionNumber version=currentVersion) : "
                   "couldn't open File for input\n"
                   "in stream " << info() << Foam::endl;
        }

        setBad();
    }
    else
    {
        setOpened();
    }
    
    lineNumber_ = 1;
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

OFstreamMod::~OFstreamMod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void OFstreamMod::print(Ostream& os) const
{
    // Print File data
    os  << "    OFstreamMod: ";
    OSstream::print(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
