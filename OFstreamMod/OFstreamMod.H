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

Class
    OFstreamMod

Description
    Output to file stream.

SourceFiles
    OFstreamMod.C

\*---------------------------------------------------------------------------*/

#ifndef OFstreamMod_H
#define OFstreamMod_H

#include "OSstream.H"
#include "fileName.H"
#include "className.H"

#include <fstream>
using std::ofstream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OFstreamMod;

/*---------------------------------------------------------------------------*\
                       Class OFstreamModAllocator Declaration
\*---------------------------------------------------------------------------*/

class OFstreamModAllocator
{
    friend class OFstreamMod;

    ostream* ofPtr_;

    // Constructors

        //- Construct from pathname
        OFstreamModAllocator
        (
            const fileName& pathname,
            ios_base::openmode mode,
            IOstream::compressionType compression
        );

    // Destructor

        ~OFstreamModAllocator();
};


/*---------------------------------------------------------------------------*\
                           Class OFstreamMod Declaration
\*---------------------------------------------------------------------------*/

class OFstreamMod
:
    private OFstreamModAllocator,
    public OSstream
{
    // Private data

        fileName pathname_;


public:

    // Declare name of the class and it's debug switch
    ClassName("OFstreamMod");


    // Constructors

        //- Construct from pathname
        OFstreamMod
        (
            const fileName& pathname,
            ios_base::openmode mode = ios_base::out|ios_base::trunc,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        );


    // Destructor

        ~OFstreamMod();


    // Member functions

        // Access

            //- Return the name of the stream
            const fileName& name() const
            {
                return pathname_;
            }

            //- Return non-const access to the name of the stream
            fileName& name()
            {
                return pathname_;
            }


        // Print

            //- Print description of IOstream to Ostream
            void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global predefined null output stream

//extern OFstreamMod Snull;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
