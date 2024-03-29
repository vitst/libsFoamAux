/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "field2tiff.H"
#include "interpolation.H"
#include "interpolationCellPoint.H"
#include "triSurfaceSearch.H"
#include "meshSearch.H"

#include <hdf5/serial/hdf5.h>
#include <hdf5/serial/hdf5_hl.h>
//#include <hdf5.h>
//#include <hdf5_hl.h>

#include <assert.h>
// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::field2tiff::foundObject
(
    const word& name,
    const bool verbose
) const
{
    if (fvMeshFunctionObject::foundObject<Type>(name))
    {
        return true;
    }
    else
    {
        if (debug || verbose)
        {
            Warning
                << "    functionObjects::" << type() << " " << this->name()
                << " cannot find required object " << name << " of type "
                << Type::typeName << endl;
        }

        return false;
    }
}

template<class Type>
void Foam::functionObjects::field2tiff::convertFields
(
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const fieldType& field = lookupObject<fieldType>(fieldName);
        autoPtr< interpolation<Type> > interpolator = 
                interpolation<Type>::New("cellPoint", field);
        meshSearch searchEngine(mesh_);

        //word fileName = fieldName + ".tif";
        fileName current_postPr_dir = 
                 "postProcessing/field2tiff" / mesh_.time().timeName();
        if ( !isDir(current_postPr_dir) ) mkDir(current_postPr_dir);

        //!!!fileName fN = current_postPr_dir / fieldName + ".tif";
        fileName fN = current_postPr_dir / fieldName + ".hdf5";

/*
        scalar minX = gMin( mesh_.points().component(vector::X) );
        scalar maxX = gMax( mesh_.points().component(vector::X) );
        scalar minY = gMin( mesh_.points().component(vector::Y) );
        scalar maxY = gMax( mesh_.points().component(vector::Y) );
        scalar minZ = gMin( mesh_.points().component(vector::Z) );
        scalar maxZ = gMax( mesh_.points().component(vector::Z) );

        int nX = static_cast<int>( (maxX - minX) / imageUnit_ ); 
        int nY = static_cast<int>( (maxY - minY) / imageUnit_ ); 
        int nZ = static_cast<int>( (maxZ - minZ) / imageUnit_ ); 
*/
        int nX = static_cast<int>( (maxPoint_.x() - minPoint_.x()) / imageUnit_ ); 
        int nY = static_cast<int>( (maxPoint_.y() - minPoint_.y()) / imageUnit_ ); 
        int nZ = static_cast<int>( (maxPoint_.z() - minPoint_.z()) / imageUnit_ ); 

        Info<< "file: " << fN << " to write at time: "
                        << mesh_.time().timeName() << nl;
        Info<< "Image dimensions: " << nX << "  "<< nY << "  " << nZ <<nl<<endl;

        /* !!!
        TIFF *outTiff = TIFFOpen(fN.c_str(), "w8");
        uint32 x_dim = static_cast<uint32>(nX);
        uint32 y_dim = static_cast<uint32>(nY);

        int sampleperpixel = 1; // for B/W is 1
        uint8 bitspersample  = 32;
    
        uint32 nstacks = static_cast<uint32>(nZ);
        */
        
        Info << "[----------]"<<endl;
        Info << "[" << flush;
        //!!!int nPrint = static_cast<int>(nstacks / 10.0);
        int nPrint = static_cast<int>(nZ / 10.0);
        int countPrint = 0;

        //float data[ x_dim * y_dim * z_dim ];
        double *data;
        data = (double *) malloc(nX*nY*nZ*sizeof(*data));

        // loop over tiff directions
        //for (uint32 iZ = 0; iZ < nstacks; iZ++)
        for (uint32 iZ = 0; iZ < nZ; iZ++)
        {
            /* !!!
            TIFFSetField(outTiff, TIFFTAG_IMAGELENGTH, y_dim);
            TIFFSetField(outTiff, TIFFTAG_IMAGEWIDTH, x_dim);
            TIFFSetField(outTiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    
            TIFFSetField (outTiff, TIFFTAG_ROWSPERSTRIP, 1);
    
            TIFFSetField(outTiff, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);
            TIFFSetField(outTiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
            TIFFSetField(outTiff, TIFFTAG_BITSPERSAMPLE, bitspersample);
    
            // for 32!!!!!
            if(bitspersample>16)
                TIFFSetField(outTiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    
            TIFFSetField(outTiff, TIFFTAG_ROWSPERSTRIP, 
                    TIFFDefaultStripSize(outTiff, x_dim*sampleperpixel));
    
            TIFFSetField(outTiff, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
            TIFFSetField(outTiff, TIFFTAG_PAGENUMBER, iZ, nstacks);
    
            // not sure if this is needed...
            TIFFSetField(outTiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
            TIFFSetField(outTiff, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
            //TIFFSetField(outTiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            TIFFSetField(outTiff, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
    
            float data[ x_dim * y_dim ]; moved up
            */

            double z_pos = minPoint_.z() + (iZ + 0.5) * imageUnit_;

            //for(uint32 iY=0; iY < y_dim; iY++)
            for(uint32 iY=0; iY < nY; iY++)
            {
                double y_pos = minPoint_.y() + (iY + 0.5) * imageUnit_;
                //for(uint32 iX=0; iX < x_dim; iX++)
                for(uint32 iX=0; iX < nX; iX++)
                {
                    //uint32 ixy = iY + iX * y_dim;
                    //!!!uint32 ixy = iX + iY * x_dim;
                    double x_pos = minPoint_.x() + (iX + 0.5) * imageUnit_;

                    point samp_point(x_pos, y_pos, z_pos);
                    label cellI = searchEngine.findCell( samp_point );
                    //if (cellI==-1)
                    //{
                        //cellI = searchEngine.findNearestCell( samp_point );
                    //}
                    label faceI = -1;

                    Type interp_field = pTraits<Type>::zero;
                    if (cellI != -1)
                    {
                        interp_field = interpolator().interpolate(samp_point, cellI, faceI);
                    }
                    
                    //Info << ixy << ":  " << iX << "  " << iY << "  " << iZ << nl;
                    //Info << x_pos << "  " << y_pos << "  " << z_pos << " :  " << mag( interp_field )
                    //     << "  is inside: " << searchEngine.findCell( samp_point ) << nl;
                    //!!!data[ixy] = static_cast<float>( mag( interp_field ) );
                    data[ iZ*nY*nX + iY*nX + iX ] = mag( interp_field );
                }
            }
            /*!!! 
            char* dat = (char*)data;
    
            for (uint32 row = 0; row < y_dim; row++) 
            {
                int ret = TIFFWriteScanline(outTiff, dat+((y_dim - row - 1)*x_dim)*sizeof(float), row, 0);
            }

            TIFFWriteDirectory(outTiff);
            */
            countPrint++;
            if(countPrint==nPrint)
            {
                Info << "*" << flush;
                countPrint=0;
            }
            
        }
        Info << "]" << endl;
        //!!!TIFFClose(outTiff);
        Info << "Writing HDF5..."<<endl;

        /* now the various steps involved in preparing an hdf5 file */
        hid_t file_id;
        /* at this time our basic file just has the values strung out,
           so the hdf5 rank is 1 */
        hsize_t dims[3] = {nX, nY, nZ};
        herr_t status;
      
        /* create a HDF5 file */
        file_id = H5Fcreate(fN.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        /* create and write a double type dataset named "/tomodata" */
        status = H5LTmake_dataset(file_id, "/tomodata", 3, dims,
                                  H5T_NATIVE_DOUBLE, data);
        assert(status != -1);
        /* add some hdf5 attributes with the metadata we need */
        H5LTset_attribute_int(file_id, "/tomodata", "nx", &nX, 1);
        H5LTset_attribute_int(file_id, "/tomodata", "ny", &nY, 1);
        H5LTset_attribute_int(file_id, "/tomodata", "nz", &nZ, 1);
      
        status = H5Fclose (file_id); assert(status != -1);
        assert(status != -1);

    }
}



// ************************************************************************* //
