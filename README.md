
# VertexFitter_HYDRA

Workspace for the development of the vertex fitter

Description how to run the external library:

Enviornment Variables:
export MYHADDIR=/lustre/hades/user/YOUR_USER_FOLDER/YOUR_EXTERNAL_LIB_FOLDER
. /cvmfs/hades.gsi.de/install/5.34.34/hydra2-5.3/defall.sh

(it should work to run the external library with HYDRA versions around 5.3)

To run the Vertex fit:

Place it in $MYHADDIR/lib

run make 
and make install

add the library libVertexFit.so to the rootlogon.C macro:

common_libs +="VertexFit";

You might need to remove other external libraries from rootlogon.C.

Now the macro analysisVertexFinder.C should work.

The macro needs to be run pre-compiled. To run the macro on virgo, type

root -b

and in root session

.X analysisVertexFinder.C+

How the classes work:

The hvertexfinder.cc class finds the primary vertex as well as the decay vertex. It also creates the neutral canidate. Currently the primary vertex is not used in this creation. 

The vertexfitter.cc class can perform two types of fitting procedures, first a fit where the tracks are constrained to originate from the same point, i.e. a vertex constraint. This can be done for both the primary and secondary vertices. It can also take as input this neutral candidate created in the hvertexfinder.cc class and perform a fit with a 3C constraint utilizing the information from this object. 
