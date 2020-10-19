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

Now the macro analysis.C should work.
