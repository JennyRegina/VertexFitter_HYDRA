
# VertexFitter_HYDRA

Workspace for the development of the vertex fitter

Description how to run the external library:

git clone https://github.com/JennyRegina/VertexFitter_HYDRA.git YOUR_EXTERNAL_LIB_FOLDER
(where YOUR_EXTERNAL_LIB_FOLDER is the name of the folder where you want to install your library)

cd YOUR_EXTERNAL_LIB_FOLDER

git checkout decaybuilder_jana

Open singularity container

(this gives you the correct version of the HYDRA code)

export MYHADDIR=/lustre/hades/user/YOUR_USER_FOLDER/YOUR_EXTERNAL_LIB_FOLDER 

. /cvmfs/hades.gsi.de/install/5.34.34/hydra2-5.3/defall.sh

(here you can use your preferred HYDRA version)

mkdir bin build include install lib macros

(this creates the structure you want for your external library)

Copy all the files you downloaded from git into the folder lib

cd $MYHADDIR/lib

make 

make install

copy rootlogon.C into $MYHADDIR/macros

(you might need to run: export ROOTLOGON=/PATH_TO_YOUR_ROOTLOGON/rootlogon.C depending on how your .rootrc looks)

root -b

add the library libKinFit.so to the rootlogon.C macro if it is not present:

common_libs +="KinFit";

You might need to remove other external libraries from your rootlogon.C.
Alternatively you can use the rootlogon provided here.

Now the macro analysisVertexFinder_million_Geant.C should work.

The macro needs to be run pre-compiled. To run the macro on virgo, type

root -l -b -q analysisVertexFinder_million_GeantInfo.C+

This macro analysis.C illustrates how to use the vertex finding and fitting classes together. The output is a root file called "testvertexfit.root" that contains some example histograms. This macro currently needs to be updated to work for the new structure of the vertexfitter. This macro will be updated to reflect the newest changes.

How the classes work:

analysisVertexFinder_million_GeantInfo.C is an example macro of how the vertex finding and fitting and 3C fit can be used in an exclusive event where the decay chain is known from the Geant information so it gives the best case scenario.

A more realistic analysis can be found in the macro analysisVertexFinder_million_Realistic.C. Here cuts on the vertex fit probabilities are used to select the veritices from which the neutral mother candidate is calculated. No 3C fit is done yet here. 

The hvertexfinder.cc class finds the primary vertex as well as the decay vertex. 

The class hneutralcandfinder.cc creates the neutral canidate. Currently the primary vertex has to be used in this creation. 

The vertexfitter.cc class can perform two types of fitting procedures, first a fit where the tracks are constrained to originate from the same point, i.e. a vertex constraint. This can be done for both the primary and secondary vertices. It can also take as input the neutral candidate created in the hvertexfinder.cc class together with its daughters and perform a fit with a 3C constraint utilizing the information from this object. Here the mass of the neutral candidate is fixed and its momentum is an unmeasured variable. Its direction is taken from the neutral candidate constructed by the vertex finder. The daughters should be fitted to the secondary vertex before doing the 3c fit.
