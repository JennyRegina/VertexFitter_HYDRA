#!/bin/bash

#. /cvmfs/hades.gsi.de/install/5.34.34/hydra2-4.9m/defall.sh
. /cvmfs/hades.gsi.de/install/root-5.34.34/bin/thisroot.sh 
#. /cvmfs/hades.gsi.de/install/root-6.12.06/bin/thisroot.sh

export HADDIR=/lustre/nyx/hades/user/rlalik/fwdet/install/hydra2-fwdet
export CERN_ROOT=/cvmfs/hades.gsi.de/install/cernlib_gfortran/2005/
export HGEANT_DIR=/lustre/nyx/hades/user/rlalik/fwdet/install/hgeant2-fwdet

export ORA_USER=hades_ana/hades@db-hades
export ORACLE_HOME=/cvmfs/it.gsi.de/oracle/product/12.1.2/client_x86_64_1

export LC_ALL=C

#export MYHADDIR=/lustre/nyx/hades/user/${USER}/hades/pp45/install
export MYHADDIR=/lustre/hades/user/jrieger/Hydra/myLibs
export ROOTLOGON=$MYHADDIR/macros/rootlogon.C

INSTALL_DIR=/lustre/nyx/hades/user/rlalik/hades/install/5.34.34
export PATH=${INSTALL_DIR}/bin:${MYHADDIR}/bin:${PATH}
#export LD_LIBRARY_PATH=${INSTALL_DIR}/lib:${MYHADDIR}/lib:${PLUTODIR}:${HADDIR}/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=${INSTALL_DIR}/lib:${MYHADDIR}/lib:${HADDIR}/lib:${LD_LIBRARY_PATH}
#export LD_LIBRARY_PATH=/cvmfs/it.gsi.de/oracle/product/12.1.2/client_x86_64_1/lib/:${LD_LIBRARY_PATH}
