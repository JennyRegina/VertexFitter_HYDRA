##############################################################################
#
#  Makefile for Hydra library libVertexFit.so
#
#  This makefile contains all definitions local to this module. All
#  general definitions are included from makefiles named "hydra.*.mk".
#
##############################################################################


LIB_NAME := KinFit

USES_ORACLE= : no

ifeq ($(HADDIR), /lustre/nyx/hades/user/kempter/svn/hydra_BT)
    $(info Running with Data)
    $(info $(HADDIR))
SOURCE_FILES := hdstfitter.cc \
        hdecaybuilder.cc \		
		hcovariancekinfit.cc \
        hrefitcand.cc \
		hkinfitter.cc \
		hvertexfinder.cc \
        hneutralcandfinder.cc 
else
    $(info Running with simulations)
    $(info $(HADDIR))
SOURCE_FILES := hdstfitter.cc \
        hdecaybuilder.cc \		
		hcovariancekinfit.cc \
        hrefitcand.cc \
		hkinfitter.cc \
		hvertexfinder.cc \
        hneutralcandfinder.cc 
endif

include $(HADDIR)/hades.def.mk

# set this, while debugging
SO_CXX_FLAGS += -O0

include $(HADDIR)/hades.module.mk
