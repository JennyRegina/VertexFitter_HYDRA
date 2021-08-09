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

SOURCE_FILES := hdecaybuilder.cc \
        hrefitcand.cc \
		hkinfitter.cc \
		hvertexfinder.cc \
        hneutralcandfinder.cc \
		hdecaybuilder.cc \
		hdstfitter.cc \
		hcovariancekinfit.cc \


include $(HADDIR)/hades.def.mk

# set this, while debugging
SO_CXX_FLAGS += -O0

include $(HADDIR)/hades.module.mk
