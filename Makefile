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
SOURCE_FILES := hrefitcand.cc \
		hkinfitter.cc \
		hvertexfinder.cc \
        hneutralcandfinder.cc 
    #export DATA_SYS_KINFIT=1
    #$(info $(DATA_SYS_KINFIT))
else
    $(info Running with simulations)
    $(info $(HADDIR))
SOURCE_FILES := hrefitcand.cc \
		hkinfitter.cc \
		hvertexfinder.cc \
        hneutralcandfinder.cc 
    #export SIM_SYS_KINFIT=$(HADDIR)
endif

#SOURCE_FILES := hrefitcand.cc \
#		hkinfitter.cc \
#		hvertexfinder.cc \
#        hneutralcandfinder.cc   

include $(HADDIR)/hades.def.mk

# set this, while debugging
SO_CXX_FLAGS += -O0

include $(HADDIR)/hades.module.mk
