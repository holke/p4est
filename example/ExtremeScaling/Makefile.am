
# This file is part of p4est.
# Makefile.am in example/simple
# included non-recursively from toplevel directory

#if P4EST_ENABLE_BUILD_2D
#bin_PROGRAMS += example/ExtremeScaling/p4est_bunny
#example_ExtremeScaling_p4est_bunny_SOURCES = example/ExtremeScaling/bunny.c

#LINT_CSOURCES += $(example_simple_p4est_simple_SOURCES)
#endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += \
        example/ExtremeScaling/p4est_exchange

example_ExtremeScaling_p4est_exchange_SOURCES = example/ExtremeScaling/exchange2.c

LINT_CSOURCES += \
        $(example_ExtremeScale_p4est_exchange_SOURCES)
endif

if P4EST_ENABLE_BUILD_3D
bin_PROGRAMS += \
        example/ExtremeScaling/p8est_bunny \
        example/ExtremeScaling/p8est_exchange

example_ExtremeScaling_p8est_bunny_SOURCES = example/ExtremeScaling/bunny.c
example_ExtremeScaling_p8est_exchange_SOURCES = example/ExtremeScaling/exchange3.c

LINT_CSOURCES += \
        $(example_ExtremeScale_p8est_bunny_SOURCES) \
        $(example_ExtremeScale_p8est_exchange_SOURCES)
endif
