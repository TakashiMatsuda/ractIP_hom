bin_PROGRAMS = ractip_hom
SUBDIRS = centroidalifold/contralign centroidalifold/probconsRNA
ractip_hom_SOURCES = \
	ractip_hom.cpp ip.cpp ip.h fa.cpp fa.h \
	boltzmann_param.c boltzmann_param.h \
	cmdline.c cmdline.h pf_duplex.c pf_duplex.h ushuffle.c ushuffle.h \
	contrafold/InferenceEngine.hpp contrafold/InferenceEngine.ipp \
	contrafold/DuplexEngine.hpp contrafold/DuplexEngine.ipp \
	contrafold/ParameterManager.hpp contrafold/ParameterManager.ipp \
	contrafold/SStruct.hpp contrafold/SStruct.cpp \
	contrafold/Utilities.hpp contrafold/Utilities.ipp contrafold/Utilities.cpp \
	contrafold/Defaults.ipp contrafold/LogSpace.hpp contrafold/Config.hpp contrafold/rand.h \
	centroidalifold/engine/contrafold.h centroidalifold/engine/contrafold.cpp \
	centroidalifold/engine/contrafoldm.h centroidalifold/engine/contrafoldm.cpp \
	centroidalifold/engine/mccaskill.h centroidalifold/engine/mccaskill.cpp \
	centroidalifold/engine/alifold.h centroidalifold/engine/alifold.cpp \
	centroidalifold/engine/pfold.h centroidalifold/engine/pfold.cpp \
	centroidalifold/engine/averaged.h centroidalifold/engine/averaged.cpp \
	centroidalifold/engine/mixture.h centroidalifold/engine/mixture.cpp \
	centroidalifold/engine/aux.h centroidalifold/engine/aux.cpp \
	centroidalifold/folding_engine.h centroidalifold/folding_engine.cpp \
	centroidalifold/folding_engine.cpp centroidalifold/folding_engine.h \
	centroidalifold/bp.cpp centroidalifold/bp.h \
	centroidalifold/aln.h centroidalifold/aln.cpp \
	centroidalifold/cyktable.h \
	centroidalifold/centroid.h centroidalifold/centroid.cpp \
	centroidalifold/mea.h centroidalifold/mea.cpp \
	centroidalifold/diana.h centroidalifold/diana.cpp \
	centroidalifold/ps_plot.h centroidalifold/ps_plot.cpp \
	centroidalifold/rule.h \
	centroidalifold/contrafold/contrafold.h \
	centroidalifold/contrafold/wrapper.cpp \
	centroidalifold/contrafold/Defaults.ipp \
	centroidalifold/engine/contrafoldhom.h centroidalifold/engine/contrafoldhom.cpp \
	centroidalifold/engine/mccaskillhom.h centroidalifold/engine/mccaskillhom.cpp \
	centroidalifold/util.h 
ractip_hom_LDADD = centroidalifold/contralign/libcontralign.a centroidalifold/probconsRNA/libprobcons.a
