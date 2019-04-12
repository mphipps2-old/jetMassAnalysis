ROOUNFOLDDIR = /home/mike/Desktop/JetAnalysis/RooUnfold-1.1.1

CXX      = g++
CXXFLAGS = -Wall `root-config --cflags` -g -o2
LDFLAGS  = -Wall `root-config --glibs` -Wl,-rpath=.$(ROOUNFOLDDIR)

INCLUDEDIRS = -I./include -I$(ROOUNFOLDDIR)/src

LIBDIRS     = -L$(ROOUNFOLDDIR)
LIBS        = -lRooUnfold

#OBJ         = analysis2
OBJ         = analysis
#SOURCES     = getMcWeights.C
#SOURCES     = getMCSpectra.C
#SOURCES     = getDataSpectra.C
#SOURCES     = getCrossSec.C
#SOURCES     = getCrossSecPerEta.C
#SOURCES      = getJES.C
#SOURCES      = getJES2.C
#SOURCES      = getJES2_mpt.C
#SOURCES      = getUnfoldingClosure.C
#SOURCES      = getUnfolding.C
#SOURCES      = getResponseMatrix.C
SOURCES      = getRadim.C
#SOURCES      = getUnfoldingMinimization.C
#SOURCES      = getReweighting.C
#SOURCES      = getJESYongsun.C
#SOURCES      = PlotMassVsPt.C
#SOURCES     = makeUnfMatrix2D.cxx
#INCLUDES    = DiJetAnalysis.h

#DSOURCES    = AtlasStyle.C 
#DINCLUDES   = AtlasStyle.h

CSOURCES    = $(addprefix $(DSOURCES),$(SOURCES))
CINCLUDES   = $(addprefix $(DINCLUDES),$(INCLUDES))

all: $(CSOURCES) $(CINCLUDES)
	$(CXX) $(CXXFLAGS) $(CSOURCES) $(INCLUDEDIRS) \
	$(LIBDIRS) $(LIBS) $(LDFLAGS) -o $(OBJ)

clean :
	rm $(OBJ) 
