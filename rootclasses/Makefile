ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
EVENTLIB      = $(EVENTO)
OutPutOpt     = -o
 
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
 
# Solaris
CXX           = CC
CXXFLAGS      = -O -KPIC
LD            = CC
LDFLAGS       = -g
SOFLAGS       = -G
 
 
CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) -lcomplex
GLIBS         = $(ROOTGLIBS) -lcomplex
 
#------------------------------------------------------------------------------

HDRS = TMCdata.h THitsData.h TReconstruction.h TResolver.h TChannel.h

OBJS =  TMCdata.$(ObjSuf) TReconstruction.$(ObjSuf) THitsData.$(ObjSuf) \
        TChannel.$(ObjSuf) TResolver.$(ObjSuf) TTable.$(ObjSuf) \
        RootDictionary.$(ObjSuf)

LIBRARY = rootlib.$(DllSuf)

.SUFFIXES: .$(SrcSuf)

all : $(LIBRARY)

tester : tester.$(ObjSuf) $(OBJS)
	@echo "Linking $@ ..." 
	@CC -o $@ $^ $(CXXFLAGS) $(GLIBS)
	@echo "done"

myroot : main.$(ObjSuf) $(OBJS)
	@echo "Linking $@ ..." 
	@$(LD) $(LDFLAGS) $< $(OBJS) $(GLIBS) $(OutPutOpt) $@
	@echo "done"

$(LIBRARY) : $(OBJS)
	@echo "Building shared library ..."
	@$(LD) $(SOFLAGS) $(LDFLAGS) $(OBJS) $(OutPutOpt) $@
	@echo "done"

testTable : TTable.$(SrcSuf)
	@echo "Linking $@ ..."
	@CC -D VERIFIER $^ $(CXXFLAGS) $(GLIBS) $(OutPutOpt) $@
	@echo "done"

clean:
	@rm -f $(OBJS) RootDictionary.* core

RootDictionary.$(ObjSuf) : $(HDRS)
RootDictionary.$(SrcSuf) : $(HDRS) LinkDef.h
	@echo "Generating $@ from $(HDRS)"
	@$(ROOTSYS)/bin/rootcint -f $@ -c $^

TMCdata.$(ObjSuf) : TMCdata.h

THitsData.$(ObjSuf) : THitsData.h

TReconstruction.$(ObjSuf) : TReconstruction.h TChannel.h

TChannel.$(ObjSuf) : TChannel.h

TPi0s.$(ObjSuf) : TPi0s.h

TResolver.$(ObjSuf) : TResolver.h TChannel.h

.$(SrcSuf).$(ObjSuf):
	@echo "Making $@"
	@$(CXX) $(CXXFLAGS) -c $<

