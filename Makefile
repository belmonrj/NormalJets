FFLAGS = `fastjet-config --cxxflags --libs --plugins`

PFLAGS = `pythia8-config --cxxflags --libs`

RFLAGS = `root-config --cflags --libs`

PROGRAMS = fj_test fj_pythia_test

all: $(PROGRAMS)

clean:
	rm -f $(PROGRAMS)

fj_test: fj_test.C
	g++ -o fj_test fj_test.C -Wall $(FFLAGS)

fj_pythia_test: fj_pythia_test.C
	g++ -o fj_pythia_test fj_pythia_test.C -Wall $(FFLAGS) $(PFLAGS) $(RFLAGS)

