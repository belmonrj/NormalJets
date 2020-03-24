// This is a simple executable C++ meant to apply FastJet to Pythia8 events

#include "TROOT.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"

#include "Pythia8/Pythia.h"

#include "fastjet/ClusterSequence.hh"

//Added this to write .off file
#include <fstream>
#include <iostream>

using namespace fastjet;
using namespace Pythia8;
using namespace std;



//const double pi = 3.14159265358979323;

int main()
{

  // initialize a new ROOT histogram to fill with the loop
  TH1D* pTr = new TH1D("pTr", "Transverse Momentum", 100, 0, 100);
  TH1D* pTj = new TH1D("pTj", "Jet Transverse Momentum", 100, 0, 100);


  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  pythia.init();

  // -------------------
  // --- Main event loop
  // -------------------

  for ( int iEvent = 0; iEvent < 2; ++iEvent )
    {

      // --- get ready for the next event
      if ( !pythia.next() ) continue;

      // --- generate the event
      Event& event = pythia.event;

      // --- make the object for FastJet
      vector<PseudoJet> particles;

      // --- loop over the particles in the event
      for (int i = 0; i < event.size(); ++i)
        {
          // Particle object contains various particle-based information
          Particle& p = event[i];

          if ( !p.isFinal() ) continue; // only want final state particles
          // --- these may be useful at a later time
          // bool charge = p.isCharged();
          // double phi = p.phi();
          // double eta = p.eta();
          double pT  = p.pT();
	  pTr->Fill( pT );

          // --- double check these
          double px = p.px();
          double py = p.py();
          double pz = p.pz();
	  double E = p.e();


          // add the particles to the FastJet PseudoJet object
          particles.push_back( PseudoJet( px, py, pz, E) );

        } // end loop over particles

      // ----------------------------------------------------------
      // --- Done collecting particles from event, now do jet stuff
      // ----------------------------------------------------------

      // --- choose a jet definition
      double R = 0.7;
      JetDefinition jet_def(antikt_algorithm, R);

      // --- run the clustering, extract the jets
      ClusterSequence cs(particles, jet_def);
      vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

      // --- print out some infos
      cout << "Clustering with " << jet_def.description() << endl;

      // --- print the jets
      cout <<   "        pt y phi" << endl;
      for (unsigned i = 0; i < jets.size(); i++)
        {
          cout << "jet " << i << ": "<< jets[i].pt() << " " << jets[i].rap() << " " << jets[i].phi() << endl;
          vector<PseudoJet> constituents = jets[i].constituents();
	  pTj->Fill(jets[i].pt());
          for (unsigned j = 0; j < constituents.size(); j++)
            {
              cout << "    constituent " << j << "'s pt: " << constituents[j].pt() << endl;
            } // loop over jet constituents
        } // loop over jets

      // --- all done

    } // end of loop over events

  pythia.stat(); // tell about some statistics for this run



  //Tfile for I/O stuff
  TFile* JetHistFile = new TFile("testout.root","recreate");
  pTr->Write();
  pTj->Write();
  JetHistFile->Close();

  return 0;

} // end of int main()
