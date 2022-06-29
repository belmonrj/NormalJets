// This is a simple executable C++ meant to apply FastJet to Pythia8 events

#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
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

  const int n_jr = 6;
  //const int n_alg = 1;
  double jet_radii[n_jr] = {0.2,0.3,0.4,0.5,0.6,0.7};

  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  pythia.readString("Beams:eCM = 8000.");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  pythia.init();

  // ----------------------
  // --- generate the event
  // ----------------------

  ofstream fout;
  fout.open("list_of_particles.txt");

  // --- get ready for the next event
  if ( !pythia.next() )
    {
      cout << "Um, oops?" << endl;
    }

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
      //double pT  = p.pT();


      // --- double check these
      double px = p.px();
      double py = p.py();
      double pz = p.pz();
      //double E = p.e(); // MC only
      double E = sqrt(px*px + py*py + pz*pz); // can do this, assuming m=0, or assume m=m_pion

      //double y = 0.5*log((E+pz)/(E-pz));

      // add the particles to the FastJet PseudoJet object
      particles.push_back( PseudoJet( px, py, pz, E) );

      //fout << px << " " << py << " " << pz << endl;
      fout << "Event particle " << i << " px py pz " << px << " " << py << " " << pz << endl;

    } // end loop over particles

  fout.close();

  // ----------------------------------------------------------
  // --- Done collecting particles from event, now do jet stuff
  // ----------------------------------------------------------

  JetDefinition jet_def(antikt_algorithm, 0.7);
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  // ---

  int number_of_jets = jets.size();

  // ----------------------------------------------------------------------------------------------------


  TH2D* jet_map_etaphi[number_of_jets];

  fout.open("from_fastjet.txt");
  // --- print out some infos
  fout << "Clustering with " << jet_def.description() << endl;
  // --- print the jets
  fout <<   "        pt y phi" << endl;
  for (int i = 0; i < number_of_jets; i++)
    {
      fout << "jet " << i << ": "<< jets[i].pt() << " " << jets[i].rap() << " " << jets[i].phi() << endl;

      vector<PseudoJet> constituents = jets[i].constituents();

      jet_map_etaphi[i] = new TH2D(Form("th2d_etaphimap_jet%0d",i),"",524,-6,6,256,0,2*pi);

      int number_of_constituents = constituents.size();
      for (int j = 0; j < number_of_constituents; j++)
        {
          fout << "    constituent " << j << " px py pz " << constituents[j].px() << " " << constituents[j].py() << " " <<  constituents[j].pz() << endl;
          jet_map_etaphi[i]->Fill(constituents[j].eta(),constituents[j].phi());
        } // loop over jet constituents

    } // loop over jets

  fout.close();




  // ----------------------------------------------------------------------------------------------------
  // --- different jet radii...

  // --- choose a jet definition
  JetDefinition jet_def_antikt_R0(antikt_algorithm, jet_radii[0]);
  JetDefinition jet_def_antikt_R1(antikt_algorithm, jet_radii[1]);
  JetDefinition jet_def_antikt_R2(antikt_algorithm, jet_radii[2]);
  JetDefinition jet_def_antikt_R3(antikt_algorithm, jet_radii[3]);
  JetDefinition jet_def_antikt_R4(antikt_algorithm, jet_radii[4]);
  JetDefinition jet_def_antikt_R5(antikt_algorithm, jet_radii[5]);

  // --- run the clustering
  ClusterSequence cs_antikt_R0(particles, jet_def_antikt_R0);
  ClusterSequence cs_antikt_R1(particles, jet_def_antikt_R1);
  ClusterSequence cs_antikt_R2(particles, jet_def_antikt_R2);
  ClusterSequence cs_antikt_R3(particles, jet_def_antikt_R3);
  ClusterSequence cs_antikt_R4(particles, jet_def_antikt_R4);
  ClusterSequence cs_antikt_R5(particles, jet_def_antikt_R5);

  // --- get the jets from the clusters
  vector<PseudoJet> jets_antikt_R0 = sorted_by_pt(cs_antikt_R0.inclusive_jets());
  vector<PseudoJet> jets_antikt_R1 = sorted_by_pt(cs_antikt_R1.inclusive_jets());
  vector<PseudoJet> jets_antikt_R2 = sorted_by_pt(cs_antikt_R2.inclusive_jets());
  vector<PseudoJet> jets_antikt_R3 = sorted_by_pt(cs_antikt_R3.inclusive_jets());
  vector<PseudoJet> jets_antikt_R4 = sorted_by_pt(cs_antikt_R4.inclusive_jets());
  vector<PseudoJet> jets_antikt_R5 = sorted_by_pt(cs_antikt_R5.inclusive_jets());


  // --- all done



  pythia.stat(); // tell about some statistics for this run



  //Tfile for I/O stuff
  TFile* JetHistFile = new TFile("testout.root","recreate");

  for (int i = 0; i < number_of_jets; i++)
    {
      jet_map_etaphi[i]->Write();
    }

  JetHistFile->Close();

  return 0;

} // end of int main()
