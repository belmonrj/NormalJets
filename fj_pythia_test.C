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

  const int n_jr = 6;
  const int n_alg = 1;
  double jet_radii[n_jr] = {0.2,0.3,0.4,0.5,0.6,0.7};

  TH1D* pTr = new TH1D("pTr", "Transverse Momentum", 100, 0, 100);
  TH1D* pTj = new TH1D("pTj", "Jet Transverse Momentum", 100, 0, 100);
  // = new TH1D("pTr", "Transverse Momentum", 100, 0, 100);
  TH1D* pTj_set[n_jr][n_alg];
  for ( int i = 0; i < n_jr; ++i )
    {
      for ( int j = 0; j < n_alg; ++j )
        {
          pTj_set[i][j] = new TH1D(Form("pTj_set_R%d_A%d",i,j), "", 100, 0, 100);
        }
    }

  TH1D* numConst[n_jr][n_alg]; // = new TH1D("numConst","Constiuets per Jet Radii", 100, 0, 100);
    for (int k = 0; k < n_jr; ++k)
    {
      for ( int j = 0; j < n_alg; ++j)
	{
	  numConst[k][j] = new TH1D(Form("numConst_R%d",k),"", 100, 0, 100);
	}
    }

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

  int n_events = 2;

  for ( int iEvent = 0; iEvent < n_events; ++iEvent )
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
          pTr->Fill( pT ); // pT distribution of single particles

          // --- double check these
          double px = p.px();
          double py = p.py();
          double pz = p.pz();
          //double E = p.e(); // MC only
          double E = sqrt(px*px + py*py + pz*pz); // can do this, assuming m=0, or assume m=m_pion

          // add the particles to the FastJet PseudoJet object
          particles.push_back( PseudoJet( px, py, pz, E) );

        } // end loop over particles

      // ----------------------------------------------------------
      // --- Done collecting particles from event, now do jet stuff
      // ----------------------------------------------------------

      JetDefinition jet_def(antikt_algorithm, 0.7);
      ClusterSequence cs(particles, jet_def);
      vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

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

      // --- fill the jet pt histograms
      for ( unsigned int i = 0; i < jets_antikt_R0.size(); ++i ) pTj_set[0][0]->Fill(jets_antikt_R0[i].pt());
      for ( unsigned int i = 0; i < jets_antikt_R1.size(); ++i ) pTj_set[1][0]->Fill(jets_antikt_R1[i].pt());
      for ( unsigned int i = 0; i < jets_antikt_R2.size(); ++i ) pTj_set[2][0]->Fill(jets_antikt_R2[i].pt());
      for ( unsigned int i = 0; i < jets_antikt_R3.size(); ++i ) pTj_set[3][0]->Fill(jets_antikt_R3[i].pt());
      for ( unsigned int i = 0; i < jets_antikt_R4.size(); ++i ) pTj_set[4][0]->Fill(jets_antikt_R4[i].pt());
      for ( unsigned int i = 0; i < jets_antikt_R5.size(); ++i ) pTj_set[5][0]->Fill(jets_antikt_R5[i].pt());

      // --- fill numConst histograms
      for ( unsigned int i = 0; i < jets_antikt_R0.size(); i++ ) numConst[0][0]->Fill(jets_antikt_R0[i].constituents().size());
      for ( unsigned int i = 0; i < jets_antikt_R1.size(); i++ ) numConst[1][0]->Fill(jets_antikt_R1[i].constituents().size());
      for ( unsigned int i = 0; i < jets_antikt_R2.size(); i++ ) numConst[2][0]->Fill(jets_antikt_R2[i].constituents().size());
      for ( unsigned int i = 0; i < jets_antikt_R3.size(); i++ ) numConst[3][0]->Fill(jets_antikt_R3[i].constituents().size());
      for ( unsigned int i = 0; i < jets_antikt_R4.size(); i++ ) numConst[4][0]->Fill(jets_antikt_R4[i].constituents().size());
      for ( unsigned int i = 0; i < jets_antikt_R5.size(); i++ ) numConst[5][0]->Fill(jets_antikt_R5[i].constituents().size());

      // ----------------------------------------------------------------------------------------------------

      if ( n_events < 6 )
        {
          // --- print out some infos
          cout << "Clustering with " << jet_def.description() << endl;
          // --- print the jets
          cout <<   "        pt y phi" << endl;
          for (unsigned i = 0; i < jets.size(); i++)
            {
              cout << "jet " << i << ": "<< jets[i].pt() << " " << jets[i].rap() << " " << jets[i].phi() << endl;
              vector<PseudoJet> constituents = jets[i].constituents();
              pTj->Fill(jets[i].pt());
              unsigned int number_of_constituents = constituents.size();
              for (unsigned j = 0; j < number_of_constituents; j++)
                {
                  cout << "    constituent " << j << "'s pt: " << constituents[j].pt() << endl;
                } // loop over jet constituents
            } // loop over jets
        }

      // ----------------------------------------------------------------------------------------------------

      // --- all done

    } // end of loop over events

  pythia.stat(); // tell about some statistics for this run



  //Tfile for I/O stuff
  TFile* JetHistFile = new TFile("testout.root","recreate");
  pTr->Write();
  pTj->Write();
  for ( int i = 0; i < n_jr; ++i )
    {
      for ( int j = 0; j < n_alg; ++j )
        {
          pTj_set[i][j]->Write();
        }
    }

for (int i=0; i < n_jr; ++i )
  {
    for (int j=0; j < n_alg; ++j)
      {
	numConst[i][j]->Write();
      }
  }
  JetHistFile->Close();

  return 0;

} // end of int main()
