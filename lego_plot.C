void lego_plot()
{

  TCanvas c1("c1","");

  // there's a smarter way to do this but let's be stupid for now
  const int how_many = 10;
  TH2D* jetmap[how_many];
  TFile* fin = TFile::Open("testout.root");
  for ( int i = 0; i < how_many; ++i )
    {
      jetmap[i] = (TH2D*)fin->Get(Form("th2d_etaphimap_jet%d",i));
    }

  // --- if the first one is null we can't do anything
  if ( jetmap[0] == NULL )
    {
      cout << "You're gonna have a bad time" << endl;
      return;
    }

  // --- the plots have too many bins so let's try rebinning
  // --- they also have too big a range but we'll worry about that later
  for ( int i = 0; i < how_many; ++i )
    {
      jetmap[i]->RebinX(10);
      jetmap[i]->RebinY(10);
    }


  // have the maps, now let's draw them

  jetmap[0]->SetLineColor(2);
  jetmap[0]->Draw("lego");
  c1.Print("map_test_lego_one.png");
  for ( int i = 1; i < how_many; ++i )
    {
      jetmap[i]->SetLineColor(i+2);
      jetmap[i]->Draw("same lego");
    }

  c1.Print("map_test_lego_some.png");

}
