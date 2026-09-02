// Example: root -l 'inspect.C("clumps.root")'
void inspect(const char *fname = "clumps.root")
{
   TFile *f = TFile::Open(fname);
   TTree *t = (TTree *)f->Get("clumps");

   t->Print();                                  // branch listing
   printf("entries = %lld\n", t->GetEntries());

   t->Scan("long:lat:d", "", "", 10);           // first 10 clumps
   t->Scan("name:prof:J:Mtid:Dgal", "", "", 10);

   // Which profiles are in the file, and their integer codes?
   if (TList *m = (TList *)f->Get("prof_codes"))
      for (auto o : *m) printf("prof %-14s -> code %s\n", o->GetName(), ((TNamed *)o)->GetTitle());

   printf("extended clumps: %lld / %lld\n",
          t->GetEntries("is_extended"), t->GetEntries());
   printf("clumps with an unusable tidal radius: %lld\n", t->GetEntries("!tid_valid"));

   TCanvas *c = new TCanvas("c", "clumps", 1200, 800);
   c->Divide(2, 2);

   c->cd(1); gPad->SetLogx(); gPad->SetLogy();
   t->Draw("J", "J>0");                                    // J-factor distribution

   c->cd(2); gPad->SetLogx(); gPad->SetLogy();
   t->Draw("Mtid:J", "J>0 && Mtid>0", "COLZ");             // Mtid vs J

   c->cd(3);
   t->Draw("lat:long", "", "COLZ");                        // sky map (glat:glon also works)

   c->cd(4); gPad->SetLogy();
   t->Draw("f_aperture", "is_extended");                   // J fraction inside the pixel

   c->SaveAs("clumps_overview.png");
}
