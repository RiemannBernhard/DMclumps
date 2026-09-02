// Print everything a CLUMPY-derived ROOT file knows about itself.
//   root -l 'describe.C("clumps.root")'          -- full dictionary
//   root -l 'describe.C("clumps.root","J_tot")'  -- one branch
void describe(const char *fname = "clumps.root", const char *branch = "")
{
   TFile *f = TFile::Open(fname);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", fname); return; }

   if (branch && branch[0]) {                       // single-branch lookup
      TTree *t = nullptr; f->GetObject("clumps", t);
      if (t) if (TList *ui = t->GetUserInfo())
         if (TObject *o = ui->FindObject(branch)) {
            printf("\n%s\n  %s\n\n", o->GetName(), ((TNamed *)o)->GetTitle());
            f->Close(); return;
         }
      printf("no dictionary entry for '%s'\n", branch);
      f->Close(); return;
   }

   if (TMacro *rd = (TMacro *)f->Get("README")) rd->Print();
   else printf("this file carries no README\n");

   if (TList *pc = (TList *)f->Get("prof_codes")) {
      printf("\nPROFILE CODES\n");
      for (auto o : *pc) printf("  %-14s -> %s\n", o->GetName(), ((TNamed *)o)->GetTitle());
   }
   f->Close();
}
