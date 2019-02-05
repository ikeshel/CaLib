/*************************************************************************
 * Author: Dominik Werthmueller, Thomas Strub
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// CheckScR.C                                                           //
//                                                                      //
// Compare the number of scaler reads in a file to the value stored in  //
// the calibration database.
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
void CheckScR(const Char_t* loc)
{
    // Main method.

    // load CaLib
    gSystem->Load("libCaLib.so");
    TCMySQLManager::GetManager();

    // try to get directory content
    TSystemDirectory dir("rawdir", loc);
    TList* list = dir.GetListOfFiles();
    if (!list)
    {
        Error("CheckScR", "'%s' is not a directory!", loc);
        gSystem->Exit(1);
    }

    // sort files
    list->Sort();

    // loop over directory content
    TIter next(list);
    TSystemFile* f;
    while ((f = (TSystemFile*)next()))
    {
        // look for ROOT files
        TString str(f->GetName());
        if (str.EndsWith(".root"))
        {
            Int_t n_hist = 0;
            Int_t n_tree = 0;
            Int_t n_db = 0;

            // full path
            TString full = TString::Format("%s/%s", loc, f->GetName());

            // extract run number
            TString tmp(f->GetName());
            tmp.ReplaceAll("ARHistograms_", "");
            tmp.ReplaceAll("CBTaggTAPS_", "");
            tmp.ReplaceAll("CBTaggTAPSPed_", "");
            tmp.ReplaceAll("TaggEff_", "");
            tmp.ReplaceAll(".root", "");
            Int_t run = tmp.Atoi();

            // open the file
            TFile file(full.Data());

            // try to load EventInfo histogram
            TH1* h = (TH1*) file.Get("EventInfo");
            if (h)
                n_hist = h->GetBinContent(TCConfig::kNScREventHBin);

            // try to load ScalerEvents tree
            TTree* t = (TTree*) file.Get("ScalerEvents");
            if (t)
                n_tree = t->GetEntries();

            // read from db
            n_db = TCMySQLManager::GetManager()->GetRunNScR(run);

            // debug
            //printf("Run %5d  EventInfo: %4d  ScalerEvents: %4d  Database: %4d\n",
            //       run, n_hist, n_tree, n_db);

            // compare numbers
            if (n_hist != n_db)
                Error("CheckScR", "Run %5d : Different number of scaler reads in EventInfo and DB! (%4d, %4d, diff = %d)",
                      run, n_hist, n_db, n_db-n_hist);
        }
    }

    gSystem->Exit(0);
}

