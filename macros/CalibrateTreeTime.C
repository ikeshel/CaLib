/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// CalibrateTreeTime.C                                                  //
//                                                                      //
// Calibrate detector timings using an event tree.                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include "TCMySQLManager.h"

Double_t* gGain;
Long64_t gEntries;
Short_t* gElem_1;
Short_t* gElem_2;
Float_t* gRaw_Time_1;
Float_t* gRaw_Time_2;
Bool_t gIsTAPS;
Int_t gMaxTAPS;

Double_t CostFunction(const Double_t* par);
TH2* CreateHisto(const Char_t* name, Int_t nPar, const Double_t* par);

//______________________________________________________________________________
void CalibrateTreeTime(Int_t run)
{
    // Main method.

    //
    // configuration
    //

    // CB
    //const Char_t* treeName    = "CaLib_CB_Time_Tree";
    //const Int_t nPar          = 720;
    //const Char_t* data_off    = "Data.CB.T0";
    //const Char_t* data_gain   = 0;
    //const Char_t* calibration = "LH2_Dec_14";
    gIsTAPS                   = kFALSE;

    // TAPS
    const Char_t* treeName    = "CaLib_TAPS_Time_Tree";
    const Int_t nPar          = 438;
    const Char_t* data_off    = "Data.TAPS.T0";
    const Char_t* data_gain   = "Data.TAPS.T1";
    const Char_t* calibration = "LH2_Dec_14";
    gIsTAPS                   = kTRUE;

    // read TAPS config from config file
    gMaxTAPS = TCReadConfig::GetReader()->GetConfigInt("TAPS.Elements");

    // read calibration parameters
    Double_t par[nPar];
    gGain = new Double_t[nPar];
    Int_t set = TCMySQLManager::GetManager()->GetSetForRun(data_off, calibration, run);
    printf("Reading set %d for run %d of data '%s' for calibration '%s'\n", set, run, data_off, calibration);
    TCMySQLManager::GetManager()->ReadParameters(data_off, calibration, set, par, nPar);
    if (data_gain)
    {
        TCMySQLManager::GetManager()->ReadParameters(data_gain, calibration, set, gGain, nPar);
    }
    else
    {
        for (Int_t i = 0; i < nPar; i++)
            gGain[i] = 0.11771;
    }

    // debug
    //for (Int_t i = 0; i < nPar; i++)
    //    printf("%3d  %f   %f\n", i, par[i], gGain[i]);

    // open the file
    Int_t r[1] = { run };
    TCARFileLoader file_loader(1, r);
    if (!file_loader.LoadFiles())
        gSystem->Exit(1);
    TFile* f = file_loader.GetFiles()[0];
    if (!f)
        gSystem->Exit(1);

    // load tree
    TTree* tree = (TTree*) f->Get(treeName);
    gEntries = tree->GetEntries();

    // reserve memory
    gElem_1 = new Short_t[gEntries];
    gElem_2 = new Short_t[gEntries];
    gRaw_Time_1 = new Float_t[gEntries];
    gRaw_Time_2 = new Float_t[gEntries];

    // create tree reader
    TTreeReader reader(tree);

    // configure reader
    TTreeReaderValue<Short_t> elem_1(reader, "elem_1");
    TTreeReaderValue<Short_t> elem_2(reader, "elem_2");
    TTreeReaderValue<Float_t> raw_time_1(reader, "raw_time_1");
    TTreeReaderValue<Float_t> raw_time_2(reader, "raw_time_2");

    // read tree
    printf("Reading tree\n");
    Long64_t n = 0;
    while (reader.Next())
    {
        // user info
        if (n && n % 100000 == 0)
            printf("%lld events processed (%.1f%%)\n", n, 100*(Double_t)n/(Double_t)gEntries);

        // skip PWO
        if (gIsTAPS)
        {
            if (TCUtils::IsTAPSPWO(*elem_1, gMaxTAPS) ||
                TCUtils::IsTAPSPWO(*elem_2, gMaxTAPS)) continue;
        }

        gElem_1[n] = *elem_1;
        gElem_2[n] = *elem_2;
        gRaw_Time_1[n] = *raw_time_1;
        gRaw_Time_2[n] = *raw_time_2;
        n++;
    }

    // create minimizer
    ROOT::Minuit2::Minuit2Minimizer min(ROOT::Minuit2::kMigrad);
    min.SetPrintLevel(1);
    min.SetMaxFunctionCalls(1000000);
    min.SetMaxIterations(100000);
    min.SetTolerance(0.01);

    // set minimizer function
    ROOT::Math::Functor functor(&CostFunction, nPar);
    min.SetFunction(functor);

    // set fit parameters
    for (Int_t i = 0; i < nPar; i++)
        min.SetVariable(i, TString::Format("Par_%d", i).Data(), par[i], 1);

    // perform minimization
    printf("Minimizing\n");
    TStopwatch watch;
    watch.Start();
    min.Minimize();
    watch.Stop();
    watch.Print();

    // create histogram
    TH2* hBefore = CreateHisto("before", nPar, par);
    TH2* hAfter = CreateHisto("after", nPar, min.X());

    // save histograms
    TFile* fout = new TFile(TString::Format("tree_calib_%d.root", run).Data(), "recreate");
    hBefore->Write();
    hAfter->Write();
    delete fout;

    // save parameters
    FILE* fout_a = fopen(TString::Format("tree_calib_%d.dat", run).Data(), "w");
    for (Int_t i = 0; i < nPar; i++)
        fprintf(fout_a, "%lf\n", par[i]);
    fclose(fout_a);

    // clean-up
    //delete f;
    //delete tree;
    delete [] gGain;

    gSystem->Exit(0);
}

//______________________________________________________________________________
Double_t CostFunction(const Double_t* par)
{
    // read tree
    Double_t sum = 0;
    for (Long64_t i = 0; i < gEntries; i++)
    {
        const Double_t time_1 = gGain[gElem_1[i]] * (gRaw_Time_1[i] - par[gElem_1[i]]);
        const Double_t time_2 = gGain[gElem_2[i]] * (gRaw_Time_2[i] - par[gElem_2[i]]);
        sum += TMath::Abs(time_1 - time_2);
        sum += TMath::Abs(time_2 - time_1);
    }

    return sum;
}

//______________________________________________________________________________
TH2* CreateHisto(const Char_t* name, Int_t nPar, const Double_t* par)
{
    TH2* h;
    if (gIsTAPS)
        h = new TH2F(name, name, 4000, -100, 100, nPar, 0, nPar);
    else
        h = new TH2F(name, name, 2000, -100, 100, nPar, 0, nPar);

    // fill histo
    for (Long64_t i = 0; i < gEntries; i++)
    {
        const Double_t time_1 = gGain[gElem_1[i]] * (gRaw_Time_1[i] - par[gElem_1[i]]);
        const Double_t time_2 = gGain[gElem_2[i]] * (gRaw_Time_2[i] - par[gElem_2[i]]);
        h->Fill(time_1 - time_2, gElem_1[i]);
        h->Fill(time_2 - time_1, gElem_2[i]);
    }

    return h;
}

