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
#include "TCARFileLoader.h"
#include "TCUtils.h"
#include "TCReadConfig.h"

Double_t* gOff;
Double_t* gGain;
Long64_t gEntries;
Short_t* gElem_1;
Short_t* gElem_2;
Float_t* gRaw_Time_1;
Float_t* gRaw_Time_2;
Bool_t gIsTAPS;
Int_t gMaxTAPS;

void CalibrateMinimize(Int_t nPar);
void CalibrateFit(Int_t nPar, Double_t time_limit, Double_t bad_frac);
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
    //const Bool_t doFit        = kTRUE;
    //const Char_t* treeName    = "CaLib_CB_Time_Tree";
    //const Int_t nPar          = 720;
    //const Char_t* data_off    = "Data.CB.T0";
    //const Char_t* data_gain   = 0;
    //const Char_t* calibration = "LH2_Dec_14";
    //const Double_t time_limit = 0.1;
    //const Double_t bad_frac   = 0.05;
    //gIsTAPS                   = kFALSE;

    // TAPS
    const Bool_t doFit        = kTRUE;
    const Char_t* treeName    = "CaLib_TAPS_Time_Tree";
    const Int_t nPar          = 438;
    const Char_t* data_off    = "Data.TAPS.T0";
    const Char_t* data_gain   = "Data.TAPS.T1";
    const Char_t* calibration = "LH2_Dec_14";
    const Double_t time_limit = 0.1;
    const Double_t bad_frac   = 0.05;
    gIsTAPS                   = kTRUE;

    // read TAPS config from config file
    gMaxTAPS = TCReadConfig::GetReader()->GetConfigInt("TAPS.Elements");

    // read calibration parameters
    gOff = new Double_t[nPar];
    gGain = new Double_t[nPar];
    Int_t set = TCMySQLManager::GetManager()->GetSetForRun(data_off, calibration, run);
    printf("Reading set %d for run %d of data '%s' for calibration '%s'\n", set, run, data_off, calibration);
    TCMySQLManager::GetManager()->ReadParameters(data_off, calibration, set, gOff, nPar);
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
    //    printf("%3d  %f   %f\n", i, gOff[i], gGain[i]);

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

    // update number of entries
    gEntries = n;

    // create histogram
    TH2* hBefore = CreateHisto("before", nPar, gOff);

    // perform calibration
    TStopwatch watch;
    watch.Start();
    if (doFit)
    {
        printf("Calibrating using fitting\n");
        CalibrateFit(nPar, time_limit, bad_frac);
    }
    else
    {
        printf("Calibrating using minimization\n");
        CalibrateMinimize(nPar);
    }
    watch.Stop();
    watch.Print();

    // create histogram
    TH2* hAfter = CreateHisto("after", nPar, gOff);

    // save histograms
    TFile* fout = new TFile(TString::Format("tree_calib_%d.root", run).Data(), "recreate");
    hBefore->Write();
    hAfter->Write();
    delete fout;

    // save parameters
    FILE* fout_a = fopen(TString::Format("tree_calib_%d.dat", run).Data(), "w");
    for (Int_t i = 0; i < nPar; i++)
        fprintf(fout_a, "%lf\n", gOff[i]);
    fclose(fout_a);

    // clean-up
    //delete f;
    //delete tree;
    delete [] gOff;
    delete [] gGain;

    gSystem->Exit(0);
}

//______________________________________________________________________________
void CalibrateMinimize(Int_t nPar)
{
    // Calibrate using minimization.

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
        min.SetVariable(i, TString::Format("Par_%d", i).Data(), gOff[i], 1);

    // perform minimization
    min.Minimize();

    // copy parameters
    for (Int_t i = 0; i < nPar; i++)
        gOff[i] = min.X()[i];
}

//______________________________________________________________________________
void CalibrateFit(Int_t nPar, Double_t time_limit, Double_t bad_frac)
{
    // Calibrate using iterative fitting.

    const Double_t convergence_factor = 0.05;

    // iterate calibrations
    for (Int_t i = 0; i >= 0; i++)
    {
        // create new histogram
        TH2* h = CreateHisto("fit_histo", nPar, gOff);

        // user info
        printf("Calibration iteration %d\n", i+1);

        // loop over elements
        Int_t n = 0;
        Int_t outside_n = 0;
        for (Int_t j = 0; j < nPar; j++)
        {
            // create projection
            TH1* hProj = (TH1*) h->ProjectionX(TString::Format("Proj_%d", j).Data(), j+1, j+1, "e");

            // check for filled histograms
            if (hProj->GetEntries())
            {
                // create fitting function
                TF1* func = new TF1(TString::Format("Func_%d", j).Data(), "gaus", -5, 5);
                func->SetParameter(0, 1);
                func->SetParameter(1, hProj->GetXaxis()->GetBinCenter(hProj->GetMaximumBin()));
                func->SetParameter(2, 0.1);

                // fit histogram
                hProj->GetXaxis()->SetRangeUser(-5, 5);
                hProj->Fit(func, "0Q");
                Double_t mean = func->GetParameter(1);

                // update offset
                gOff[j] = gOff[j] + convergence_factor * mean / gGain[j];
                n++;
                if (TMath::Abs(mean) > time_limit)
                    outside_n++;

                // clean-up
                delete func;
            }

            // clean-up
            delete hProj;
        }

        // clean-up
        delete h;

        // user info
        Double_t outside_frac = (Double_t)outside_n/(Double_t)n;
        printf("Element outside limits: %.1f%%\n", outside_frac*100.);
        if (outside_frac < bad_frac)
            break;
    }
}

//______________________________________________________________________________
Double_t CostFunction(const Double_t* par)
{
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

