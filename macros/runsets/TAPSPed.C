/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAPSPed.C                                                            //
//                                                                      //
// Make run sets depending on the stability in time of a calibration.   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TList.h"
#include "TCReadARCalib.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TSystem.h"
#include "TCMySQLManager.h"
#include "TLine.h"
#include "TGraph.h"
#include "TOSUtils.h"
#endif

class TCReadARCalib;
TCReadARCalib* gReadAR;
TList* gFiles;

//______________________________________________________________________________
void CheckPedestals(const Char_t* loc, const Char_t* filePre)
{
    // Some old pedestal checking method.

    Char_t t[256];

    // number of runs
    Int_t nRuns = gFiles->GetSize();

    // number of channels
    Int_t nCh = gReadAR->GetNelements();

    // create arrays
    Double_t** pedPos = new Double_t*[nCh];
    Double_t* runNumbersD = new Double_t[nRuns];
    for (Int_t i = 0; i < nCh; i++) pedPos[i] = new Double_t[nRuns];

    // open the output files
    TFile* fROOTout = new TFile("/tmp/ped.root", "RECREATE");

    // create directories
    for (Int_t i = 0; i < nCh; i++)
    {
        sprintf(t, "%03d_%s", i, gReadAR->GetElement(i)->GetADC());
        fROOTout->mkdir(t);
    }

    TF1* func = new TF1("gausfunc", "gaus", 0 , 200);

    // loop over runs
    for (Int_t i = 0; i < nRuns; i++)
    {
        // get the file
        TFile* f = (TFile*) gFiles->At(i);

        // extract run number
        Int_t runNumber;
        sprintf(t, "%s/%s_%%d.root", loc, filePre);
        strcpy(t, TOSUtils::ExpandPath(t));
        sscanf(f->GetName(), t, &runNumber);
        runNumbersD[i] = (Double_t)runNumber;

        printf("Processing run %d (%d/%d)\n", runNumber, i+1, nRuns);

        fROOTout->cd();

        // loop over ADCs
        for (Int_t j = 0; j < nCh; j++)
        {
            // load histogram
            sprintf(t, "%03d_%s", j, gReadAR->GetElement(j)->GetADC());
            fROOTout->cd(t);
            sprintf(t, "ADC%s", gReadAR->GetElement(j)->GetADC());
            TH1D* h = (TH1D*) f->Get(t);
            if (!h) continue;

            // fit gaussian to pedestal
            h->GetXaxis()->SetRange(90, 107);
            //else h->GetXaxis()->SetRange(90, 115);
            //if (j == 27112) h->GetXaxis()->SetRange(60, 110);
            //else h->GetXaxis()->SetRange(60, 99);
            //h->GetXaxis()->SetRangeUser(40, 120);

            Double_t maxPos = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
            Double_t MaxCont = h->GetBinContent(h->GetMaximumBin());
            if(h->GetEntries()<5)
            {
                maxPos = 0;
                MaxCont = 0;
            }
            else if(h->GetEntries()>50)
            {
                // find first big jump
                for(int u=10; u<200; u++)
                {
                    Double_t BinMin = 10;
                    Double_t BinDiff = 5;
                    if(h->GetEntries()<150)
                    {
                        BinMin = 2;
                        BinDiff =1;

                    }
                    else if(h->GetEntries()<2000)
                    {
                        BinMin = 3;
                        BinDiff =2;
                    }

                    if(h->GetBinContent(u) > BinMin &&
                       h->GetBinContent(u) > BinDiff*h->GetBinContent(u-3) &&
                       h->GetBinContent(u) > BinDiff*h->GetBinContent(u+3))
                    {
                        //Printf("%i", i);
                        Double_t fTotMaxPos = u;
                        maxPos = h->GetBinCenter( fTotMaxPos );
                        MaxCont= h->GetBinContent(fTotMaxPos );
                        break;
                    }
                }
                if(TMath::Abs(maxPos-(h->GetXaxis()->GetBinCenter(h->GetMaximumBin())))<5.) maxPos = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
                //Printf("%5.3f", maxPos);
                h->GetXaxis()->SetRangeUser(40, 150);
                func->SetParameters(MaxCont, maxPos, 0.1);
                func->SetRange(maxPos - 3, maxPos + 3);
                h->Fit(func, "RBQ");
                maxPos = func->GetParameter(1);
            }
            // save position in file and memory
            pedPos[j][i] = maxPos;

            sprintf(t, "Run_%d", runNumber);
            TCanvas* c = new TCanvas(t, t);
            h->Draw();

            TLine* tline = new TLine(maxPos, 0, maxPos, func->GetParameter(0));
            tline->SetLineColor(kRed);
            tline->SetLineWidth(2);
            tline->Draw();

            c->Write(c->GetName(), TObject::kOverwrite);

            delete h;
            delete c;
            delete tline;
        }
    }

    // create pedestal evolution graphs
    fROOTout->cd();

    // loop over channels
    for (Int_t j = 0; j < nCh; j++)
    {
        printf("Creating pedestal graph for channel %d\n", j);

        TGraph* g = new TGraph(nRuns, runNumbersD, pedPos[j]);
        sprintf(t, "Overview_%03d_%s", j, gReadAR->GetElement(j)->GetADC());
        g->SetName(t);
        g->SetTitle(t);
        g->GetYaxis()->SetRangeUser(80, 120);
        g->Write(g->GetName());
        //TString nn(g->GetName());
        //if (nn.Contains("_29")) g->Write(g->GetName(), TObject::kOverwrite);
        //for(Int_t i = 0; i<g->GetN(); i++)
        //{
        //    if(TMath::Abs((g->GetY()[i]-g->GetY()[i-1]))>5. && g->GetY()[i]>0. && g->GetY()[i-1]>0 && i>0)
        //        printf("-%i \n %i", g->GetX()[i-1],g->GetX()[i]);
        //        //cout << "-"<< g->GetX()[i-1] << endl << g->GetX()[i];
        //}
        //printf("\n");

        delete g;
    }

    printf("Saving output file\n");

    delete fROOTout;

    // cleanup
    for (Int_t i = 0; i < nCh; i++) delete [] pedPos[i];
    delete [] pedPos;
    delete [] runNumbersD;
}

//______________________________________________________________________________
void TAPSPed()
{
    // Main method.

    Char_t tmp[256];

    // load CaLib
    gSystem->Load("libCaLib.so");

    // general configuration
    Bool_t watch = kFALSE;
    //const Char_t* data = "Data.TAPS.LG.E0";
    //const Char_t* data = "Data.TAPS.SG.E0";
    const Char_t* data = "Data.Veto.E0";
    const Char_t* elemDesc = "Element:";
    //const Char_t* elemDesc = "TAPSSG:";
    Double_t yMin = 110;
    Double_t yMax = 160;
    const Char_t calibration[] = "LH2_May_18";
    const Char_t* fLoc = "$HOME/loc/presort/data/May_18";
    //const Char_t* fAR = "/home/werthm/src/ROOT/acqu/acqu_user/data/May_18/TAPS/BaF2_PWO.dat";
    const Char_t* fAR = "/home/werthm/src/ROOT/acqu/acqu_user/data/May_18/TAPS/Veto.dat";
    const Char_t* filePre = "ARHistograms_CBTaggTAPS";

    // read the calibration file with the correct element identifier
    gReadAR = new TCReadARCalib(fAR, kFALSE, elemDesc);

    // user info
    printf("Found %d elements in '%s'\n", gReadAR->GetNelements(), fAR);

    // get number of sets
    Int_t nSets = TCMySQLManager::GetManager()->GetNsets(data, calibration);

    // file array
    gFiles = new TList();
    gFiles->SetOwner(kTRUE);

    // loop over sets
    for (Int_t i = 0; i < nSets; i++)
    {
        // get runs of set
        Int_t nRuns;
        Int_t* runs = TCMySQLManager::GetManager()->GetRunsOfSet(data, calibration, i, &nRuns);

        // loop over runs
        for (Int_t j = 0; j < nRuns; j++)
        {
            // load ROOT file
            sprintf(tmp, "%s/%s_%d.root", fLoc, filePre, runs[j]);
            strcpy(tmp, TOSUtils::ExpandPath(tmp));
            TFile* f = new TFile(tmp);

            // check file
            if (!f) continue;
            if (f->IsZombie()) continue;

            // save file
            gFiles->Add(f);
        }

        // clean-up
        delete runs;
    }

    // check pedestals
    CheckPedestals(fLoc, filePre);

    printf("%d runs analyzed.\n", gFiles->GetSize());

    // clean-up
    delete gFiles;

    gSystem->Exit(0);
}

