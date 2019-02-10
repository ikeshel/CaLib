/*************************************************************************
 * Author: Irakli Keshelashvili, Dominik Werthmueller, Thomas Strub
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibCBTimeWalk                                                    //
//                                                                      //
// Calibration module for the CB time walk.                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCLine.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TSystem.h"

#include "TCCalibCBTimeWalk.h"
#include "TCConfig.h"
#include "TCMySQLManager.h"
#include "TCFileManager.h"
#include "TCReadConfig.h"
#include "TCUtils.h"

ClassImp(TCCalibCBTimeWalk)

//______________________________________________________________________________
TCCalibCBTimeWalk::TCCalibCBTimeWalk()
    : TCCalib("CB.TimeWalk", "CB time walk calibration", "Data.CB.Walk.Par0", TCConfig::kMaxCB)
{
    // Empty constructor.

    // init members
    fFileManager = 0;
    fPar0 = 0;
    fPar1 = 0;
    fPar2 = 0;
    fPar3 = 0;
    fGFitPoints = 0;
    fTimeProj = 0;
    fLine = 0;
    fDelay = 0;
    fUseEnergyWeight = kTRUE;
    fUsePointDensityWeight = kTRUE;
    fWalkType = kDefault;
}

//______________________________________________________________________________
TCCalibCBTimeWalk::~TCCalibCBTimeWalk()
{
    // Destructor.

    if (fFileManager) delete fFileManager;
    if (fPar0) delete [] fPar0;
    if (fPar1) delete [] fPar1;
    if (fPar2) delete [] fPar2;
    if (fPar3) delete [] fPar3;
    if (fGFitPoints) delete fGFitPoints;
    if (fTimeProj) delete fTimeProj;
    if (fLine) delete fLine;
}

//______________________________________________________________________________
void TCCalibCBTimeWalk::Init()
{
    // Init the module.

    // init members
    fFileManager = new TCFileManager(fData.Data(), fCalibration.Data(), fNset, fSet);
    fPar0 = new Double_t[fNelem];
    fPar1 = new Double_t[fNelem];
    fPar2 = new Double_t[fNelem];
    fPar3 = new Double_t[fNelem];
    fGFitPoints = 0;
    fTimeProj = 0;
    fLine =  new TCLine();
    fDelay = 0;

    // configure line
    fLine->SetLineColor(4);
    fLine->SetLineWidth(3);

    // get histogram name
    if (!TCReadConfig::GetReader()->GetConfig("CB.TimeWalk.Histo.Fit.Name"))
    {
        Error("Init", "Histogram name was not found in configuration!");
        return;
    }
    else fHistoName = *TCReadConfig::GetReader()->GetConfig("CB.TimeWalk.Histo.Fit.Name");

    // get projection fit display delay
    fDelay = TCReadConfig::GetReader()->GetConfigInt("CB.TimeWalk.Fit.Delay");

    // init weigths
    fUseEnergyWeight = kTRUE;
    fUsePointDensityWeight = kTRUE;

    // get correction type
    TString* type = TCReadConfig::GetReader()->GetConfig("CB.TimeWalk.Type");
    if (!type)
    {
        fWalkType = kDefault;
        Info("Init", "Using default walk correction");
    }
    else
    {
        TString t(*type);
        t.ToLower();
        if (t == "default")
        {
            fWalkType = kDefault;
            Info("Init", "Using default walk correction");
        }
        else if (t == "strub")
        {
            fWalkType = kStrub;
            Info("Init", "Using Strub walk correction");
        }
        else
        {
            fWalkType = kDefault;
            Info("Init", "Using default walk correction");
        }
    }

    // read old parameters (only from first set)
    TCMySQLManager::GetManager()->ReadParameters("Data.CB.Walk.Par0", fCalibration.Data(), fSet[0], fPar0, fNelem);
    TCMySQLManager::GetManager()->ReadParameters("Data.CB.Walk.Par1", fCalibration.Data(), fSet[0], fPar1, fNelem);
    TCMySQLManager::GetManager()->ReadParameters("Data.CB.Walk.Par2", fCalibration.Data(), fSet[0], fPar2, fNelem);
    TCMySQLManager::GetManager()->ReadParameters("Data.CB.Walk.Par3", fCalibration.Data(), fSet[0], fPar3, fNelem);

    // draw main histogram
    fCanvasFit->Divide(1, 2, 0.001, 0.001);
    fCanvasFit->cd(1)->SetLogz();

    // draw the overview histogram
    fCanvasResult->cd();
}

//______________________________________________________________________________
void TCCalibCBTimeWalk::Fit(Int_t elem)
{
    // Perform the fit of the element 'elem'.

    Char_t tmp[256];

    // get configuration
    Double_t lowLimit, highLimit;
    TCReadConfig::GetReader()->GetConfigDoubleDouble("CB.TimeWalk.Histo.Fit.Xaxis.Range", &lowLimit, &highLimit);

    // create histogram name
    sprintf(tmp, "%s_%03d", fHistoName.Data(), elem);

    // get histogram
    if (!fIsReFit)
    {
        // delete old histogram
        if (fMainHisto) delete fMainHisto;

         // get new
        fMainHisto = fFileManager->GetHistogram(tmp);
    }

    if (!fMainHisto)
    {
        Error("Fit", "Main histogram does not exist!");
        return;
    }

    // draw main histogram
    if (!fIsReFit)
    {
        if (fMainHisto->GetEntries() > 0) fCanvasFit->cd(1)->SetLogz(1);
        else fCanvasFit->cd(1)->SetLogz(0);
        TCUtils::FormatHistogram(fMainHisto, "CB.TimeWalk.Histo.Fit");
        fMainHisto->Draw("colz");
        fCanvasFit->Update();
    }

    // check for sufficient statistics
    if (fMainHisto->GetEntries() < 1000 && !TCUtils::IsCBHole(elem))
    {
        Error("Fit", "Not enough statistics!");
        return;
    }

    // copy old points
    TGraphErrors oldp;
    if (fIsReFit && fGFitPoints) oldp = TGraphErrors(*fGFitPoints);

    // create new graph (fit data)
    if (fGFitPoints) delete fGFitPoints;
    fGFitPoints = new TGraphErrors();
    fGFitPoints->SetMarkerStyle(20);
    fGFitPoints->SetMarkerColor(kBlack);

    // prepare stuff for adding
    Double_t low_e = 0;
    Int_t added = 0;
    Double_t added_e = 0;
    Double_t added_w = 0;

    // get bins for fitting range
    Int_t startBin = fMainHisto->GetXaxis()->FindBin(lowLimit);
    Int_t endBin = fMainHisto->GetXaxis()->FindBin(highLimit);

    // loop over energy bins
    for (Int_t i = startBin; i <= endBin; i++)
    {
        // create time projection
        sprintf(tmp, "ProjTime_%d_%d", elem, i);
        TH1* proj = (TH1D*) ((TH2*) fMainHisto)->ProjectionY(tmp, i, i, "e");

        // skip empty projections
        if (proj->GetEntries() == 0)
        {
            delete proj;
            continue;
        }

        // first loop (after fit)
        if (added == 0)
        {
            low_e = fMainHisto->GetXaxis()->GetBinLowEdge(i);
            if (fTimeProj) delete fTimeProj;
            fTimeProj = (TH1*) proj->Clone("TimeProjection");

            // reset values
            added_e = 0;
            added_w = 0;
        }
        else
        {
            fTimeProj->Add(proj);
        }

        // add up bin contribution
        Double_t weight = fUseEnergyWeight ? proj->GetEntries() : 1.;
        added_w += weight;
        added_e += fMainHisto->GetXaxis()->GetBinCenter(i) * weight;
        added++;

        // delete projection
        delete proj;

        // calc energy interval
        Double_t e_int = fMainHisto->GetXaxis()->GetBinUpEdge(i) - low_e;

        // minimum energy interval
        if (e_int < 0.5) continue;

        // check if projection has enough entries
        if (fTimeProj->GetEntries() < 300 && i < endBin)
            continue;

        // calculate mean energy
        Double_t energy = added_e/added_w;

        //
        // fit time projection
        //

        // create fitting function
        if (fFitFunc) delete fFitFunc;
        sprintf(tmp, "fWalkProfile_%i", elem);
        fFitFunc = new TF1(tmp, "gaus(0)");
        fFitFunc->SetLineColor(kBlue);

        // prepare fitting function
        Int_t maxbin = fTimeProj->GetMaximumBin();
        Double_t peak = fTimeProj->GetBinCenter(maxbin);
        Double_t max = fTimeProj->GetMaximum();

        // find old (changed) peak pos
        Bool_t skip = kFALSE;
        if (fIsReFit)
        {
            // init
            skip = kTRUE;

            // loop over old points
            for (Int_t j = 0; j < oldp.GetN(); j++)
            {
                // find point
                if (low_e < oldp.GetX()[j] && oldp.GetX()[j] < low_e+e_int)
                {
                    skip = kFALSE;
                    peak = oldp.GetY()[j];
                    break;
                }
            }
        }

        // skip if no old (deleted) point found
        if (skip)
        {
            added = 0;
            continue;
        }

        // prepare fit
        fFitFunc->SetRange(peak - 20, peak + 20);
        fFitFunc->SetParameters(max, peak, 1.);
        fFitFunc->SetParLimits(0, max*0.5, max*1.5); // peak height
        fFitFunc->SetParLimits(1, peak - 5, peak + 5); // peak position
        fFitFunc->SetParLimits(2, 0.5, 20.0); // sigma

        // perform fit
        fTimeProj->Fit(fFitFunc, "RBQ0");

        // get parameters
        Double_t mean = fFitFunc->GetParameter(1);
        Double_t error = fFitFunc->GetParError(1);

        // format line
        fLine->SetPos(mean);

        // check fit error
        if (error < 1.)
        {
            Int_t n = fGFitPoints->GetN();
            fGFitPoints->SetPoint(n, energy, mean);
            fGFitPoints->SetPointError(n, 0., error);
        }

        // plot projection fit
        if (fDelay > 0)
        {
            fCanvasFit->cd(2);
            fTimeProj->GetXaxis()->SetRangeUser(mean - 30, mean + 30);
            fTimeProj->Draw("hist");
            fFitFunc->Draw("same");
            fLine->Draw();
            fCanvasFit->Update();
            gSystem->Sleep(fDelay);
        }

        // reset value
        added = 0;

    } // for: loop over energy bins

    //
    // fit profile
    //

    // check for enough points
    if (fGFitPoints->GetN() < 3)
    {
        Error("Fit", "Not enough fit points!");
        return;
    }

    // create fitting function
    sprintf(tmp, "fTWalk_%d", elem);
    if (fFitFunc) delete fFitFunc;
    if (fWalkType == kDefault)
        fFitFunc = new TF1(tmp, "[0] + [1] / TMath::Power(x + [2], [3])", lowLimit, highLimit);
    else if (fWalkType == kStrub)
        fFitFunc = new TF1(tmp, "[0] + [1] / TMath::Power(x - 1, [3]) + x*[2]", lowLimit, highLimit);
    fFitFunc->SetLineColor(kBlue);
    fFitFunc->SetNpx(2000);

    // weight point errors by point density
    if (fUsePointDensityWeight)
    {
        for (Int_t i = 0; i < fGFitPoints->GetN(); i++)
        {
            Double_t* x = fGFitPoints->GetX();
            //Double_t* y = fGFitPoints->GetY();

            Double_t diff = 0.;

            // get distance to previous point (use only x-coord)
            if (i > 0)
            {
                Double_t xdiff = x[i]-x[i-1];
                Double_t ydiff = 0;//y[i]-y[i-1];
                diff  += TMath::Sqrt(xdiff*xdiff + ydiff*ydiff);
            }

            // get distance to next point
            if (i < fGFitPoints->GetN()-1)
            {
                Double_t xdiff = x[i+1]-x[i];
                Double_t ydiff = 0;//y[i+1]-y[i];
                diff  += TMath::Sqrt(xdiff*xdiff + ydiff*ydiff);
            }

            // get mean distance
            if (i > 0 && i < fGFitPoints->GetN()-1)
                diff /= 2.;

            // set new error = old error / sqrt(dist)
            fGFitPoints->SetPointError(i, fGFitPoints->GetEX()[i], 1./TMath::Sqrt(diff) * fGFitPoints->GetEY()[i]);
        }
    }


    // prepare fitting function
    fFitFunc->SetParameters(-50, 60, 0.2, 0.3);
    fFitFunc->SetParLimits(1, -5000, 5000);
    if (fWalkType == kDefault)
        fFitFunc->SetParLimits(2, -1, 10);
    fFitFunc->SetParLimits(3, 0, 1);

    // perform fit
    for (Int_t i = 0; i < 10; i++)
        if (!fGFitPoints->Fit(fFitFunc, "RB0Q")) break;

    // read parameters
    fPar0[elem] = fFitFunc->GetParameter(0);
    fPar1[elem] = fFitFunc->GetParameter(1);
    fPar2[elem] = fFitFunc->GetParameter(2);
    fPar3[elem] = fFitFunc->GetParameter(3);

    // draw energy projection and fit
    fFitFunc->SetRange(-fFitFunc->GetParameter(2), 1000);
    fCanvasResult->cd();
    fGFitPoints->Draw("ap");
    fGFitPoints->GetXaxis()->SetLimits(lowLimit, highLimit);
    fFitFunc->Draw("same");
    fCanvasResult->Update();

    fCanvasFit->cd(1);
    fFitFunc->Draw("same");
    fCanvasFit->Update();
}

//______________________________________________________________________________
void TCCalibCBTimeWalk::Calculate(Int_t elem)
{
    // Calculate the new value of the element 'elem'.

    Bool_t noval = kFALSE;

    // check no fits
    if (fPar0[elem] == 0 && fPar1[elem] == 0 &&
        fPar2[elem] == 0 && fPar3[elem] == 0) noval = kTRUE;

    // user information
    printf("Element: %03d    Par0: %12.8f    "
           "Par1: %12.8f    Par2: %12.8f    Par3: %12.8f",
           elem, fPar0[elem], fPar1[elem], fPar2[elem], fPar3[elem]);
    if (noval) printf("    -> no fit");
    if (TCUtils::IsCBHole(elem)) printf(" (hole)");
    printf("\n");
}

//______________________________________________________________________________
void TCCalibCBTimeWalk::PrintValues()
{
    // Print out the old and new values for all elements.

    // loop over elements
    for (Int_t i = 0; i < fNelem; i++)
    {
        printf("Element: %03d    Par0: %12.8f    "
               "Par1: %12.8f    Par2: %12.8f    Par3: %12.8f\n",
               i, fPar0[i], fPar1[i], fPar2[i], fPar3[i]);
    }
}

//______________________________________________________________________________
void TCCalibCBTimeWalk::WriteValues()
{
    // Write the obtained calibration values to the database.

    // write values to database
    for (Int_t i = 0; i < fNset; i++)
    {
        TCMySQLManager::GetManager()->WriteParameters("Data.CB.Walk.Par0", fCalibration.Data(), fSet[i], fPar0, fNelem);
        TCMySQLManager::GetManager()->WriteParameters("Data.CB.Walk.Par1", fCalibration.Data(), fSet[i], fPar1, fNelem);
        TCMySQLManager::GetManager()->WriteParameters("Data.CB.Walk.Par2", fCalibration.Data(), fSet[i], fPar2, fNelem);
        TCMySQLManager::GetManager()->WriteParameters("Data.CB.Walk.Par3", fCalibration.Data(), fSet[i], fPar3, fNelem);
    }
}

