/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibVetoEnergy                                                    //
//                                                                      //
// Calibration module for the Veto energy.                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TH2.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"

#include "TCCalibVetoEnergy.h"
#include "TCLine.h"
#include "TCReadConfig.h"
#include "TCFileManager.h"
#include "TCMySQLManager.h"
#include "TCUtils.h"

ClassImp(TCCalibVetoEnergy)

//______________________________________________________________________________
TCCalibVetoEnergy::TCCalibVetoEnergy()
    : TCCalib("Veto.Energy", "Veto energy calibration", "Data.Veto.E1",
              TCReadConfig::GetReader()->GetConfigInt("Veto.Elements"))
{
    // Empty constructor.

    // init members
    fFileManager = 0;
    fPeak = 0;
    fPeakMC = 0;
    fLine = 0;
    fLineMC = 0;
    fMCHisto = 0;
    fFitHistoMC = 0;
    fFitFuncMC = 0;
    fMCFile = 0;
}

//______________________________________________________________________________
TCCalibVetoEnergy::~TCCalibVetoEnergy()
{
    // Destructor.

    if (fFileManager) delete fFileManager;
    if (fLine) delete fLine;
    if (fLineMC) delete fLineMC;
    if (fMCHisto) delete fMCHisto;
    if (fFitHistoMC) delete fFitHistoMC;
    if (fFitFuncMC) delete fFitFuncMC;
    if (fMCFile) delete fMCFile;
}

//______________________________________________________________________________
void TCCalibVetoEnergy::Init()
{
    // Init the module.

    // init members
    fFileManager = new TCFileManager(fData, fCalibration.Data(), fNset, fSet);
    fPeak = 0;
    fPeakMC = 0;
    fLine = new TCLine();
    fLineMC = new TCLine();
    fMCHisto = 0;
    fMCFile = 0;

    // configure line
    fLine->SetLineColor(4);
    fLine->SetLineWidth(3);
    fLineMC->SetLineColor(4);
    fLineMC->SetLineWidth(3);

    // get histogram name
    if (!TCReadConfig::GetReader()->GetConfig("Veto.Energy.Histo.Fit.Name"))
    {
        Error("Init", "Histogram name was not found in configuration!");
        return;
    }
    else fHistoName = *TCReadConfig::GetReader()->GetConfig("Veto.Energy.Histo.Fit.Name");

    // get MC histogram file
    TString fileMC;
    if (!TCReadConfig::GetReader()->GetConfig("Veto.Energy.MC.File"))
    {
        Error("Init", "MC file name was not found in configuration!");
        return;
    }
    else fileMC = *TCReadConfig::GetReader()->GetConfig("Veto.Energy.MC.File");

    // read old parameters (only from first set)
    TCMySQLManager::GetManager()->ReadParameters(fData, fCalibration.Data(), fSet[0], fOldVal, fNelem);

    // copy to new parameters
    for (Int_t i = 0; i < fNelem; i++) fNewVal[i] = fOldVal[i];

    // create the overview histogram
    fOverviewHisto = new TH1F("Overview", ";Element;proton peak position [MeV]", fNelem, 0, fNelem);
    fOverviewHisto->SetMarkerStyle(2);
    fOverviewHisto->SetMarkerColor(4);

    // draw main histogram
    fCanvasFit->Divide(2, 2, 0.001, 0.001);

    // open the MC file
    fMCFile = new TFile(fileMC.Data());
    if (!fMCFile)
    {
        Error("Init", "Could not open MC file!");
        return;
    }
    if (fMCFile->IsZombie())
    {
        Error("Init", "Could not open MC file!");
        return;
    }

    // draw the overview histogram
    fCanvasResult->cd();
    TCUtils::FormatHistogram(fOverviewHisto, "Veto.Energy.Histo.Overview");
    fOverviewHisto->Draw("P");
}

//______________________________________________________________________________
void TCCalibVetoEnergy::FitSlice(TH2* h, Int_t elem, Bool_t isMC)
{
    // Fit the energy slice of the dE vs E histogram 'h' of the element 'elem'.

    Char_t tmp[256];

    // get configuration
    Double_t lowLimit, highLimit;
    TCReadConfig::GetReader()->GetConfigDoubleDouble("Veto.Energy.Fit.Range", &lowLimit, &highLimit);
    Int_t firstBin = h->GetXaxis()->FindBin(lowLimit);
    Int_t lastBin = h->GetXaxis()->FindBin(highLimit);

    // create projection
    TH1* fitHisto;
    if (isMC)
    {
        sprintf(tmp, "Proj_MC_%d", elem);
        if (fFitHistoMC) delete fFitHistoMC;
        fFitHistoMC = (TH1D*) h->ProjectionY(tmp, firstBin, lastBin, "e");
        fitHisto = fFitHistoMC;
    }
    else
    {
        sprintf(tmp, "Proj_%d", elem);
        if (fFitHisto) delete fFitHisto;
        fFitHisto = (TH1D*) h->ProjectionY(tmp, firstBin, lastBin, "e");
        fitHisto = fFitHisto;
    }

    // estimate peak position
    TSpectrum s;
    if (isMC) s.Search(fitHisto, 5, "goff nobackground", 0.03);
    else s.Search(fitHisto, 5, "goff nobackground", 0.05);
    Double_t peak = TMath::MaxElement(s.GetNPeaks(), s.GetPositionX());

    // create fitting function
    TF1* fitfunc;
    if (isMC)
    {
        if (fFitFuncMC) delete fFitFuncMC;
        sprintf(tmp, "fFunc_MC_%s", h->GetName());
        fFitFuncMC = new TF1(tmp, "expo(0)+gaus(2)");
        fFitFuncMC->SetLineColor(2);
        fitfunc = fFitFuncMC;
    }
    else
    {
        if (fFitFunc) delete fFitFunc;
        sprintf(tmp, "fFunc_%s", h->GetName());
        fFitFunc = new TF1(tmp, "expo(0)+gaus(2)");
        fFitFunc->SetLineColor(2);
        fitfunc = fFitFunc;
    }

    // apply re-fit
    if (fIsReFit)
    {
        if (isMC)
            peak = fLineMC->GetPos();
        else
            peak = fLine->GetPos();
    }

    // prepare fitting function
    Double_t range = 30./lowLimit+0.3;
    Double_t peak_range = 0.2;
    fitfunc->SetRange(peak - range, peak + range*2);
    fitfunc->SetParameter(2, fitHisto->GetXaxis()->FindBin(peak));
    fitfunc->SetParLimits(2, 0, 100000);
    fitfunc->SetParameter(3, peak);
    fitfunc->SetParLimits(3, peak - peak_range, peak + peak_range);
    fitfunc->SetParameter(4, 1);
    fitfunc->SetParLimits(4, 0.1, 10);

    // perform first fit
    fitHisto->Fit(fitfunc, "RB0Q");

    // adjust fitting range
    Double_t sigma = fitfunc->GetParameter(4);
    fitfunc->SetRange(peak - 3*sigma, peak + range+3*sigma);

    // perform second fit
    fitHisto->Fit(fitfunc, "RB0Q");

    // format and plot stuff
    if (isMC)
    {
        fPeakMC = fitfunc->GetParameter(3);
        fLineMC->SetPos(fPeakMC);
        fCanvasFit->cd(2);
    }
    else
    {
        fPeak = fitfunc->GetParameter(3);
        fLine->SetPos(fPeak);
        fCanvasFit->cd(4);
    }

    // draw histogram
    fitHisto->Draw("hist");
    fitfunc->Draw("same");
    if (isMC)
        fLineMC->Draw();
    else
        fLine->Draw();

    fCanvasFit->Update();
}

//______________________________________________________________________________
void TCCalibVetoEnergy::Fit(Int_t elem)
{
    // Perform the fit of the element 'elem'.

    Char_t tmp[256];

    // create histogram name
    sprintf(tmp, "%s_%03d", fHistoName.Data(), elem);

    // get histogram
    if (fMainHisto) delete fMainHisto;
    fMainHisto = (TH2*) fFileManager->GetHistogram(tmp);
    if (!fMainHisto)
    {
        Error("Init", "Main histogram does not exist!\n");
        return;
    }
    fMainHisto->SetTitle(TString::Format("%s (data)", fMainHisto->GetTitle()).Data());

    // get MC histogram
    if (fMCHisto) delete fMCHisto;
    fMCHisto = (TH2*) fMCFile->Get(tmp);
    if (!fMCHisto)
    {
        Error("Init", "MC histogram does not exist!\n");
        return;
    }
    fMCHisto->SetTitle(TString::Format("%s (MC)", fMCHisto->GetTitle()).Data());

    // draw main histogram
    fCanvasFit->cd(3)->SetLogz();
    TCUtils::FormatHistogram(fMainHisto, "Veto.Energy.Histo.Fit");
    fMainHisto->Draw("colz");
    fCanvasFit->cd(1)->SetLogz();
    TCUtils::FormatHistogram(fMCHisto, "Veto.Energy.Histo.Fit");
    fMCHisto->Draw("colz");
    fCanvasFit->Update();

    // check for sufficient statistics
    if (fMainHisto->GetEntries())
    {
        // fit the energy slice
        FitSlice((TH2*)fMCHisto, elem, kTRUE);
        FitSlice((TH2*)fMainHisto, elem, kFALSE);
    }

    fCanvasResult->cd();
    fOverviewHisto->Draw("E1");
    fCanvasResult->Update();
}

//______________________________________________________________________________
void TCCalibVetoEnergy::Calculate(Int_t elem)
{
    // Calculate the new value of the element 'elem'.

    Bool_t unchanged = kFALSE;

    // check if fit was performed
    if (fMainHisto->GetEntries())
    {
        // check if line position was modified by hand
        if (fLine->GetPos() != fPeak) fPeak = fLine->GetPos();
        if (fLineMC->GetPos() != fPeakMC) fPeakMC = fLineMC->GetPos();

        // calculate the new gain
        fNewVal[elem] = fOldVal[elem] * (fPeakMC / fPeak);

        // if new value is negative take old
        if (fNewVal[elem] < 0)
        {
            fNewVal[elem] = fOldVal[elem];
            unchanged = kTRUE;
        }

        // update overview histogram
        fOverviewHisto->SetBinContent(elem+1, fPeak);
        fOverviewHisto->SetBinError(elem+1, 0.0000001);
    }
    else
    {
        // do not change old value
        fNewVal[elem] = fOldVal[elem];
        unchanged = kTRUE;
    }

    // user information
    printf("Element: %03d    peak: %.3f    peak MC: %.3f"
           "old gain: %12.8f    new gain: %12.8f    diff: %6.2f %%",
           elem, fPeak, fPeakMC, fOldVal[elem], fNewVal[elem],
           TCUtils::GetDiffPercent(fOldVal[elem], fNewVal[elem]));
    if (unchanged) printf("    -> unchanged");
    printf("\n");
}

