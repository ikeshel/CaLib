/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibPeakFit                                                       //
//                                                                      //
// Generic peak fitting calibration module class.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"

#include "TCCalibPeakFit.h"
#include "TCFileManager.h"
#include "TCUtils.h"
#include "TCLine.h"
#include "TCReadConfig.h"

ClassImp(TCCalibPeakFit)

//______________________________________________________________________________
TCCalibPeakFit::TCCalibPeakFit(const Char_t* name, const Char_t* title, const Char_t* data,
                               Int_t nElem)
    : TCCalib(name, title, data, nElem)
{
    // Constructor.

    // init members
    fMean = 0;
    fLine = 0;
}

//______________________________________________________________________________
TCCalibPeakFit::~TCCalibPeakFit()
{
    // Destructor.

    if (fLine) delete fLine;
}

//______________________________________________________________________________
void TCCalibPeakFit::Init()
{
    // Init the module.

    Char_t tmp[256];

    // init members
    fMean = 0;
    fLine = new TCLine();

    // configure line
    fLine->SetLineColor(4);
    fLine->SetLineWidth(3);

    // get histogram name
    sprintf(tmp, "%s.Histo.Fit.Name", GetName());
    if (!TCReadConfig::GetReader()->GetConfig(tmp))
    {
        Error("Init", "Histogram name was not found in configuration!");
        return;
    }
    else fHistoName = *TCReadConfig::GetReader()->GetConfig(tmp);

    // sum up all files contained in this runset
    TCFileManager f(fData, fCalibration.Data(), fNset, fSet);

    // get the main calibration histogram
    fMainHisto = f.GetHistogram(fHistoName.Data());
    if (!fMainHisto)
    {
        Error("Init", "Main histogram does not exist!\n");
        return;
    }

    // create the overview histogram
    fOverviewHisto = new TH1F("Overview", ";Element;Peak position", fNelem, 0, fNelem);
    fOverviewHisto->SetMarkerStyle(2);
    fOverviewHisto->SetMarkerColor(4);

    // draw main histogram
    fCanvasFit->Divide(1, 2, 0.001, 0.001);
    fCanvasFit->cd(1)->SetLogz();
    sprintf(tmp, "%s.Histo.Fit", GetName());
    TCUtils::FormatHistogram(fMainHisto, tmp);
    fMainHisto->Draw("colz");

    // draw the overview histogram
    fCanvasResult->cd();
    sprintf(tmp, "%s.Histo.Overview", GetName());
    TCUtils::FormatHistogram(fOverviewHisto, tmp);
    fOverviewHisto->Draw("P");
}

//______________________________________________________________________________
void TCCalibPeakFit::Fit(Int_t elem)
{
    // Perform the fit of the element 'elem'.

    Char_t tmp[256];

    // init positon
    fMean = 0;

    // create histogram projection for this element
    sprintf(tmp, "ProjHisto_%i", elem);
    TH2* h2 = (TH2*) fMainHisto;
    if (fFitHisto) delete fFitHisto;
    fFitHisto = (TH1D*) h2->ProjectionX(tmp, elem+1, elem+1, "e");

    // draw histogram
    fFitHisto->SetFillColor(35);
    fCanvasFit->cd(2);
    sprintf(tmp, "%s.Histo.Fit", GetName());
    TCUtils::FormatHistogram(fFitHisto, tmp);
    //fFitHisto->Draw("hist");
    fFitHisto->Draw();

    // read fit config
    Double_t lowLimit = fFitHisto->GetXaxis()->GetXmin();
    Double_t highLimit = fFitHisto->GetXaxis()->GetXmax();
    sprintf(tmp, "%s.Histo.Fit.Range", GetName());
    TCReadConfig::GetReader()->GetConfigDoubleDouble(tmp, &lowLimit, &highLimit);

    // check for sufficient statistics
    if (fFitHisto->GetEntries() && !IsIgnored(elem))
    {
        // delete old function
        if (fFitFunc) delete fFitFunc;
        sprintf(tmp, "fFunc_%i", elem);

        // get important parameter positions
        Double_t fMean = fFitHisto->GetXaxis()->GetBinCenter(fFitHisto->GetMaximumBin());
        Double_t max = fFitHisto->GetBinContent(fFitHisto->GetMaximumBin());

        // check for refit
        if (fIsReFit)
        {
            fMean = fLine->GetPos();
            lowLimit = fMean - 2;
            highLimit = fMean + 2;
        }

        // the fit function
        fFitFunc = new TF1("fFitFunc", "gaus(0)+pol2(3)", lowLimit, highLimit);
        //fFitFunc = new TF1("fFitFunc", "gaus(0)", lowLimit, highLimit);
        fFitFunc->SetLineColor(2);

        // configure fitting function
        fFitFunc->SetParameters(max, fMean, 3, 1, 0.1, 0.1);
        //fFitFunc->SetParLimits(0, 0.1, max*10);

        // second iteration
        for (Int_t i = 0; i < 10; i++)
            if (!fFitHisto->Fit(fFitFunc, "RQ0")) break;

        // final results
        fMean = fFitFunc->GetParameter(1);

        // correct bad position
        if (fMean < fFitHisto->GetXaxis()->GetXmin() || fMean > fFitHisto->GetXaxis()->GetXmax())
            fMean = 0.5 * (fFitHisto->GetXaxis()->GetXmin() + fFitHisto->GetXaxis()->GetXmax());

        // draw mean indicator line
        fLine->SetPos(fMean);

        // draw fitting function
        if (fFitFunc) fFitFunc->Draw("same");

        // draw indicator line
        fLine->Draw();
    }

    // update canvas
    fCanvasFit->Update();

    // update overview
    if (elem % 20 == 0)
    {
        fCanvasResult->cd();
        fOverviewHisto->Draw("E1");
        fCanvasResult->Update();
    }
}

//______________________________________________________________________________
void TCCalibPeakFit::Calculate(Int_t elem)
{
    // Just print the peak position in this class.

    Char_t tmp[256];
    Bool_t unchanged = kFALSE;

    // check if fit was performed
    if (fFitHisto->GetEntries() && !IsIgnored(elem))
    {
        // check if line position was modified by hand
        if (fLine->GetPos() != fMean) fMean = fLine->GetPos();
        fNewVal[elem] = fMean;

        // update overview histogram
        fOverviewHisto->SetBinContent(elem + 1, fMean);
        fOverviewHisto->SetBinError(elem + 1, 0.0000001);

        // update average calculation
        fAvr += fMean;
        fAvrDiff += TMath::Abs(fMean);
        fNcalc++;
    }
    else
    {
        // do not change old value
        unchanged = kTRUE;
    }

    // user information
    printf("Element: %03d    Peak: %12.8f", elem, fMean);
    if (unchanged) printf("    -> unchanged");

    if (this->InheritsFrom("TCCalibPeakFitCB"))
    {
        if (TCUtils::IsCBHole(elem)) printf(" (hole)");
    }
    printf("\n");

    // show average
    if (elem == fNelem-1)
    {
        fAvr /= (Double_t)fNcalc;
        fAvrDiff /= (Double_t)fNcalc;
        printf("Average center          : %.3f ns\n", fAvr);
        printf("Average difference to 0 : %.3f ns\n", fAvrDiff);
    }
}

//______________________________________________________________________________
void TCCalibPeakFit::PrintValues()
{
    // Print out the values for all elements.

    // loop over elements
    printf("\n");
    for (Int_t i = 0; i < fNelem; i++)
        printf("Element: %03d    peak position: %12.8f\n", i, fNewVal[i]);
    printf("\n");
}

//______________________________________________________________________________
void TCCalibPeakFit::WriteValues()
{
    // Write the peak positions to an ASCII file.

    // open the file
    FILE* fout = fopen("peak_pos.dat", "w");

    // loop over elements
    for (Int_t i = 0; i < fNelem; i++)
        fprintf(fout, "%f\n", fNewVal[i]);

    // close the file
    fclose(fout);

    Info("WriteValues", "Peak positions were written to peak_pos.dat");
}

