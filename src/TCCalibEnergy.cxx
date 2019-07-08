/*************************************************************************
 * Author: Irakli Keshelashvili, Dominik Werthmueller, Thomas Strub
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibEnergy                                                        //
//                                                                      //
// Base energy calibration module class.                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TLine.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF1.h"

#include "TCCalibEnergy.h"
#include "TCMySQLManager.h"
#include "TCFileManager.h"
#include "TCUtils.h"
#include "TCFitUtils.h"
#include "TCLine.h"

ClassImp(TCCalibEnergy)

//______________________________________________________________________________
TCCalibEnergy::TCCalibEnergy(const Char_t* name, const Char_t* title, const Char_t* data,
                             Int_t nElem)
    : TCCalib(name, title, data, nElem)
{
    // Constructor.

    // init members
    fPi0Pos = 0;
    fLine = 0;
}

//______________________________________________________________________________
TCCalibEnergy::~TCCalibEnergy()
{
    // Destructor.

    if (fLine) delete fLine;
}

//______________________________________________________________________________
void TCCalibEnergy::Init()
{
    // Init the module.

    Char_t tmp[256];

    // init members
    fPi0Pos = 0;
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

    // read old parameters (only from first set)
    TCMySQLManager::GetManager()->ReadParameters(fData, fCalibration.Data(), fSet[0], fOldVal, fNelem);

    // copy to new parameters
    for (Int_t i = 0; i < fNelem; i++) fNewVal[i] = fOldVal[i];

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
    fOverviewHisto = new TH1F("Overview", ";Element;2#gamma inv. mass [MeV]", fNelem, 0, fNelem);
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
void TCCalibEnergy::Fit(Int_t elem)
{
    // Perform the fit of the element 'elem'.

    Char_t tmp[256];

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
    fFitHisto->Draw("hist");

    // check for sufficient statistics
    if (fFitHisto->Integral() > 100 && !IsIgnored(elem))
    {
        // delete old function
        if (fFitFunc) delete fFitFunc;
        sprintf(tmp, "fEnergy_%i", elem);
        fFitFunc = new TF1(tmp, "gaus(0)+pol3(3)");
        fFitFunc->SetLineColor(2);

        // set peak position
        if (fIsReFit)
        {
            fPi0Pos = fLine->GetPos();
        }
        else
        {
            // estimate peak position
            fPi0Pos = fFitHisto->GetBinCenter(fFitHisto->GetMaximumBin());
            if (fPi0Pos < 100 || fPi0Pos > 160) fPi0Pos = 135;
        }

        // configure fitting function
        if (this->InheritsFrom("TCCalibCBEnergy"))
        {
            fFitFunc->SetRange(fPi0Pos - 10, fPi0Pos + 10);
            fFitFunc->SetParameters(fFitHisto->GetMaximum(), fPi0Pos, 11, 1, 1, 1, 0.1);
            fFitFunc->SetParLimits(1, 130, 140);
            fFitFunc->SetParLimits(2, 7, 18);
        }
        else if (this->InheritsFrom("TCCalibTAPSEnergyLG"))
        {
            fFitFunc->SetRange(60, 200);
            fFitFunc->SetParameters(fFitHisto->GetMaximum(), fPi0Pos, 10, 1, 1, 1, 0.1);
            fFitFunc->SetParLimits(0, 1, fFitHisto->GetMaximum()*1.5);
            fFitFunc->SetParLimits(1, 115, 140);
            fFitFunc->SetParLimits(2, 5, 15);
            fFitFunc->FixParameter(6, 0);
        }

        // set +/- 3% peak position limits
        if (fIsReFit) fFitFunc->SetParLimits(1, (1. - 0.03)*fPi0Pos, (1. + 0.03)*fPi0Pos);

        // fit
        TCFitUtils::ReFit(fFitHisto, fFitFunc, "RBQ0", 10);

        // final results
        fPi0Pos = fFitFunc->GetParameter(1);

        // check if mass is in normal range
        if (!fIsReFit &&
            (fPi0Pos < fFitHisto->GetXaxis()->GetXmin() || fPi0Pos > fFitHisto->GetXaxis()->GetXmax())) fPi0Pos = 135;

        // set indicator line
        fLine->SetPos(fPi0Pos);

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
void TCCalibEnergy::ReCalculateAll()
{
    // Alternative calculation of new gains. No need for a convergence factor.
    //
    // Invariant mass m for two elements i != j with gains g_i, g_j and adc
    // values A_i, A_j (i.e., energie E = g*A):
    //
    //     m   \propto \sqrt(E_i * E_j) = \sqrt(g_i*A_i * g_j*A_j)
    //  => m^2 \propto g_i*A_i * g_j*A_j
    //
    // Let m_{mean} = mean value of peak pos for j != i. Then
    //
    //     g_i(new) = g_i(old) * m0^2/m^2 * m_{mean}/m0
    //
    // Special cases:
    //
    //    1) m_{mean} == m0:    g_i(new) = g_i(old) * m0^2/m^2
    //    2) m_{mean} == m:     g_i(new) = g_i(old) * m0  /m
    //
    // Case #2 matches standart TAPS energy calibration (1 CB hit, 1 TAPS hit)
    // because the photon in CB is already calibrated.

    // get sum of pi0 positions
    Double_t sum = 0.;
    for (Int_t i = 0; i < fNelem; i++)
    {
        if (fOverviewHisto->GetBinContent(i+1) > 0.)
            sum += fOverviewHisto->GetBinContent(i+1);
        else
            sum += TCConfig::kPi0Mass;
    }

    // loop over all elements and set new gain
    for (Int_t i = 0; i < fNelem; i++)
    {
        Double_t pi0pos = fOverviewHisto->GetBinContent(i+1);

        if (pi0pos > 0.)
        {
            // get average of all except i
            Double_t mean = (sum - pi0pos) / (fNelem - 1);

            // calculate
            fNewVal[i] = fOldVal[i] * (TCConfig::kPi0Mass * TCConfig::kPi0Mass) / (pi0pos * pi0pos) * (mean / TCConfig::kPi0Mass);
        }
    }
}

//______________________________________________________________________________
void TCCalibEnergy::Calculate(Int_t elem)
{
    // Calculate the new value of the element 'elem'.

    Bool_t unchanged = kFALSE;

    // check if fit was performed
    if (fFitHisto->Integral() > 100 && !IsIgnored(elem))
    {
        // check if line position was modified by hand
        if (fLine->GetPos() != fPi0Pos) fPi0Pos = fLine->GetPos();

        // calculate the new offset
        Double_t a = TCConfig::kPi0Mass / fPi0Pos;
        //Double_t a = TCConfig::kPi0Mass*TCConfig::kPi0Mass / fPi0Pos / fPi0Pos;

        Double_t b = a - 1.;
        fNewVal[elem] = fOldVal[elem] + fOldVal[elem]*b*fConvergenceFactor;

        // if new value is negative take old
        if (fNewVal[elem] < 0)
        {
            fNewVal[elem] = fOldVal[elem];
            unchanged = kTRUE;
        }

        // update overview histogram
        fOverviewHisto->SetBinContent(elem+1, fPi0Pos);
        fOverviewHisto->SetBinError(elem+1, 0.0000001);

        // update average calculation
        fAvr += fPi0Pos;
        fAvrDiff += TMath::Abs(fPi0Pos - TCConfig::kPi0Mass);
        fNcalc++;
    }
    else
    {
        // do not change old value
        fNewVal[elem] = fOldVal[elem];
        unchanged = kTRUE;
    }

    // user information
    printf("Element: %03d    Pi0: %12.8f    "
           "old gain: %12.8f    new gain: %12.8f    diff: %6.2f %%",
           elem, fPi0Pos, fOldVal[elem], fNewVal[elem],
           TCUtils::GetDiffPercent(fOldVal[elem], fNewVal[elem]));
    if (unchanged) printf("    -> unchanged");
    if (this->InheritsFrom("TCCalibCBEnergy"))
    {
        if (TCUtils::IsCBHole(elem)) printf(" (hole)");
    }
    printf("\n");

    // show average
    if (elem == fNelem-1)
    {
        fAvr /= (Double_t)fNcalc;
        fAvrDiff /= (Double_t)fNcalc;
        printf("Average pi0 position           : %.3f MeV\n", fAvr);
        printf("Average difference to pi0 mass : %.3f MeV\n", fAvrDiff);
    }
}

//______________________________________________________________________________
void TCCalibEnergy::WriteValues()
{
    // Write the obtained calibration values to the database.

    // recalculate all gains
    if (this->InheritsFrom("TCCalibCBEnergy")) ReCalculateAll();

    // call parent's function
    TCCalib::WriteValues();
}
