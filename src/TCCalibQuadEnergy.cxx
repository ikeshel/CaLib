/*************************************************************************
 * Author: Dominik Werthmueller, Thomas Strub
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibQuadEnergy                                                    //
//                                                                      //
// Base quadratic energy correction module class.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TH2.h"
#include "TF1.h"
#include "TCLine.h"
#include "TCanvas.h"
#include "TMath.h"

#include "TCCalibQuadEnergy.h"
#include "TCMySQLManager.h"
#include "TCFileManager.h"
#include "TCUtils.h"

ClassImp(TCCalibQuadEnergy)

//______________________________________________________________________________
TCCalibQuadEnergy::TCCalibQuadEnergy(const Char_t* name, const Char_t* title, const Char_t* data,
                                     Int_t nElem)
    : TCCalib(name, title, data, nElem)
{
    // Constructor.

    // init members
    fPar0Old = 0;
    fPar1Old = 0;
    fPar0New = 0;
    fPar1New = 0;
    fMainHisto2 = 0;
    fMainHisto3 = 0;
    fFitHisto1b = 0;
    fFitHisto2 = 0;
    fFitHisto3 = 0;
    fFitFunc1b = 0;
    fPi0Pos = 0;
    fEtaPos = 0;
    fPi0MeanE = 0;
    fEtaMeanE = 0;
    fLinePi0 = 0;
    fLineEta = 0;
    fLineMeanEPi0 = 0;
    fLineMeanEEta = 0;
    fPi0PosHisto = 0;
    fEtaPosHisto = 0;
    fPi0MeanEHisto = 0;
    fEtaMeanEHisto = 0;
}

//______________________________________________________________________________
TCCalibQuadEnergy::~TCCalibQuadEnergy()
{
    // Destructor.

    if (fPar0Old) delete [] fPar0Old;
    if (fPar1Old) delete [] fPar1Old;
    if (fPar0New) delete [] fPar0New;
    if (fPar1New) delete [] fPar1New;
    if (fMainHisto2) delete fMainHisto2;
    if (fMainHisto3) delete fMainHisto3;
    if (fFitHisto1b) delete fFitHisto1b;
    if (fFitHisto2) delete fFitHisto2;
    if (fFitHisto3) delete fFitHisto3;
    if (fFitFunc1b) delete fFitFunc1b;
    if (fLinePi0) delete fLinePi0;
    if (fLineEta) delete fLineEta;
    if (fLineMeanEPi0) delete fLineMeanEPi0;
    if (fLineMeanEEta) delete fLineMeanEEta;
    if (fPi0PosHisto) delete fPi0PosHisto;
    if (fEtaPosHisto) delete fEtaPosHisto;
    if (fPi0MeanEHisto) delete fPi0MeanEHisto;
    if (fEtaMeanEHisto) delete fEtaMeanEHisto;
}

//______________________________________________________________________________
void TCCalibQuadEnergy::Init()
{
    // Init the module.

    Char_t tmp[256];

    // init members
    fPar0Old = new Double_t[fNelem];
    fPar1Old = new Double_t[fNelem];
    fPar0New = new Double_t[fNelem];
    fPar1New = new Double_t[fNelem];
    fFitHisto1b = 0;
    fFitHisto2 = 0;
    fFitHisto3 = 0;
    fFitFunc1b = 0;
    fPi0Pos = 0;
    fEtaPos = 0;
    fPi0MeanE = 0;
    fEtaMeanE = 0;
    fLinePi0 = new TCLine();
    fLineEta = new TCLine();
    fLineMeanEPi0 = new TCLine();
    fLineMeanEEta = new TCLine();

    // configure lines
    fLinePi0->SetLineColor(4);
    fLinePi0->SetLineWidth(3);
    fLineEta->SetLineColor(4);
    fLineEta->SetLineWidth(3);
    fLineMeanEPi0->SetLineColor(4);
    fLineMeanEPi0->SetLineWidth(3);
    fLineMeanEEta->SetLineColor(4);
    fLineMeanEEta->SetLineWidth(3);

    // get main histogram name
    sprintf(tmp, "%s.Histo.Fit.Name", GetName());
    if (!TCReadConfig::GetReader()->GetConfig(tmp))
    {
        Error("Init", "Histogram name was not found in configuration!");
        return;
    }
    else fHistoName = *TCReadConfig::GetReader()->GetConfig(tmp);

    // get mean pi0 energy histogram name
    TString hMeanPi0Name;
    sprintf(tmp, "%s.Histo.MeanE.Pi0.Name", GetName());
    if (!TCReadConfig::GetReader()->GetConfig(tmp))
    {
        Error("Init", "Histogram name was not found in configuration!");
        return;
    }
    else hMeanPi0Name = *TCReadConfig::GetReader()->GetConfig(tmp);

    // get mean eta energy histogram name
    TString hMeanEtaName;
    sprintf(tmp, "%s.Histo.MeanE.Eta.Name", GetName());
    if (!TCReadConfig::GetReader()->GetConfig(tmp))
    {
        Error("Init", "Histogram name was not found in configuration!");
        return;
    }
    else hMeanEtaName = *TCReadConfig::GetReader()->GetConfig(tmp);

    // read old parameters (only from first set)
    if (this->InheritsFrom("TCCalibCBQuadEnergy"))
    {
        TCMySQLManager::GetManager()->ReadParameters("Data.CB.Energy.Quad.Par0", fCalibration.Data(), fSet[0], fPar0Old, fNelem);
        TCMySQLManager::GetManager()->ReadParameters("Data.CB.Energy.Quad.Par1", fCalibration.Data(), fSet[0], fPar1Old, fNelem);
    }
    else if (this->InheritsFrom("TCCalibTAPSQuadEnergy"))
    {
        TCMySQLManager::GetManager()->ReadParameters("Data.TAPS.Energy.Quad.Par0", fCalibration.Data(), fSet[0], fPar0Old, fNelem);
        TCMySQLManager::GetManager()->ReadParameters("Data.TAPS.Energy.Quad.Par1", fCalibration.Data(), fSet[0], fPar1Old, fNelem);
    }

    // copy parameters
    for (Int_t i = 0; i < fNelem; i++)
    {
        fPar0New[i] = fPar0Old[i];
        fPar1New[i] = fPar1Old[i];
    }

    // sum up all files contained in this runset
    TCFileManager f(fData, fCalibration.Data(), fNset, fSet);

    // get the main calibration histogram
    fMainHisto = f.GetHistogram(fHistoName.Data());
    if (!fMainHisto)
    {
        Error("Init", "Main histogram does not exist!\n");
        return;
    }

    // get the pi0 mean energy histogram
    fMainHisto2 = (TH2*) f.GetHistogram(hMeanPi0Name.Data());
    if (!fMainHisto2)
    {
        Error("Init", "Pi0 mean energy histogram does not exist!\n");
        return;
    }

    // get the eta mean energy histogram
    fMainHisto3 = (TH2*) f.GetHistogram(hMeanEtaName.Data());
    if (!fMainHisto3)
    {
        Error("Init", "Eta mean energy histogram does not exist!\n");
        return;
    }

    // create the pi0 overview histogram
    fPi0PosHisto = new TH1F("Pi0 position overview", ";Element;#pi^{0} peak position [MeV]", fNelem, 0, fNelem);
    fPi0PosHisto->SetMarkerStyle(2);
    fPi0PosHisto->SetMarkerColor(4);
    sprintf(tmp, "%s.Histo.Overview.Pi0", GetName());
    TCUtils::FormatHistogram(fPi0PosHisto, tmp);

    // create the eta overview histogram
    fEtaPosHisto = new TH1F("Eta position overview", ";Element;#eta peak position [MeV]", fNelem, 0, fNelem);
    fEtaPosHisto->SetMarkerStyle(2);
    fEtaPosHisto->SetMarkerColor(4);
    sprintf(tmp, "%s.Histo.Overview.Eta", GetName());
    TCUtils::FormatHistogram(fEtaPosHisto, tmp);

    fPi0MeanEHisto = new TH1F("Pi0MeanE", ";Element;Mean photon energy [MeV]", fNelem, 0, fNelem);
    fEtaMeanEHisto = new TH1F("EtaMeanE", ";Element;Mean photon energy [MeV]", fNelem, 0, fNelem);

    // prepare fit histogram canvas
    fCanvasFit->Divide(1, 4, 0.001, 0.001);

    // draw the overview histograms
    fCanvasResult->Divide(1, 2, 0.001, 0.001);
    fCanvasResult->cd(1);
    fPi0PosHisto->Draw("P");
    fCanvasResult->cd(2);
    fEtaPosHisto->Draw("P");
}

//______________________________________________________________________________
void TCCalibQuadEnergy::Fit(Int_t elem)
{
    // Perform the fit of the element 'elem'.

    Char_t tmp[256];

    // get the 2g invariant mass histograms
    sprintf(tmp, "ProjHisto_%d", elem);
    TH2* h2 = (TH2*) fMainHisto;
    if (fFitHisto) delete fFitHisto;
    if (fFitHisto1b) delete fFitHisto1b;
    fFitHisto = (TH1*) h2->ProjectionX(tmp, elem+1, elem+1, "e");
    sprintf(tmp, "ProjHisto_%db", elem);
    fFitHisto1b = (TH1*) fFitHisto->Clone(tmp);
    sprintf(tmp, "%s.Histo.Fit.Pi0.IM", GetName());
    TCUtils::FormatHistogram(fFitHisto, tmp);
    sprintf(tmp, "%s.Histo.Fit.Eta.IM", GetName());
    TCUtils::FormatHistogram(fFitHisto1b, tmp);

    // get pi0 mean energy projection
    sprintf(tmp, "ProjHistoMeanPi0_%d", elem);
    h2 = (TH2*) fMainHisto2;
    if (fFitHisto2) delete fFitHisto2;
    fFitHisto2 = h2->ProjectionX(tmp, elem+1, elem+1, "e");
    sprintf(tmp, "%s.Histo.Fit.Pi0.MeanE", GetName());
    TCUtils::FormatHistogram(fFitHisto2, tmp);

    // get eta mean energy projection
    sprintf(tmp, "ProjHistoMeanEta_%d", elem);
    h2 = (TH2*) fMainHisto3;
    if (fFitHisto3) delete fFitHisto3;
    fFitHisto3 = h2->ProjectionX(tmp, elem+1, elem+1, "e");
    sprintf(tmp, "%s.Histo.Fit.Eta.MeanE", GetName());
    TCUtils::FormatHistogram(fFitHisto3, tmp);

    // draw pi0
    fCanvasFit->cd(1);
    fFitHisto->SetFillColor(35);
    fFitHisto->Draw("hist");

    // draw eta
    fCanvasFit->cd(2);
    fFitHisto1b->SetFillColor(35);
    fFitHisto1b->Draw("hist");

    // draw pi0 mean energy
    fCanvasFit->cd(3);
    fFitHisto2->SetFillColor(35);
    fFitHisto2->Draw("hist");

    // draw eta mean energy
    fCanvasFit->cd(4);
    fFitHisto3->SetFillColor(35);
    fFitHisto3->Draw("hist");

    // check for sufficient statistics
    if (fFitHisto->GetEntries() && !IsIgnored(elem))
    {
        // delete old functions
        if (fFitFunc) delete fFitFunc;
        if (fFitFunc1b) delete fFitFunc1b;

        // create pi0 fitting function
        sprintf(tmp, "fPi0_%i", elem);
        fFitFunc = new TF1(tmp, "gaus(0)+pol2(3)", 100, 170);
        fFitFunc->SetLineColor(2);

        // create eta fitting function
        sprintf(tmp, "fEta_%i", elem);
        fFitFunc1b = new TF1(tmp, "gaus(0)+pol1(3)", 350, 650);
        fFitFunc1b->SetLineColor(2);

        // get x-axis range
        Double_t xmin = fFitHisto1b->GetXaxis()->GetBinCenter(fFitHisto1b->GetXaxis()->GetFirst());
        Double_t xmax = fFitHisto1b->GetXaxis()->GetBinCenter(fFitHisto1b->GetXaxis()->GetLast());

        // set new range & get the peak position of eta
        fFitHisto1b->GetXaxis()->SetRangeUser(500, 600);
        Double_t fMaxEta = fFitHisto1b->GetBinCenter(fFitHisto1b->GetMaximumBin());
        fFitHisto1b->GetXaxis()->SetRangeUser(xmin, xmax);

        // init peak positions
        fPi0Pos = 135.;
        fEtaPos = fMaxEta;

        if (fIsReFit)
        {
            fPi0Pos = fLinePi0->GetPos();
            fEtaPos = fLineEta->GetPos();
        }

        // configure fitting functions
        // pi0
        fFitFunc->SetParameters(fFitHisto->GetMaximum(), fPi0Pos, 10, 1, 1, 1);
        fFitFunc->SetParLimits(0, 0.1*fFitHisto->GetMaximum(), 1.5*fFitHisto->GetMaximum());
        fFitFunc->SetParLimits(1, fPi0Pos - 15., fPi0Pos+15.);
        fFitFunc->SetParLimits(2, 2, 40);

        // eta
        fFitFunc1b->SetParameters(fFitHisto1b->GetMaximum(), fEtaPos, 15, 1, 1, 1, 0.1);
        fFitFunc1b->SetParLimits(0, 0.1*fFitHisto1b->GetMaximum(), 1.5*fFitHisto1b->GetMaximum());
        fFitFunc1b->SetParLimits(1, fEtaPos - 30, fEtaPos + 30);
        fFitFunc1b->SetParLimits(2, 10, 50);
        //fFitFunc1b->SetParLimits(3, 0, 100);
        //fFitFunc1b->SetParLimits(4, -1, 0);
        //fFitFunc1b->SetParLimits(5, -1, 0);//0, 50

        // set strict 3% limits for refitting
        if (fIsReFit)
        {
            fFitFunc->SetParLimits(1, (1 - 0.03)*fPi0Pos, (1 + 0.03)*fPi0Pos);
            fFitFunc1b->SetParLimits(1, (1 - 0.03)*fEtaPos, (1 + 0.03)*fEtaPos);
        }

        // fit peaks
        for (Int_t i = 0; i < 10; i++)
            if (!fFitHisto->Fit(fFitFunc, "RBQ0")) break;
        for (Int_t i = 0; i < 10; i++)
            if (!fFitHisto1b->Fit(fFitFunc1b, "RBQ0W")) break;

        // get results
        fPi0Pos = fFitFunc->GetParameter(1);
        fEtaPos = fFitFunc1b->GetParameter(1);
        fPi0MeanE = fFitHisto2->GetMean();
        fEtaMeanE = fFitHisto3->GetMean();

        // check if mass is in normal range
        if (!fIsReFit)
        {
            if (fPi0Pos < 80 || fPi0Pos > 200) fPi0Pos = 135;
            if (fEtaPos < 450 || fEtaPos > 650) fEtaPos = 547;
        }

        // set indicator lines
        fLinePi0->SetPos(fPi0Pos);
        fLineEta->SetPos(fEtaPos);

        // set lines
        fLineMeanEPi0->SetPos(fPi0MeanE);
        fLineMeanEEta->SetPos(fEtaMeanE);

        // draw pi0
        fCanvasFit->cd(1);
        if (fFitFunc) fFitFunc->Draw("same");
        fLinePi0->Draw();

        // draw eta
        fCanvasFit->cd(2);
        if (fFitFunc1b) fFitFunc1b->Draw("same");
        fLineEta->Draw();

        // draw pi0 mean energy
        fCanvasFit->cd(3);
        fLineMeanEPi0->Draw();

        // draw eta mean energy
        fCanvasFit->cd(4);
        fLineMeanEEta->Draw();
    }

    // update canvas
    fCanvasFit->Update();

    // update overview
    if (elem % 20 == 0)
    {
        fCanvasResult->cd(1);
        fPi0PosHisto->Draw("E1");
        fCanvasResult->cd(2);
        fEtaPosHisto->Draw("E1");
        fCanvasResult->Update();
    }
}

//______________________________________________________________________________
void TCCalibQuadEnergy::CalculateNewPar(Double_t& par0, Double_t& par1,
                                        Double_t pi0Pos, Double_t pi0_mean_e,
                                        Double_t etaPos, Double_t eta_mean_e,
                                        Double_t pi0_factor /*= 1.*/,
                                        Double_t eta_factor /*= 1.*/,
                                        Double_t conv_factor /*= 1.*/)
{
    // Calculates new parameters from old ones and returns them via the first
    // two arguments.

    // Denote
    //   E0          energy w/o quadratic correction
    //   m0          theoretical mass of meson
    //   m           mass of meson with old quad. corr
    //   E , a , b   energy and parameters of old quad. corr
    //   E', a', b'  energy and parameters of new quad. corr
    //
    //   E  = a  E0 + b  E0^2
    //   E' = a' E0 + b' E0^2
    //
    //   if b == 0  =>  E0 = E/a
    //   if b != 0  =>  E0 = (sqrt(a^2 + 4 b E) - a) / (2 b)
    //
    // Then,
    //   a  E0 + b  E0^2  propto  m^2
    //   a' E0 + b' E0^2  propto  m0^2
    //
    // gives together
    //   a' E0 + b' E0^2 = (a E0 + b E0^2) (m0/m)^2
    //   a'    + b' E0   = E/E0 * (m0/m)^2
    //
    // Define
    //   R :=  E/E0 * (m0/m)^2
    //
    //   a'    + b' E0  = R
    //
    // For two equations, e.g. eta & pion,
    // solve:
    //
    //   a' +  b' E0_1 = R_1      (1)
    //   a' +  b' E0_2 = R_2      (2)
    //
    // Solve (2) for b gives
    //  b' = (R_2 - a') / E0_2   (2')
    //
    // Insert (2') into (1)
    //  a' = [ R1 - (E0_1/E0_2)*R2 ] / [ 1 - E0_1/E0_2 ]


    // get mean energy of photons (E from above)
    Double_t e0pi0 = pi0_mean_e;
    Double_t e0eta = eta_mean_e;

    // undo quadratic energy correction (E0 from above)
    if (par1 == 0.)
    {
        e0pi0 = pi0_mean_e / par0;
        e0eta = eta_mean_e / par0;
    }
    else
    {
        e0pi0 = (TMath::Sqrt(par0*par0 + 4.*par1*pi0_mean_e) - par0) / (2.* par1);
        e0eta = (TMath::Sqrt(par0*par0 + 4.*par1*eta_mean_e) - par0) / (2.* par1);
    }

    // energy ratio (E0_1/E0_2 from above):
    Double_t mean_E_ratio = e0eta / e0pi0;

    // calculate ratios (R1 and R2 from above)
    Double_t pion_im_ratio = (TCConfig::kPi0Mass * TCConfig::kPi0Mass) / (pi0Pos * pi0Pos) * (pi0_mean_e / e0pi0) * pi0_factor;
    Double_t eta_im_ratio =  (TCConfig::kEtaMass * TCConfig::kEtaMass) / (etaPos * etaPos) * (eta_mean_e / e0eta) * eta_factor;

    // set new parameters (a', b' from above)
    Double_t par0new = (eta_im_ratio - mean_E_ratio*pion_im_ratio) / (1. - mean_E_ratio);
    Double_t par1new = (pion_im_ratio - par0new) / e0pi0;

    // check for convergence factor and non-zero old parameters
    if (conv_factor != 1.0 && par0 != 0 && par1 != 0)
    {
        // calculate product correction factor
        Double_t f0 = par0new / par0;
        Double_t f1 = par1new / par1;

        // calculate update correction factor
        Double_t b0 = f0 - 1;
        Double_t b1 = f1 - 1;

        // calculate new parameters including convergence factor
        par0new = par0 + par0*b0*conv_factor;
        par1new = par1 + par1*b1*conv_factor;
    }

    // return variables
    par0 = par0new;
    par1 = par1new;

    return;
}

//______________________________________________________________________________
void TCCalibQuadEnergy::ReCalculateAll()
{
    // Re-calculates all elements taking in account the peak position mean
    // values of all elements. Results in good convergence of the calibration.

    // init mean positions
    Double_t pi0_aver = 0;
    Double_t eta_aver = 0;

    // loop over all elements
    for (Int_t i = 0; i < fNelem; i++)
    {
        // calc mean values
        if (fPi0PosHisto->GetBinContent(i+1) != 0.)
            pi0_aver += fPi0PosHisto->GetBinContent(i+1);
        else
            pi0_aver += TCConfig::kPi0Mass;

        if (fEtaPosHisto->GetBinContent(i+1) != 0.)
            eta_aver += fEtaPosHisto->GetBinContent(i+1);
        else
            eta_aver += TCConfig::kEtaMass;
    }

    // calculate mean value
    pi0_aver /= fNelem;
    eta_aver /= fNelem;

    // loop over all elements
    for (Int_t i = 0; i < fNelem; i++)
    {
        // check whether element was calibrated
        if (fPi0PosHisto->GetBinContent(i+1) == 0.) continue;

        // get (old) parameter
        Double_t par0 = fPar0Old[i];
        Double_t par1 = fPar1Old[i];

        // get fit values
        Double_t pi0pos  = fPi0PosHisto->GetBinContent(i+1);
        Double_t etapos  = fEtaPosHisto->GetBinContent(i+1);
        Double_t pi0mean_e  = fPi0MeanEHisto->GetBinContent(i+1);
        Double_t etamean_e  = fEtaMeanEHisto->GetBinContent(i+1);

        // calculate convergence factors
        Double_t pi0_factor = (pi0_aver * fNelem - fPi0PosHisto->GetBinContent(i+1)) / (fNelem - 1) / TCConfig::kPi0Mass;
        Double_t eta_factor = (eta_aver * fNelem - fEtaPosHisto->GetBinContent(i+1)) / (fNelem - 1) / TCConfig::kEtaMass;

        // re-calculate element
        CalculateNewPar(par0, par1, pi0pos, pi0mean_e, etapos, etamean_e, pi0_factor, eta_factor, fConvergenceFactor);

        // set result
        fPar0New[i] = par0;
        fPar1New[i] = par1;
    }
}

//______________________________________________________________________________
void TCCalibQuadEnergy::Calculate(Int_t elem)
{
    // Calculate the new value of the element 'elem'.

    Bool_t no_corr = kFALSE;

    // check if fit was performed
    if (fFitHisto->GetEntries() && !IsIgnored(elem))
    {
        // check if pi0 line position was modified by hand
        if (fLinePi0->GetPos() != fPi0Pos) fPi0Pos = fLinePi0->GetPos();
        if (fLineMeanEPi0->GetPos() != fPi0MeanE) fPi0MeanE = fLineMeanEPi0->GetPos();

        // check if etaline position was modified by hand
        if (fLineEta->GetPos() != fEtaPos) fEtaPos = fLineEta->GetPos();
        if (fLineMeanEEta->GetPos() != fEtaMeanE) fEtaMeanE = fLineMeanEEta->GetPos();

        /* obsolete
        // calculate quadratic correction factors
        Double_t mean_E_ratio = fEtaMeanE / fPi0MeanE;
        Double_t pion_im_ratio = TCConfig::kPi0Mass / fPi0Pos;
        Double_t eta_im_ratio = TCConfig::kEtaMass / fEtaPos;
        fPar0[elem] = (eta_im_ratio - mean_E_ratio*pion_im_ratio) / (1. - mean_E_ratio);
        fPar1[elem] = (pion_im_ratio - fPar0[elem]) / fPi0MeanE;
        */

        // correct for wrong initial db value
        if (fPar0Old[elem] == 0.)
        {
            fPar0Old[elem] = 1.;
            fPar0New[elem] = 1.;
        }

        // get old parameters
        Double_t par0 = fPar0Old[elem];
        Double_t par1 = fPar1Old[elem];

        // calculate new values (returned via first two args) --> c.f. also ReCalculateAll()
        if (this->InheritsFrom("TCCalibCBQuadEnergy"))
        {
            // standard calculation for an IM mass from two CB photons (first iteration)
            CalculateNewPar(par0, par1, fPi0Pos, fPi0MeanE, fEtaPos, fEtaMeanE, fPi0Pos/TCConfig::kPi0Mass, fEtaPos/TCConfig::kEtaMass, fConvergenceFactor);
        }
        if (this->InheritsFrom("TCCalibTAPSQuadEnergy"))
        {
            // standard calculation for an IM mass from one CB photon and one TAPS photon (CB already calibrated)
            CalculateNewPar(par0, par1, fPi0Pos, fPi0MeanE, fEtaPos, fEtaMeanE, 1, 1, fConvergenceFactor);
        }

        // set result
        fPar0New[elem] = par0;
        fPar1New[elem] = par1;

        // check values
        if (TMath::IsNaN(fPar0New[elem]) || TMath::IsNaN(fPar1New[elem]))
        {
            Warning("Calculate", "One newly calculated parameter is NaN!");
            fPar0New[elem] = fPar0Old[elem];
            fPar1New[elem] = fPar1Old[elem];
            no_corr = kTRUE;
        }
        else
        {
            // update histograms
            fPi0PosHisto->SetBinContent(elem+1, fPi0Pos);
            fPi0PosHisto->SetBinError(elem+1, 0.0000001);
            fEtaPosHisto->SetBinContent(elem+1, fEtaPos);
            fEtaPosHisto->SetBinError(elem+1, 0.0000001);

            fPi0MeanEHisto->SetBinContent(elem+1, fPi0MeanE);
            fEtaMeanEHisto->SetBinContent(elem+1, fEtaMeanE);
        }
    }
    else //needed???
    {
        fPar0New[elem] = fPar0Old[elem];
        fPar1New[elem] = fPar1Old[elem];
        no_corr = kTRUE;
    }

    // user information
    printf("Element: %03d    Pi0 Pos.: %6.2f    Pi0 ME: %6.2f    "
           "Eta Pos.: %6.2f    Eta ME: %6.2f    Par0 old: %7.6f  new: %7.6f    Par1 old: %+4.3e  new: %+4.3e",
           elem, fPi0Pos, fPi0MeanE, fEtaPos, fEtaMeanE, fPar0Old[elem], fPar0New[elem], fPar1Old[elem], fPar1New[elem]);
    if (no_corr) printf("    -> no correction");
    if (this->InheritsFrom("TCCalibCBQuadEnergy"))
    {
        if (TCUtils::IsCBHole(elem)) printf(" (hole)");
    }
    printf("\n");
}

//______________________________________________________________________________
void TCCalibQuadEnergy::PrintValues()
{
    // Print out the old and new values for all elements.

    // loop over elements
    for (Int_t i = 0; i < fNelem; i++)
    {
        printf("Element: %03d    Par0Old: %7.6f    Par0New: %7.6f    "
                                "Par1Old: %+4.3e    Par1New: %+4.3e\n",
               i, fPar0Old[i], fPar0New[i], fPar1Old[i], fPar1New[i]);
    }
}

//______________________________________________________________________________
void TCCalibQuadEnergy::WriteValues()
{
    // Write the obtained calibration values to the database.

    // re-calculate gains for all elements
    if (this->InheritsFrom("TCCalibCBQuadEnergy")) ReCalculateAll();

    // write values to database
    for (Int_t i = 0; i < fNset; i++)
    {
        if (this->InheritsFrom("TCCalibCBQuadEnergy"))
        {
            TCMySQLManager::GetManager()->WriteParameters("Data.CB.Energy.Quad.Par0", fCalibration.Data(), fSet[i], fPar0New, fNelem);
            TCMySQLManager::GetManager()->WriteParameters("Data.CB.Energy.Quad.Par1", fCalibration.Data(), fSet[i], fPar1New, fNelem);
        }
        else if (this->InheritsFrom("TCCalibTAPSQuadEnergy"))
        {
            TCMySQLManager::GetManager()->WriteParameters("Data.TAPS.Energy.Quad.Par0", fCalibration.Data(), fSet[i], fPar0New, fNelem);
            TCMySQLManager::GetManager()->WriteParameters("Data.TAPS.Energy.Quad.Par1", fCalibration.Data(), fSet[i], fPar1New, fNelem);
        }
    }

    // save overview canvas
    SaveCanvas(fCanvasResult, "Overview");
}

