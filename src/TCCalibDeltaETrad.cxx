/*************************************************************************
 * Author: Dominik Werthmueller, Thomas Strub
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibDeltaETrad                                                    //
//                                                                      //
// Calibration module for DeltaE calibrations (traditional method).     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TMath.h"

#include "TCCalibDeltaETrad.h"
#include "TCLine.h"
#include "TCFileManager.h"
#include "TCUtils.h"
#include "TCReadConfig.h"
#include "TCMySQLManager.h"

ClassImp(TCCalibDeltaETrad)

// init static members
Double_t TCCalibDeltaETrad::fgGNP = 10;

//______________________________________________________________________________
TCCalibDeltaETrad::TCCalibDeltaETrad(const Char_t* name, const Char_t* title, const Char_t* data,
                                     Int_t nElem)
    : TCCalib(name, title, data, nElem)
{
    // Constructor.

    // init members
    fFileManager = 0;
    fPedOld = 0;
    fPedNew = 0;
    fGainOld = 0;
    fGainNew = 0;
    fPionMC = 0;
    fProtonMC = 0;
    fPionData = 0;
    fProtonData = 0;
    fPionPos = 0;
    fProtonPos = 0;
    fLinePion = 0;
    fLineProt = 0;
    fLinePionMC = 0;
    fLineProtMC = 0;
    fMCFile = 0;
    fMainHistoMC = 0;
    fFitHistoMC = 0;
    fFitFuncMC = 0;
}

//______________________________________________________________________________
TCCalibDeltaETrad::~TCCalibDeltaETrad()
{
    // Destructor.

    if (fFileManager) delete fFileManager;
    if (fPedOld) delete [] fPedOld;
    if (fPedNew) delete [] fPedNew;
    if (fGainOld) delete [] fGainOld;
    if (fGainNew) delete [] fGainNew;
    if (fPionPos) delete fPionPos;
    if (fProtonPos) delete fProtonPos;
    if (fLinePion) delete fLinePion;
    if (fLineProt) delete fLineProt;
    if (fLinePionMC) delete fLinePionMC;
    if (fLineProtMC) delete fLineProtMC;
    if (fMCFile) delete fMCFile;
    if (fMainHistoMC) delete fMainHistoMC;
    if (fFitHistoMC) delete fFitHistoMC;
    if (fFitFuncMC) delete fFitFuncMC;
}

//______________________________________________________________________________
void TCCalibDeltaETrad::Init()
{
    // Init the module.

    Char_t tmp[256];

    // init members
    fFileManager = new TCFileManager(fData, fCalibration.Data(), fNset, fSet);
    fPedOld = new Double_t[fNelem];
    fPedNew = new Double_t[fNelem];
    fGainOld = new Double_t[fNelem];
    fGainNew = new Double_t[fNelem];
    fPionMC = 0;
    fProtonMC = 0;
    fPionData = 0;
    fProtonData = 0;
    fLinePion = new TCLine();
    fLineProt = new TCLine();
    fLinePionMC = new TCLine();
    fLineProtMC = new TCLine();
    fMCFile = 0;
    fMainHistoMC = 0;
    fFitHistoMC = 0;
    fFitFuncMC = 0;

    // configure line
    fLinePion->SetLineColor(kBlue);
    fLinePion->SetLineWidth(3);
    fLineProt->SetLineColor(kGreen);
    fLineProt->SetLineWidth(3);
    fLinePionMC->SetLineColor(kBlue);
    fLinePionMC->SetLineWidth(3);
    fLinePionMC->SetLineStyle(2);
    fLineProtMC->SetLineColor(kGreen);
    fLineProtMC->SetLineWidth(3);
    fLineProtMC->SetLineStyle(2);

    // get histogram name
    sprintf(tmp, "%s.Histo.Fit.Name", GetName());
    if (!TCReadConfig::GetReader()->GetConfig(tmp))
    {
        Error("Init", "Histogram name was not found in configuration!");
        return;
    }
    else fHistoName = *TCReadConfig::GetReader()->GetConfig(tmp);

    // get MC histogram file
    TString fileMC;
    sprintf(tmp, "%s.MC.File", GetName());
    if (!TCReadConfig::GetReader()->GetConfig(tmp))
    {
        Error("Init", "MC file name was not found in configuration!");
        return;
    }
    else fileMC = *TCReadConfig::GetReader()->GetConfig(tmp);

    // read old parameters (only from first set)
    TString peds = fData;
    peds.ReplaceAll("E1", "E0");
    TCMySQLManager::GetManager()->ReadParameters(peds.Data(), fCalibration.Data(), fSet[0], fPedOld, fNelem);
    TCMySQLManager::GetManager()->ReadParameters(fData.Data(), fCalibration.Data(), fSet[0], fGainOld, fNelem);

    // copy parameters
    memcpy(fPedNew, fPedOld, fNelem * sizeof(Double_t));
    memcpy(fGainNew, fGainOld, fNelem * sizeof(Double_t));

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

    // create the pion position overview histogram
    fPionPos = new TH1F("Pion position overview", ";Element;pion peak position difference [%]", fNelem, 0, fNelem);
    fPionPos->SetMarkerStyle(2);
    fPionPos->SetMarkerColor(4);

    // create the proton position overview histogram
    fProtonPos = new TH1F("Proton position overview", ";Element;proton peak position difference [%]", fNelem, 0, fNelem);
    fProtonPos->SetMarkerStyle(2);
    fProtonPos->SetMarkerColor(4);

    // draw the overview histograms
    fCanvasResult->Divide(1, 2, 0.001, 0.001);
    fCanvasResult->cd(1);
    fPionPos->Draw("P");
    fCanvasResult->cd(2);
    fProtonPos->Draw("P");
}

//______________________________________________________________________________
Double_t TCCalibDeltaETrad::FitFunc(Double_t* x, Double_t* par)
{
   // Source: $ROOTSYS/tutorials/fit/langaus.C
   //
   //Fit parameters:
   //par[0]=exp constant
   //par[1]=exp slope
   //par[2]=Width (scale) parameter of Landau density
   //par[3]=Most Probable (MP, location) parameter of Landau density
   //par[4]=Total area (integral -inf to inf, normalization constant)
   //par[5]=Width (sigma) of convoluted Gaussian function
   //par[6]=gauss constant
   //par[7]=gauss mean
   //par[8]=gauss sigma
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t sc =   4.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[3] - mpshift * par[2];

      // Range of convolution integral
      xlow = x[0] - sc * par[5];
      xupp = x[0] + sc * par[5];

      step = (xupp-xlow) / fgGNP;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=fgGNP/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[2]);
         sum += fland * TMath::Gaus(x[0],xx,par[5]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[2]);
         sum += fland * TMath::Gaus(x[0],xx,par[5]);
      }

    Double_t out = (par[4] * step * sum * invsq2pi / par[3]);
    out += TMath::Exp(par[0] + par[1]*x[0]);
    out += par[6] * TMath::Gaus(x[0], par[7], par[8]);
    return out;
}

//______________________________________________________________________________
void TCCalibDeltaETrad::FitSlice(TH2* h, Bool_t isMC)
{
    // Fit the energy slice of the dE vs E histogram 'h'.

    Char_t tmp[256];

    // only refit if line was moved
    if (fIsReFit)
    {
        if (isMC)
        {
            if (fLinePionMC->GetPos() == fPionMC && fLineProtMC->GetPos() == fProtonMC)
                return;
        }
        else
        {
            if (fLinePion->GetPos() == fPionData && fLineProt->GetPos() == fProtonData)
                return;
         }
    }

    // get configuration
    Double_t lowLimit, highLimit;
    sprintf(tmp, "%s.Fit.Range.ClEnergy", GetName());
    TCReadConfig::GetReader()->GetConfigDoubleDouble(tmp, &lowLimit, &highLimit);
    Int_t firstBin = h->GetXaxis()->FindBin(lowLimit);
    Int_t lastBin = h->GetXaxis()->FindBin(highLimit);

    // create projection
    TH1* fitHisto;
    if (isMC)
    {
        sprintf(tmp, "%s_Proj_MC", h->GetName());
        if (fFitHistoMC) delete fFitHistoMC;
        fFitHistoMC = (TH1D*) h->ProjectionY(tmp, firstBin, lastBin, "e");
        fitHisto = fFitHistoMC;
    }
    else
    {
        sprintf(tmp, "%s_Proj", h->GetName());
        if (fFitHisto) delete fFitHisto;
        fFitHisto = (TH1D*) h->ProjectionY(tmp, firstBin, lastBin, "e");
        fitHisto = fFitHisto;
    }

    // rebin histogram
    sprintf(tmp, "%s.Histo.Fit.Rebin", GetName());
    Int_t rebin = TCReadConfig::GetReader()->GetConfigInt(tmp);
    if (rebin)
        fitHisto->Rebin(rebin);

    // look for peaks
    TSpectrum s;
    s.Search(fitHisto, 10, "goff nobackground", 0.05);
    Double_t pionPos = TMath::MinElement(s.GetNPeaks(), s.GetPositionX());
    Double_t protonPos = TMath::MaxElement(s.GetNPeaks(), s.GetPositionX());

    // create fitting function
    TF1* fitfunc;
    if (isMC)
    {
        if (fFitFuncMC) delete fFitFuncMC;
        sprintf(tmp, "fFunc_MC_%s", h->GetName());
        fFitFuncMC = new TF1(tmp, FitFunc, 0.4, protonPos+5, 9);
        fFitFuncMC->SetLineColor(2);
        fitfunc = fFitFuncMC;
    }
    else
    {
        if (fFitFunc) delete fFitFunc;
        sprintf(tmp, "fFunc_%s", h->GetName());
        fFitFunc = new TF1(tmp, FitFunc, 0.2*pionPos, protonPos+5, 9);
        fFitFunc->SetLineColor(2);
        fitfunc = fFitFunc;
    }

    // apply re-fit
    if (fIsReFit)
    {
        if (isMC)
        {
            pionPos = fLinePionMC->GetPos();
            protonPos = fLineProtMC->GetPos();
        }
        else
        {
            pionPos = fLinePion->GetPos();
            protonPos = fLineProt->GetPos();
        }
        fgGNP = 100;
    }
    else
    {
        fgGNP = 10;
    }

    // prepare fitting function
    fitfunc->SetParameters(9, -3.76050e-01,
                           10, pionPos, 100, 0.3,
                           10, protonPos, 0.4);
    fitfunc->SetParLimits(2, 0, 1e6);
    fitfunc->SetParLimits(3, 0.85*pionPos, 1.15*pionPos);
    fitfunc->SetParLimits(4, 10, 1e6);
    fitfunc->SetParLimits(5, 0.1, 2);
    fitfunc->SetParLimits(6, 0, 1e5);
    fitfunc->SetParLimits(7, 0.85*protonPos, 1.15*protonPos);
    fitfunc->SetParLimits(8, 0.1, 2);

    // fit
    for (Int_t i = 0; i < 10; i++)
        if (!fitHisto->Fit(fitfunc, "RB0Q")) break;

    // reset range for second fit
    Double_t start = isMC ? fitfunc->GetParameter(3)-2.5*fitfunc->GetParameter(5)
                          : fitfunc->GetParameter(3)-3*fitfunc->GetParameter(5);
    fitfunc->SetRange(start, fitfunc->GetParameter(7) + 4*fitfunc->GetParameter(8));

    // second fit
    for (Int_t i = 0; i < 10; i++)
        if (!fitHisto->Fit(fitfunc, "RB0Q")) break;

    // get positions
    pionPos = fitfunc->GetParameter(3);
    protonPos = fitfunc->GetParameter(7);

    // correct weird fits
    if (pionPos < fitHisto->GetXaxis()->GetXmin() || pionPos > fitHisto->GetXaxis()->GetXmax())
        pionPos = (fitHisto->GetXaxis()->GetXmin() + fitHisto->GetXaxis()->GetXmax()) / 2;
    if (protonPos < fitHisto->GetXaxis()->GetXmin() || protonPos > fitHisto->GetXaxis()->GetXmax())
        protonPos = (fitHisto->GetXaxis()->GetXmin() + fitHisto->GetXaxis()->GetXmax()) / 2;

    // format and plot stuff
    if (isMC)
    {
        fPionMC = pionPos;
        fProtonMC = protonPos;
        fLinePionMC->SetPos(pionPos);
        fLineProtMC->SetPos(protonPos);
        fCanvasFit->cd(2);
    }
    else
    {
        fPionData = pionPos;
        fProtonData = protonPos;
        fLinePion->SetPos(pionPos);
        fLineProt->SetPos(protonPos);
        fCanvasFit->cd(4);
    }

    // draw histogram
    fitHisto->GetXaxis()->SetRangeUser(0, fitfunc->GetParameter(7) + 4*fitfunc->GetParameter(8));
    fitHisto->Draw("hist");
    fitfunc->Draw("same");
    if (isMC)
    {
        fLinePionMC->Draw();
        fLineProtMC->Draw();
    }
    else
    {
        fLinePion->Draw();
        fLineProt->Draw();
    }

    fCanvasFit->Update();
}

//______________________________________________________________________________
void TCCalibDeltaETrad::Fit(Int_t elem)
{
    // Perform the fit of the element 'elem'.

    Char_t tmp[256];

    // create histogram name
    sprintf(tmp, "%s_%03d", fHistoName.Data(), elem);

    // get histogram
    TH1* hmain = (TH1*) fFileManager->GetHistogram(tmp);
    if (!hmain)
    {
        Error("Init", "Main histogram does not exist!\n");
        return;
    }

    // load the MC histogram
    TH1* hmainMC = (TH1*) fMCFile->Get(tmp);
    if (!hmainMC)
    {
        Error("Init", "Could not open MC histogram!");
        return;
    }

    // get configuration
    Double_t lowLimit, highLimit;
    sprintf(tmp, "%s.Fit.Range.DetPos", GetName());
    TCReadConfig::GetReader()->GetConfigDoubleDouble(tmp, &lowLimit, &highLimit);

    // create 2D projection if necessary (data)
    if (fMainHisto) delete fMainHisto;
    if (hmain->InheritsFrom("TH3"))
    {
        sprintf(tmp, "%02d_Data_yxe", elem);
        hmain->GetZaxis()->SetRangeUser(lowLimit, highLimit);
        fMainHisto = (TH2D*) ((TH3*)hmain)->Project3D(tmp);
        fMainHisto->SetTitle(tmp);
        delete hmain;
    }
    else if (hmain->InheritsFrom("TH2"))
    {
        fMainHisto = hmain;
    }
    else
    {
        Error("Init", "Main histogram has to be either TH3 or TH2!\n");
        delete hmain;
        return;
    }

    // create 2D projection if necessary (MC)
    if (fMainHistoMC) delete fMainHistoMC;
    if (hmainMC->InheritsFrom("TH3"))
    {
        sprintf(tmp, "%02d_MC_yx", elem);
        hmainMC->GetZaxis()->SetRangeUser(lowLimit, highLimit);
        fMainHistoMC = (TH2D*) ((TH3*)hmainMC)->Project3D(tmp);
        fMainHistoMC->SetTitle(tmp);
        delete hmainMC;
    }
    else if (hmainMC->InheritsFrom("TH2"))
    {
        fMainHistoMC = hmainMC;
    }
    else
    {
        Error("Init", "Main MC histogram has to be either TH3 or TH2!\n");
        delete hmainMC;
        return;
    }

    // draw main histogram
    fCanvasFit->cd(3)->SetLogz();
    sprintf(tmp, "%s.Histo.Fit", GetName());
    TCUtils::FormatHistogram(fMainHisto, tmp);
    fMainHisto->Draw("colz");
    fCanvasFit->cd(1)->SetLogz();
    TCUtils::FormatHistogram(fMainHistoMC, tmp);
    fMainHistoMC->Draw("colz");
    fCanvasFit->Update();

    // check for sufficient statistics
    if (fMainHisto->GetEntries())
    {
        // fit the energy slices
        FitSlice((TH2*)fMainHistoMC, kTRUE);
        FitSlice((TH2*)fMainHisto, kFALSE);
    }

    // update overview
    fCanvasResult->cd(1);
    fPionPos->Draw("E1");
    fCanvasResult->cd(2);
    fProtonPos->Draw("E1");
    fCanvasResult->Update();
}

//______________________________________________________________________________
void TCCalibDeltaETrad::Calculate(Int_t elem)
{
    // Calculate the new value of the element 'elem'.

    Bool_t noval = kFALSE;

    // check if fit was performed
    if (fMainHisto->GetEntries())
    {
        // check if line position was modified by hand
        if (fLinePion->GetPos() != fPionData) fPionData = fLinePion->GetPos();
        if (fLineProt->GetPos() != fProtonData) fProtonData = fLineProt->GetPos();
        if (fLinePionMC->GetPos() != fPionMC) fPionMC = fLinePionMC->GetPos();
        if (fLineProtMC->GetPos() != fProtonMC) fProtonMC = fLineProtMC->GetPos();

        // calculate adc values of data fits
        Double_t adc_pion = fPionData/fGainOld[elem] + fPedOld[elem];
        Double_t adc_proton = fProtonData/fGainOld[elem] + fPedOld[elem];

        // calculate new values
        fPedNew[elem] = (adc_pion - (fPionMC/fProtonMC)*adc_proton)/(1. - (fPionMC/fProtonMC));
        fGainNew[elem] = fProtonMC/(adc_proton - fPedNew[elem]);

        // update overview histograms
        fPionPos->SetBinContent(elem+1, 100*(fPionData-fPionMC)/fPionMC );
        fPionPos->SetBinError(elem+1, 0.0000001);
        fProtonPos->SetBinContent(elem+1, 100*(fProtonData-fProtonMC)/fProtonMC);
        fProtonPos->SetBinError(elem+1, 0.0000001);
    }
    else
    {
        noval = kTRUE;
    }

    // user information
    printf("Element: %03d    Pion: %.2f (MC: %.2f)   Proton: %.2f (MC: %.2f)   Pedestal: %12.8f (%2.1f%%)  Gain: %12.8f (%2.1f%%)",
           elem, fPionData, fPionMC, fProtonData, fProtonMC,
           fPedNew[elem], TCUtils::GetDiffPercent(fPedOld[elem], fPedNew[elem]),
           fGainNew[elem], TCUtils::GetDiffPercent(fGainOld[elem], fGainNew[elem]));
    if (noval) printf("    -> no fit");
    printf("\n");
}

//______________________________________________________________________________
void TCCalibDeltaETrad::PrintValues()
{
    // Print out the old and new values for all elements.

    // loop over elements
    for (Int_t i = 0; i < fNelem; i++)
    {
        printf("Element: %03d    Pedestal: %12.8f    Gain: %12.8f\n",
               i, fPedNew[i], fGainNew[i]);
    }
}

//______________________________________________________________________________
void TCCalibDeltaETrad::WriteValues()
{
    // Write the obtained calibration values to the database.

    // name of pedestal data type
    TString peds = fData;
    peds.ReplaceAll("E1", "E0");

    // write values to database
    for (Int_t i = 0; i < fNset; i++)
    {

        TCMySQLManager::GetManager()->WriteParameters(peds.Data(), fCalibration.Data(), fSet[i], fPedNew, fNelem);
        TCMySQLManager::GetManager()->WriteParameters(fData.Data(), fCalibration.Data(), fSet[i], fGainNew, fNelem);
    }

    // save overview canvas
    SaveCanvas(fCanvasResult, "Overview");
}

