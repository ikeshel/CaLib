// SVN Info: $Id$

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


#ifndef TCCALIBVETOENERGY_H
#define TCCALIBVETOENERGY_H

#include "TCanvas.h"
#include "TH2.h"
#include "TH3.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSpectrum.h"

#include "TCCalib.h"
#include "TCFileManager.h"
#include "TCUtils.h"


class TCCalibVetoEnergy : public TCCalib
{

private:
    TCFileManager* fFileManager;        // file manager
    Double_t* fPed;                     // pedestal values
    Double_t* fGain;                    // gain values
    TGraph* fLinPlot;                   // linear fitting histogram
    Int_t fNpeak;                       // number of proton peaks
    Double_t* fPeak;                    //[fNpeak] proton peak positions
    Double_t* fPeakMC;                  //[fNpeak] proton MC peak positions
    TLine* fLine;                       // mean indicator line
    Int_t fDelay;                       // projection fit display delay
    TH2* fMCHisto;                      // MC histogram
    TFile* fMCFile;                     // MC ROOT file

    virtual void Init();
    virtual void Fit(Int_t elem);
    virtual void Calculate(Int_t elem);
    
    void FitSlices(TH2* h);

public:
    TCCalibVetoEnergy();
    virtual ~TCCalibVetoEnergy();

    virtual void Write();
    virtual void PrintValues();


    ClassDef(TCCalibVetoEnergy, 0) // Veto energy calibration
};

#endif
