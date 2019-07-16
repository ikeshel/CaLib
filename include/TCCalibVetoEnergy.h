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

#include "TCCalib.h"

class TCFileManager;
class TCLine;
class TFile;
class TH2;

class TCCalibVetoEnergy : public TCCalib
{

private:
    TCFileManager* fFileManager;        // file manager
    Double_t fPeak;                     // proton peak position
    Double_t fPeakMC;                   // proton MC peak position
    TCLine* fLine;                      // mean indicator line
    TCLine* fLineMC;                    // mean indicator line
    TH2* fMCHisto;                      // MC histogram
    TH1* fFitHistoMC;                   // MC fit histogram
    TF1* fFitFuncMC;                    // MC fitting function
    TFile* fMCFile;                     // MC ROOT file

    virtual void Init();
    virtual void Fit(Int_t elem);
    virtual void Calculate(Int_t elem);

    void FitSlice(TH2* h, Int_t elem, Bool_t isMC);

public:
    TCCalibVetoEnergy();
    virtual ~TCCalibVetoEnergy();

    ClassDef(TCCalibVetoEnergy, 0) // Veto energy calibration
};

#endif

