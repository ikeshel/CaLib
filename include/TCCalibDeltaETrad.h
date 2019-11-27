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


#ifndef TCCALIBDELTAETRAD_H
#define TCCALIBDELTAETRAD_H

#include "TCCalib.h"
#include "TCConfig.h"
#include "TCReadConfig.h"

class TCLine;
class TH2;
class TFile;
class TCFileManager;

class TCCalibDeltaETrad : public TCCalib
{

private:
    TCFileManager* fFileManager;        // file manager
    Double_t* fPedOld;                  // old pedestal values
    Double_t* fPedNew;                  // new pedestal values
    Double_t* fGainOld;                 // old gain values
    Double_t* fGainNew;                 // new gain values
    Double_t fPionMC;                   // pion position in simulation
    Double_t fProtonMC;                 // proton position in simulation
    Double_t fPionData;                 // pion position in data
    Double_t fProtonData;               // proton position in data
    TH1* fPionPos;                      // pion peak position histogram
    TH1* fProtonPos;                    // proton peak position histogram
    TCLine* fLinePion;                  // mean indicator line
    TCLine* fLineProt;                  // mean indicator line
    TCLine* fLinePionMC;                // mean indicator line
    TCLine* fLineProtMC;                // mean indicator line
    TFile* fMCFile;                     // MC ROOT file
    TH1* fMainHistoMC;                  // MC histogram
    TH1* fFitHistoMC;                   // MC histogram
    TF1* fFitFuncMC;                    // MC fitting function

    virtual void Init();
    virtual void Fit(Int_t elem);
    virtual void Calculate(Int_t elem);

    static Double_t FitFunc(Double_t* x, Double_t* par);

    void FitSlice(TH2* h, Bool_t isMC);

public:
    TCCalibDeltaETrad() : TCCalib(), fFileManager(0),
                          fPedOld(0), fPedNew(0), fGainOld(0), fGainNew(0),
                          fPionMC(0), fProtonMC(0), fPionData(0), fProtonData(0),
                          fPionPos(0), fProtonPos(0),
                          fLinePion(0), fLineProt(0),
                          fLinePionMC(0), fLineProtMC(0),
                          fMCFile(0), fMainHistoMC(0), fFitHistoMC(0),
                          fFitFuncMC(0) { }
    TCCalibDeltaETrad(const Char_t* name, const Char_t* title, const Char_t* data,
                      Int_t nElem);
    virtual ~TCCalibDeltaETrad();

    virtual void WriteValues();
    virtual void PrintValues();

    static Double_t fgGNP;

    ClassDef(TCCalibDeltaETrad, 0) // DeltaE energy calibration (traditional method)
};

class TCCalibPIDEnergyTrad : public TCCalibDeltaETrad
{

public:
    TCCalibPIDEnergyTrad()
        : TCCalibDeltaETrad("PID.Energy.Trad", "PID energy calibration (traditional)",
                            "Data.PID.E1", TCConfig::kMaxPID) { }
    virtual ~TCCalibPIDEnergyTrad() { }

    ClassDef(TCCalibPIDEnergyTrad, 0) // PID energy calibration class
};

class TCCalibVetoEnergyTrad : public TCCalibDeltaETrad
{

public:
    TCCalibVetoEnergyTrad()
        : TCCalibDeltaETrad("Veto.Energy.Trad", "Veto energy calibration (traditional)",
                            "Data.Veto.E1", TCReadConfig::GetReader()->GetConfigInt("Veto.Elements")) { }
    virtual ~TCCalibVetoEnergyTrad() { }

    ClassDef(TCCalibVetoEnergyTrad, 0) // Veto energy calibration class
};

class TCCalibPizzaEnergyTrad : public TCCalibDeltaETrad
{

public:
    TCCalibPizzaEnergyTrad()
        : TCCalibDeltaETrad("Pizza.Energy.Trad", "Pizza energy calibration (traditional)",
                            "Data.Pizza.E1", TCConfig::kMaxPizza) { }
    virtual ~TCCalibPizzaEnergyTrad() { }

    ClassDef(TCCalibPizzaEnergyTrad, 0) // Pizza detector energy calibration class
};

#endif

