/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibDroop                                                         //
//                                                                      //
// Calibration module for droop corrections.                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef TCCALIBDROOP_H
#define TCCALIBDROOP_H

#include "TCCalib.h"
#include "TCConfig.h"

class TH2;
class TH3;
class TGraph;
class TCLine;
class TCFileManager;

class TCCalibDroop : public TCCalib
{

private:
    TCFileManager* fFileManager;        // file manager
    TH2* fProj2D;                       // dE vs E projection
    TGraph** fLinPlot;                  // correction graphs
    Int_t fNpeak;                       // number of proton peaks
    Int_t fNpoint;                      // number of points in graph
    Double_t* fPeak;                    //[fNpeak] proton peak positions
    Double_t* fTheta;                   //[fNpeak] theta positions
    TCLine** fLine;                     // mean indicator lines

    virtual void Init();
    virtual void Fit(Int_t elem);
    virtual void Calculate(Int_t elem);

    Bool_t FitHisto(Double_t* outPeak);
    void FitSlices(TH3* h, Int_t elem);

public:
    TCCalibDroop() : TCCalib(), fFileManager(0), fProj2D(0),
                     fLinPlot(0), fNpeak(0), fNpoint(0), fPeak(0), fTheta(0),
                     fLine(0) { }
    TCCalibDroop(const Char_t* name, const Char_t* title, const Char_t* data,
                 Int_t nElem);
    virtual ~TCCalibDroop();

    virtual void WriteValues();
    virtual void PrintValues();

    ClassDef(TCCalibDroop, 0) // droop correction
};

class TCCalibPIDDroop : public TCCalibDroop
{

public:
    TCCalibPIDDroop()
        : TCCalibDroop("PID.Droop", "PID droop correction",
                       "Data.PID.E0", TCConfig::kMaxPID) { }
    virtual ~TCCalibPIDDroop() { }

    ClassDef(TCCalibPIDDroop, 0) // PID droop correction
};

class TCCalibPizzaDroop : public TCCalibDroop
{

public:
    TCCalibPizzaDroop()
        : TCCalibDroop("Pizza.Droop", "Pizza droop correction",
                       "Data.Pizza.E0", TCConfig::kMaxPizza) { }
    virtual ~TCCalibPizzaDroop() { }

    ClassDef(TCCalibPizzaDroop, 0) // Pizza detector droop correction
};

#endif

