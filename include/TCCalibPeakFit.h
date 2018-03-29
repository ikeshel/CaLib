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


#ifndef TCCALIBPEAKFIT_H
#define TCCALIBPEAKFIT_H

#include "TCCalib.h"
#include "TCConfig.h"

class TCLine;

class TCCalibPeakFit : public TCCalib
{

private:
    Double_t fMean;                     // mean time position
    TCLine* fLine;                      // indicator line

    virtual void Init();
    virtual void Fit(Int_t elem);
    virtual void Calculate(Int_t elem);

public:
    TCCalibPeakFit() : TCCalib(), fMean(0), fLine(0) { }
    TCCalibPeakFit(const Char_t* name, const Char_t* title, const Char_t* data,
                   Int_t nElem);
    virtual ~TCCalibPeakFit();

    virtual void WriteValues();
    virtual void PrintValues();

    ClassDef(TCCalibPeakFit, 0) // Generic peak fitting calibration class
};

class TCCalibPeakFitCB : public TCCalibPeakFit
{

public:
    TCCalibPeakFitCB()
        : TCCalibPeakFit("CB.PeakFit", "CB generic peak fitting",
                         "Data.CB.T0",
                         TCConfig::kMaxCB) { }
    virtual ~TCCalibPeakFitCB() { }

    ClassDef(TCCalibPeakFitCB, 0) // CB peak fitting class
};

class TCCalibPeakFitTAPS : public TCCalibPeakFit
{

public:
    TCCalibPeakFitTAPS()
        : TCCalibPeakFit("TAPS.PeakFit", "TAPS generic peak fitting",
                         "Data.TAPS.T0",
                         TCConfig::kMaxTAPS) { }
    virtual ~TCCalibPeakFitTAPS() { }

    ClassDef(TCCalibPeakFitTAPS, 0) // CB peak fitting class
};

#endif

