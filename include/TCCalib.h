// SVN Info: $Id$

/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalib                                                              //
//                                                                      //
// Abstract calibration module class.                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef TCCALIB_H
#define TCCALIB_H

#include "TString.h"
#include "TError.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TTimer.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TStyle.h"

#include "TCMySQLManager.h"
#include "TCReadConfig.h"


class TCCalib : public TNamed
{

protected:
    CalibData_t fData;          // used calibration data
    Int_t fSet;                 // set to be calibrated
    TString fHistoName;         // name of the calibration histogram
    Int_t fNelem;               // number of calibration values
    Int_t fCurrentElem;         // number of current element
    Double_t* fOldVal;          //[fNelem] old calibration value array
    Double_t* fNewVal;          //[fNelem] new calibration value array
    
    TH1* fMainHisto;            // main histogram 
    TH1* fFitHisto;             // fitting histogram
    TF1* fFitFunc;              // fitting function

    TH1* fOverviewHisto;        // overview result histogram

    TCanvas* fCanvasFit;        // canvas containing the fits
    TCanvas* fCanvasResult;     // canvas containing the results
    
    TTimer* fTimer;             // slow-motion timer

    virtual void Init() = 0;
    virtual void Fit(Int_t elem) = 0;
    virtual void Calculate(Int_t elem) = 0;

public:
    TCCalib() : TNamed(),
                fData(kCALIB_NODATA), 
                fSet(0), fHistoName(), 
                fNelem(0), fCurrentElem(0),
                fOldVal(0), fNewVal(0),
                fMainHisto(0), fFitHisto(0), fFitFunc(0),
                fOverviewHisto(0),
                fCanvasFit(0), fCanvasResult(0), 
                fTimer(0) { }
    TCCalib(const Char_t* name, const Char_t* title, CalibData_t data,
            Int_t nElem) 
        : TNamed(name, title),
           fData(data), 
           fSet(0), fHistoName(), 
           fNelem(nElem), fCurrentElem(0),
           fOldVal(0), fNewVal(0),
           fMainHisto(0), fFitHisto(0), fFitFunc(0),
           fOverviewHisto(0),
           fCanvasFit(0), fCanvasResult(0), 
           fTimer(0) { }
    virtual ~TCCalib();
    
    virtual void Write();
    virtual void PrintValues();

    void Start(Int_t set);
    void ProcessAll(Int_t msecDelay = 0);
    void ProcessElement(Int_t elem);
    void Previous();
    void Next();
    void StopProcessing();
    
    CalibData_t GetCalibData() const { return fData; }

    ClassDef(TCCalib, 0)         // Base calibration module class
};

#endif

