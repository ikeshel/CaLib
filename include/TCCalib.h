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

#include "TNamed.h"

class TH1;
class TF1;
class TCanvas;
class TCLine;

class TCCalib : public TNamed
{

protected:
    TString fData;                  // used calibration data
    TString fCalibration;           // calibration identifier
    Int_t fNset;                    // number of sets to be calibrated
    Int_t* fSet;                    //[fNset] array of sets to be calibrated
    TString fHistoName;             // name of the calibration histogram
    Int_t fNelem;                   // number of calibration values
    Int_t fCurrentElem;             // number of current element
    Double_t* fOldVal;              //[fNelem] old calibration value array
    Double_t* fNewVal;              //[fNelem] new calibration value array
    Double_t fAvr;                  // average value
    Double_t fAvrDiff;              // average difference to aimed value
    Int_t fNcalc;                   // number of calculated elements
    Double_t fConvergenceFactor;    // factor to control convergence

    TH1* fMainHisto;                // main histogram
    TH1* fFitHisto;                 // fitting histogram
    TF1* fFitFunc;                  // fitting function
    TCLine* fLine;                  // indicator line

    TH1* fOverviewHisto;            // overview result histogram

    TCanvas* fCanvasFit;            // canvas containing the fits
    TCanvas* fCanvasResult;         // canvas containing the results

    TTimer* fTimer;                 // slow-motion timer
    Bool_t fTimerRunning;           // timer running state

    Bool_t fIsReFit;                // re-fit flag

    Int_t fNIgnore;                 // number of elements to ignore
    Int_t* fIgnore;                 // list of elements to ignore

    virtual void Init() = 0;
    virtual void Fit(Int_t elem) = 0;
    virtual void Calculate(Int_t elem) = 0;
    void SaveCanvas(TCanvas* c, const Char_t* name);
    Bool_t IsIgnored(Int_t elem);

public:
    TCCalib() : TNamed(),
                fData(),
                fCalibration(),
                fNset(0), fSet(0), fHistoName(),
                fNelem(0), fCurrentElem(0),
                fOldVal(0), fNewVal(0),
                fAvr(0), fAvrDiff(0), fNcalc(0),
                fConvergenceFactor(1),
                fMainHisto(0), fFitHisto(0), fFitFunc(0),
                fLine(0),
                fOverviewHisto(0),
                fCanvasFit(0), fCanvasResult(0),
                fTimer(0), fTimerRunning(kFALSE),
                fIsReFit(kFALSE),
                fNIgnore(0), fIgnore(0) { }
    TCCalib(const Char_t* name, const Char_t* title,
            const Char_t* data, Int_t nElem)
        : TNamed(name, title),
          fData(data),
          fCalibration(),
          fNset(0), fSet(0), fHistoName(),
          fNelem(nElem), fCurrentElem(0),
          fOldVal(0), fNewVal(0),
          fAvr(0), fAvrDiff(0), fNcalc(0),
          fMainHisto(0), fFitHisto(0), fFitFunc(0),
          fLine(0),
          fOverviewHisto(0),
          fCanvasFit(0), fCanvasResult(0),
          fTimer(0), fTimerRunning(kFALSE),
          fIsReFit(kFALSE),
          fNIgnore(0), fIgnore(0) { }
    virtual ~TCCalib();

    virtual void WriteValues();
    virtual void PrintValues();
    virtual void PrintValuesChanged();

    void Start(const Char_t* calibration, Int_t nSet, Int_t* set);
    void ProcessAll(Int_t msecDelay = 0);
    void ProcessElement(Int_t elem, Bool_t ignorePrev = kFALSE);
    void Previous();
    void Next();
    virtual void ReFit();
    void Ignore();
    void StopProcessing();

    TString GetCalibData() { return fData; }

    void EventHandler(Int_t event, Int_t ox, Int_t oy, TObject* selected);

    ClassDef(TCCalib, 0) // Base calibration module class
};

#endif

