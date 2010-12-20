// SVN Info: $Id$

/*************************************************************************
 * Author: Irakli Keshelashvili, Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibCBTimeWalk                                                    //
//                                                                      //
// Calibration module for the CB time walk.                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef TCCALIBCBTIMEWALK_H
#define TCCALIBCBTIMEWALK_H

#include "TCanvas.h"
#include "TH2.h"
#include "TLine.h"

#include "TCCalib.h"
#include "TCFileManager.h"
#include "TCUtils.h"


class TCCalibCBTimeWalk : public TCCalib
{

private:
    TCFileManager* fFileManager;        // file manager
    Double_t* fPar0;                    // time walk parameter 0
    Double_t* fPar1;                    // time walk parameter 1
    Double_t* fPar2;                    // time walk parameter 2
    Double_t* fPar3;                    // time walk parameter 3

    virtual void Init();
    virtual void Fit(Int_t elem);
    virtual void Calculate(Int_t elem);

public:
    TCCalibCBTimeWalk();
    virtual ~TCCalibCBTimeWalk();

    virtual void Write();
    virtual void PrintValues();

    ClassDef(TCCalibCBTimeWalk, 0)   // CB time walk calibration
};

#endif

