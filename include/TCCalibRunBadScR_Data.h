/************************************************************************
 * Authors: Thomas Strub                                                *
 ************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibRunBadScR_Data                                                //
//                                                                      //
// Beamtime calibration module class for run by run bad scaler reads    //
// calibration.                                                         //
//                                                                      //
// Have fun!                                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef TCCALIBRUNBADSCR_DATA_H
#define TCCALIBRUNBADSCR_DATA_H

#include "TCCalibRunBadScR.h"

class TCCalibRunBadScR_NaI : public TCCalibRunBadScR
{

public:
    TCCalibRunBadScR_NaI()
      : TCCalibRunBadScR("BadScR.NaI", "Bad scaler read calibration (NaI)", "Data.Run.BadScR.NaI", kTRUE) { }
    virtual ~TCCalibRunBadScR_NaI() { }

    //Bool_t IsTrueCalib() { return kTRUE; }

    ClassDef(TCCalibRunBadScR_NaI, 0) // NaI bad scaler read calibration class
};

class TCCalibRunBadScR_PID : public TCCalibRunBadScR
{

public:
    TCCalibRunBadScR_PID()
      : TCCalibRunBadScR("BadScR.PID", "Bad scaler read calibration (PID)", "Data.Run.BadScR.PID", kTRUE) { }
    virtual ~TCCalibRunBadScR_PID() { }

    //virtual Bool_t IsTrueCalib() const { return kTRUE; }

    ClassDef(TCCalibRunBadScR_PID, 0) // PID bad scaler read calibration class
};

class TCCalibRunBadScR_MWPC : public TCCalibRunBadScR
{

public:
    TCCalibRunBadScR_MWPC()
      : TCCalibRunBadScR("BadScR.MWPC", "Bad scaler read calibration (MWPC)", "Data.Run.BadScR.MWPC", kTRUE) { }
    virtual ~TCCalibRunBadScR_MWPC() { }

    //virtual Bool_t IsTrueCalib() const { return kTRUE; }

    ClassDef(TCCalibRunBadScR_MWPC, 0) // MWPC bad scaler read calibration class
};

class TCCalibRunBadScR_BaF2PWO : public TCCalibRunBadScR
{

public:
    TCCalibRunBadScR_BaF2PWO()
      : TCCalibRunBadScR("BadScR.BaF2PWO", "Bad scaler read calibration (BaF2PWO)", "Data.Run.BadScR.BaF2PWO", kTRUE) { }
    virtual ~TCCalibRunBadScR_BaF2PWO() { }

    //virtual Bool_t IsTrueCalib() const { return kTRUE; }

    ClassDef(TCCalibRunBadScR_BaF2PWO, 0) // BaF2PWO bad scaler read calibration class
};

class TCCalibRunBadScR_BaF2 : public TCCalibRunBadScR
{

public:
    TCCalibRunBadScR_BaF2()
      : TCCalibRunBadScR("BadScR.BaF2", "Bad scaler read calibration (BaF2)", "Data.Run.BadScR.BaF2", kTRUE) { }
    virtual ~TCCalibRunBadScR_BaF2() { }

    //virtual Bool_t IsTrueCalib() const { return kTRUE; }

    ClassDef(TCCalibRunBadScR_BaF2, 0) // BaF2 bad scaler read calibration class
};

class TCCalibRunBadScR_PWO : public TCCalibRunBadScR
{

public:
    TCCalibRunBadScR_PWO()
      : TCCalibRunBadScR("BadScR.PWO", "Bad scaler read calibration (PWO)", "Data.Run.BadScR.PWO", kTRUE) { }
    virtual ~TCCalibRunBadScR_PWO() { }

    //virtual Bool_t IsTrueCalib() const { return kTRUE; }

    ClassDef(TCCalibRunBadScR_PWO, 0) // PWO bad scaler read calibration class
};

class TCCalibRunBadScR_Veto : public TCCalibRunBadScR
{

public:
    TCCalibRunBadScR_Veto()
      : TCCalibRunBadScR("BadScR.Veto", "Bad scaler read calibration (Veto)", "Data.Run.BadScR.Veto", kTRUE) { }
    virtual ~TCCalibRunBadScR_Veto() { }

    //virtual Bool_t IsTrueCalib() const { return kTRUE; }

    ClassDef(TCCalibRunBadScR_Veto, 0) // Veto bad scaler read calibration class
};

class TCCalibRunBadScR_Pizza : public TCCalibRunBadScR
{

public:
    TCCalibRunBadScR_Pizza()
      : TCCalibRunBadScR("BadScR.Pizza", "Bad scaler read calibration (Pizza)", "Data.Run.BadScR.Pizza", kTRUE) { }
    virtual ~TCCalibRunBadScR_Pizza() { }

    //virtual Bool_t IsTrueCalib() const { return kTRUE; }

    ClassDef(TCCalibRunBadScR_Pizza, 0) // Pizza bad scaler read calibration class
};

class TCCalibRunBadScR_Ladder : public TCCalibRunBadScR
{

public:
    TCCalibRunBadScR_Ladder()
      : TCCalibRunBadScR("BadScR.Ladder", "Bad scaler read calibration (Ladder)", "Data.Run.BadScR.Ladder", kTRUE) { }
    virtual ~TCCalibRunBadScR_Ladder() { }

    //virtual Bool_t IsTrueCalib() const { return kTRUE; }

    ClassDef(TCCalibRunBadScR_Ladder, 0) // Ladder bad scaler read calibration class
};

class TCCalibRunBadScR_LadderScalers : public TCCalibRunBadScR
{

public:
    TCCalibRunBadScR_LadderScalers()
      : TCCalibRunBadScR("BadScR.LadderScalers", "Bad scaler read calibration (Ladder Scalers)", "Data.Run.BadScR.LadderScalers", kTRUE) { }
    virtual ~TCCalibRunBadScR_LadderScalers() { }

    //virtual Bool_t IsTrueCalib() const { return kTRUE; }

    ClassDef(TCCalibRunBadScR_LadderScalers, 0) // Ladder scalers bad scaler read calibration class
};

#endif

