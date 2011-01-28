// SVN Info: $Id$

/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCWriteARCalib                                                       //
//                                                                      //
// Write AcquRoot calibration files.                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TCWriteARCalib.h"

ClassImp(TCWriteARCalib)


//______________________________________________________________________________
TCWriteARCalib::TCWriteARCalib(CalibDetector_t det, const Char_t* templateFile)
{
    // Constructor.
    
    // init members
    fDetector = det;
    strcpy(fTemplate, templateFile);
}

//______________________________________________________________________________
void TCWriteARCalib::Write(const Char_t* calibFile, 
                           const Char_t* calibration, Int_t run)
{
    // Write the calibration file 'calibFile' for the run 'run' using the
    // calibration 'calibration'.
    
    Char_t line[256];
    
    // get MySQL manager
    TCMySQLManager* m = TCMySQLManager::GetManager();
    
    // read the template file
    Bool_t isTagger = kFALSE;
    if (fDetector == kDETECTOR_TAGG) isTagger = kTRUE;
    TCReadARCalib* r = new TCReadARCalib(fTemplate, isTagger);
     
    // read SG for TAPS
    TCReadARCalib* rSG = 0;
    if (fDetector == kDETECTOR_TAPS) 
        rSG = new TCReadARCalib(fTemplate, kFALSE, "TAPSSG:");

    // get the number of detectors
    Int_t nDet = r->GetNelements();
    Int_t nDetTW = r->GetNtimeWalks();
    Int_t nDetSG = 0;
    if (rSG) nDetSG = rSG->GetNelements();

    // create parameter array
    Double_t par[nDet];

    // check the detetector
    switch (fDetector)
    {
        case kDETECTOR_TAGG:
        {
            // read time offset
            if (m->ReadParametersRun(kCALIB_TAGG_T0, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetOffset(par[i]);

            break;
 
        }
        case kDETECTOR_CB:
        {
            // read time offset
            if (m->ReadParametersRun(kCALIB_CB_T0, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetOffset(par[i]);

            // read ADC gain
            if (m->ReadParametersRun(kCALIB_CB_E1, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetADCGain(par[i]);
            
            if (nDetTW)
            {
                // read time walk parameter 0
                if (m->ReadParametersRun(kCALIB_CB_WALK0, calibration, run, par, nDetTW))
                    for (Int_t i = 0; i < nDetTW; i++) r->GetTimeWalk(i)->SetPar0(par[i]);

                // read time walk parameter 1
                if (m->ReadParametersRun(kCALIB_CB_WALK1, calibration, run, par, nDetTW))
                    for (Int_t i = 0; i < nDetTW; i++) r->GetTimeWalk(i)->SetPar1(par[i]);

                // read time walk parameter 2
                if (m->ReadParametersRun(kCALIB_CB_WALK2, calibration, run, par, nDetTW))
                    for (Int_t i = 0; i < nDetTW; i++) r->GetTimeWalk(i)->SetPar2(par[i]);

                // read time walk parameter 3
                if (m->ReadParametersRun(kCALIB_CB_WALK3, calibration, run, par, nDetTW))
                    for (Int_t i = 0; i < nDetTW; i++) r->GetTimeWalk(i)->SetPar3(par[i]);
            }

            break;
        }
        case kDETECTOR_TAPS:
        {
            // read time offset
            if (m->ReadParametersRun(kCALIB_TAPS_T0, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetOffset(par[i]);

            // read TDC gain
            if (m->ReadParametersRun(kCALIB_TAPS_T1, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetTDCGain(par[i]);

            // read ADC pedestal
            if (m->ReadParametersRun(kCALIB_TAPS_LG_E0, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetPedestal(par[i]);

            // read ADC gain
            if (m->ReadParametersRun(kCALIB_TAPS_LG_E1, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetADCGain(par[i]);
            
            if (nDetSG)
            {
                // read SG ADC pedestal
                if (m->ReadParametersRun(kCALIB_TAPS_SG_E0, calibration, run, par, nDetSG))
                    for (Int_t i = 0; i < nDetSG; i++) rSG->GetElement(i)->SetPedestal(par[i]);

                // read SG ADC gain
                if (m->ReadParametersRun(kCALIB_TAPS_SG_E1, calibration, run, par, nDetSG))
                    for (Int_t i = 0; i < nDetSG; i++) rSG->GetElement(i)->SetADCGain(par[i]);
            }

            break;
        }
        case kDETECTOR_PID:
        {
            // read phi angle
            if (m->ReadParametersRun(kCALIB_PID_PHI, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetZ(par[i]);

            // read time offset
            if (m->ReadParametersRun(kCALIB_PID_T0, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetOffset(par[i]);

            // read ADC pedestal
            if (m->ReadParametersRun(kCALIB_PID_E0, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetPedestal(par[i]);

            // read ADC gain
            if (m->ReadParametersRun(kCALIB_PID_E1, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetADCGain(par[i]);

            break;
        }
        case kDETECTOR_VETO:
        {
            // read time offset
            if (m->ReadParametersRun(kCALIB_VETO_T0, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetOffset(par[i]);

            // read TDC gain
            if (m->ReadParametersRun(kCALIB_VETO_T1, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetTDCGain(par[i]);
            
            // read ADC pedestal
            if (m->ReadParametersRun(kCALIB_VETO_E0, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetPedestal(par[i]);

            // read ADC gain
            if (m->ReadParametersRun(kCALIB_VETO_E1, calibration, run, par, nDet))
                for (Int_t i = 0; i < nDet; i++) r->GetElement(i)->SetADCGain(par[i]);

            break;
        }
        case kDETECTOR_NODET:
        {
            break;
        }
    }
    
    // open the template file
    FILE* ftemp = fopen(fTemplate, "r");
    if (!ftemp)
    {   Error("Write", "Could not open template AcquRoot calibration file!");
        return;
    }

    // open the output file
    FILE* fout = fopen(calibFile, "w");
    if (!fout)
    {   Error("Write", "Could not open new AcquRoot calibration file!");
        return;
    }

    // read template file
    Char_t desc[256];
    Int_t nElem = 0;
    Int_t nElemTW = 0;
    Int_t nElemSG = 0;
    while (fgets(line, 256, ftemp))
    {   
        // read first string
        sscanf(line, "%s", desc);   
        
        // check for element line
        if(!strcmp(desc, "Element:"))
        {
            r->GetElement(nElem++)->Format(line);
            fprintf(fout, "Element: %s\n", line);
        }
        else if(!strcmp(desc, "TimeWalk:"))
        {
            r->GetTimeWalk(nElemTW++)->Format(line);
            fprintf(fout, "TimeWalk: %s\n", line);
        }
        else if(!strcmp(desc, "TAPSSG:"))
        {
            rSG->GetElement(nElemSG++)->Format(line);
            fprintf(fout, "TAPSSG: %s\n", line);
        }
        else
        {
            // write normal line
            fprintf(fout, "%s", line);
        }
    }

    // close the files
    fclose(ftemp);
    fclose(fout);

    // clean-up
    delete r;
    if (rSG) delete rSG;
}
