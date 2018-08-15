/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// WriteTDCGainsTAPS.C                                                  //
//                                                                      //
// Write TAPS TDC gain calibration from a TAPSMaintain input file to    //
// the database.                                                        //
// Assumes the TAPS 2009 setup.                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TCUtils.h"
#endif


//______________________________________________________________________________
void WriteTDCGainsTAPS()
{
    // load CaLib
    gSystem->Load("libCaLib.so");

    // configuration
    //const Char_t* inputFile   = "tdc.dat";
    //const Char_t* data        = "Data.TAPS.T1";
    //const Int_t nElem         = 438;
    const Char_t* inputFile   = "tdc.dat";
    const Char_t* data        = "Data.Veto.T1";
    const Int_t nElem         = 384;
    const Char_t* calibration = "LH2_May_18";
    const Int_t set           = 0;

    // read file via TGraph
    TGraph* dg = new TGraph(inputFile);

    // remove ring 1 and 2 elements
    for (Int_t i = dg->GetN(); i >= 0; i--)
    {
        if (TCUtils::GetTAPSRing(i, 384) <= 2)
            dg->RemovePoint(i);
    }
    printf("Processed gains for %d elements\n", dg->GetN());

    // prepare parameter array
    Double_t par[nElem];
    Int_t npar = 0;
    for (Int_t i = 0; i < nElem; i++)
    {
        if (TCUtils::GetTAPSRing(i, nElem) <= 2)
            par[i] = 0.1;
        else
            par[i] = dg->GetY()[npar++];
    }

    // show parameters
    for (Int_t i = 0; i < nElem; i++) printf("Par %03d : %lf\n", i, par[i]);

    // ask user
    Char_t answer[128];
    printf("Write to set %d of data '%s' of calibration '%s'? (yes/no) : ",
           set, data, calibration);
    scanf("%s", answer);
    if (strcmp(answer, "yes"))
    {
        printf("Aborted.\n");
        gSystem->Exit(0);
    }

    // write to database
    TCMySQLManager::GetManager()->WriteParameters(data, calibration, set, par, nElem);

    // clean-up
    delete dg;

    gSystem->Exit(0);
}

