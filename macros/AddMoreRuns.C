/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AddMoreRuns.C                                                        //
//                                                                      //
// Add more runs to an existing calibration database.                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
void AddMoreRuns()
{
    // load CaLib
    gSystem->Load("libCaLib.so");

    // macro configuration: just change here for your needs and leave
    // the other parts of the code unchanged
    const Char_t rawfilePath[]      = "/home/werthm/loc/raw/May_18/";
    const Char_t target[]           = "LH2";
    const Int_t newFirstRun         = 0;            // 0 to keep current first run
    const Int_t newLastRun          = 0;            // 0 to keep current first run
    const Char_t calibName[]        = "LH2_May_18";

    // add more raw files to the database
    TCMySQLManager::GetManager()->AddRunFiles(rawfilePath, target, "CBTaggTAPS");
    TCMySQLManager::GetManager()->AddRunFiles(rawfilePath, target, "CBTaggTAPSPed");
    TCMySQLManager::GetManager()->AddRunFiles(rawfilePath, target, "TaggEff");

    // set new run range
    TCMySQLManager::GetManager()->ChangeCalibrationRunRange(calibName, newFirstRun, newLastRun);

    gSystem->Exit(0);
}

