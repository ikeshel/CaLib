/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AddRuns.C                                                            //
//                                                                      //
// Add runs to the run database.                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
void AddRuns()
{
    // load CaLib
    gSystem->Load("libCaLib.so");

    // macro configuration: just change here for your needs and leave
    // the other parts of the code unchanged
    const Char_t rawfilePath[]      = "/home/werthm/loc/raw/Jan_19";
    const Char_t target[]           = "Solid_Jan_19";

    // add raw files to the database
    TCMySQLManager::GetManager()->AddRunFiles(rawfilePath, target, "CBTaggTAPS");
    TCMySQLManager::GetManager()->AddRunFiles(rawfilePath, target, "CBTaggTAPSPed");
    TCMySQLManager::GetManager()->AddRunFiles(rawfilePath, target, "TaggEff");

    gSystem->Exit(0);
}

