/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// WriteCFDThrTAPS.C                                                    //
//                                                                      //
// Write TAPS CFD thresholds using th current calibration and a file    //
// with raw ADC threshold positions.                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
void WriteCFDThrTAPS(const Char_t* thrFile)
{
    // load CaLib
    gSystem->Load("libCaLib.so");

    // configuration
    const Char_t* calibration = "Solid_Nov_18";
    const Int_t e_calib_set   = 0;
    const Int_t cfd_calib_set = 0;

    // try to open the file
    FILE* thr_file = fopen(thrFile, "r");
    if (!thr_file)
    {
        printf("Error: Could not open file '%s'!\n", thrFile);
        gSystem->Exit(0);
    }

    // read the file
    Double_t thr_pos[510];
    Char_t line[256];
    Int_t nElem = 0;
    while (fgets(line, 256, thr_file))
        sscanf(line, "%lf", &thr_pos[nElem++]);
    printf("Found %d element thresholds\n", nElem);

    // close input file
    fclose(thr_file);

    // read calibration from database
    Double_t ped[nElem];
    Double_t gain[nElem];
    Double_t cfd_old[nElem];
    TCMySQLManager::GetManager()->ReadParameters("Data.TAPS.LG.E0", calibration, e_calib_set, ped, nElem);
    TCMySQLManager::GetManager()->ReadParameters("Data.TAPS.LG.E1", calibration, e_calib_set, gain, nElem);
    TCMySQLManager::GetManager()->ReadParameters("Data.TAPS.CFD", calibration, e_calib_set, cfd_old, nElem);

    // calculate new thresholds
    Double_t cfd_new[nElem];
    for (Int_t i = 0; i < nElem; i++)
    {
        if (cfd_old[i] == 9999)
            cfd_new[i] = cfd_old[i];
        else
            cfd_new[i] = gain[i] * (thr_pos[i] - ped[i]);

        // user info
        printf("Elem: %03d    Old CFD: %8.3f    New CFD: %8.3f    Change: %5.1f%%\n",
               i, cfd_old[i], cfd_new[i], (cfd_new[i]-cfd_old[i])/cfd_old[i]*100.);
    }

    // ask user
    Char_t answer[128];
    printf("Write to set %d of TAPS CFD of calibration '%s'? (yes/no) : ",
           cfd_calib_set, calibration);
    scanf("%s", answer);
    if (strcmp(answer, "yes"))
    {
        printf("Aborted.\n");
        gSystem->Exit(0);
    }

    // write to database
    TCMySQLManager::GetManager()->WriteParameters("Data.TAPS.CFD", calibration, cfd_calib_set, cfd_new, nElem);

    gSystem->Exit(0);
}

