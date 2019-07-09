/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// DisplayValues.C                                                      //
//                                                                      //
// Display calibration data.                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//______________________________________________________________________________
void DisplayValues()
{
    // load CaLib
    gSystem->Load("libCaLib.so");

    // configuration
    const Int_t nSet = 1;
    const Int_t nPar = 24;
    const Char_t data[] = "Data.PID.Phi";
    //const Int_t nPar = 384;
    //const Char_t data[] = "Data.Veto.T1";
    const Char_t calibration[] = "Solid_Nov_18";

    // data array
    Double_t par[nSet][nPar];

    // read sets
    for (Int_t i = 0; i < nSet; i++)
        TCMySQLManager::GetManager()->ReadParameters(data, calibration, i, par[i], nPar);

    // display parameters
    // loop over parameters
    for (Int_t i = 0; i < nPar; i++)
    {
        // loop over sets
        for (Int_t j = 0; j < nSet; j++)
        {
            //printf("%10.3lf  ", par[j][i]);
            printf("%lf;", par[j][i]);
            //if (TCUtils::GetTAPSRing(i, 438) <= 2)
            //    printf("%lf;", 0.);
            //else
            //  printf("%lf;", par[j][i]);
        }
        printf("\n");
    }

    gSystem->Exit(0);
}

