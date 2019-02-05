/******************************************************************************
 * Author: Thomas Strub                                                       *
 ******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// WriteNScRtoDB.C                                                            //
//                                                                            //
// Gets the number of scaler reads from histogram 'EventInfo' (bin number     //
// TCConfig::kNScREventHBin) for all runs of calibration 'calibration' and,   //
// after user confirmation, writes these values to the database.              //
//                                                                            //
// NB: Needs the correct CaLib config file.                                   //
//     Mk1: around 50 ScR/GB.                                                 //
//     Mk2: around ?? ScR/GB.                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


void WriteNScRtoDB()
{
    // Main method.

    // macro configuration (to be adapted by the user)
    const Char_t* calibration = "Solid_Oct_18";

    // get configuration ///////////////////////////////////////////////////////

    // get input file pattern from config file
    TString* filepatt = TCReadConfig::GetReader()->GetConfig("File.Input.Rootfiles");
    if (filepatt)
    {
        // check file pattern
        if (!filepatt->Contains("RUN"))
        {
            printf("Error: Error in file pattern configuration!\n");
            gSystem->Exit(1);
        }
    }
    else
    {
        printf("Error: Could not load input file pattern from configuration!\n");
        gSystem->Exit(1);
    }

    // get no. of sets from data base (use data "Data.Tagger.T0")
    Int_t nsets = TCMySQLManager::GetManager()->GetNsets("Data.Tagger.T0", calibration);
    if (!nsets)
    {
        printf("Error: No sets found!\n");
        gSystem->Exit(1);
    }

    // init final list of runs variables
    Int_t nruns = 0;
    Int_t* runs = 0;

    // init list of runs variables for the individual sets
    Int_t* nruns_set = new Int_t[nsets];
    Int_t** runs_set = new Int_t*[nsets];

    // loop over sets
    for (Int_t i = 0; i < nsets; i++)
    {
        // get list of runs for this set
        runs_set[i] = TCMySQLManager::GetManager()->GetRunsOfSet("Data.Tagger.T0", calibration, i, &nruns_set[i]);

        // check result
        if (!runs_set[i])
        {
            printf("Error: Runs of set %i not found!\n", i);

            // clean up & exit
            if (nruns_set) delete [] nruns_set;
            if (runs_set)
            {
                for (Int_t i = 0; i < nsets; i++)
                    if (runs_set[i]) delete [] runs_set[i];
                delete [] runs_set;
            }
            gSystem->Exit(1);
        }

        // add number of runs of this set to total number of runs
        nruns += nruns_set[i];
    }

    // user info
    printf("Info: %d runs in %d set(s) found.\n", nruns, nsets);

    // create and fill array of run numbers
    runs = new Int_t[nruns];
    Int_t index = 0;
    for (Int_t i = 0; i < nsets; i++)
    {
        for (Int_t j = 0; j < nruns_set[i]; j++)
            runs[index + j] = runs_set[i][j];
        index += nruns_set[i];
    }

    // clean up
    if (nruns_set) delete [] nruns_set;
    if (runs_set)
    {
        for (Int_t i = 0; i < nsets; i++)
            if (runs_set[i]) delete [] runs_set[i];
        delete [] runs_set;
    }

    // create array for number of bad scaler reads
    Int_t* n_scr = new Int_t[nruns];


    // read number of scaler reads from files //////////////////////////////////

    // init error counter
    Int_t nerrors = 0;

    // loop over runs
    for (Int_t i = 0; i < nruns; i++)
    {
        // init number of scaler reads for this run
        n_scr[i] = -1;

        // construct file name
        TString filename(*filepatt);
        filename.ReplaceAll("RUN", TString::Format("%d", runs[i]));

        // open the file
        TFile* f = TFile::Open(filename.Data());

        // check file
        if (!f)
        {
            printf("Error: Could not open file '%s'\n", filename.Data());
            nerrors++;
            continue;
        }

        // check bad file
        if (f->IsZombie())
        {
            printf("Error: Could not open file '%s'\n", filename.Data());
            if (f) delete f;
            nerrors++;
            continue;
        }

        // get 'EventInfo' histo
        TH1* h = (TH1*) f->Get("EventInfo");

        // check histo
        if (h)
        {
            // set number of scaler reads for this run
            n_scr[i] = h->GetBinContent(TCConfig::kNScREventHBin);
            printf("Info: Number of scaler reads for run '%d': %d\n", runs[i], n_scr[i]);
        }
        else
        {
            printf("Error: Histogram 'EventInfo' for run '%d' not found.\n", runs[i]);
            nerrors++;
        }

        // try to get scaler read tree
        TTree* t = (TTree*)f->Get("ScalerEvents");
        if (t)
        {
            Int_t n_scr_tree = t->GetEntries();
            if (n_scr_tree != n_scr[i])
                printf("Number of scaler reads in EventInfo and ScalerEvents differ! (%d != %d)\n",
                       n_scr_tree, n_scr[i]);
        }

        // close file
        f->Close();
        if (f) delete f;
    }

    // user info
    if (nerrors) printf("Error: %d errors detected.\n", nerrors);
    else printf("Info: All histograms found.\n");


    // ask for confirmation ////////////////////////////////////////////////////

    Char_t answer[128];
    printf("Write values to database of calibration '%s'? (yes/no) : ", calibration);
    scanf("%s", answer);
    if (strcmp(answer, "yes"))
    {
        printf("Info: Aborted.\n");
        gSystem->Exit(0);
    }

    Bool_t kWriteDefault = kFALSE;
    if (nerrors)
    {
        printf("Write default value '-1' for not existing runs? (yes/no) : ");
        scanf("%s", answer);
        if (strcmp(answer, "yes"))
            kWriteDefault = kFALSE;
        else
            kWriteDefault = kTRUE;
    }


    // save values to data base ////////////////////////////////////////////////

    // reset error counter
    nerrors = 0;

    // loop over runs
    for (Int_t i = 0; i < nruns; i++)
    {
        // check writing option
        if (!kWriteDefault && n_scr[i] == -1)
        {
            printf("Info: Skipping run '%d': %d\n", runs[i], n_scr[i]);
            continue;
        }

        // save
        if(TCMySQLManager::GetManager()->ChangeRunNScR(runs[i], n_scr[i]))
        {
            printf("Info: Write to data base: Number of scaler reads for run '%d': %d\n", runs[i], n_scr[i]);
        }
        else
        {
            nerrors++;
            printf("Error: Unable to write number of scaler reads for run '%d'\n", runs[i]);
        }
    }

    // user info
    if (nerrors)
    {
        printf("Error: %d error(s) detected.\n", nerrors);
        gSystem->Exit(1);
    }
    else
    {
        printf("Info: All values successfully written.\n");
        gSystem->Exit(0);
    }
}

