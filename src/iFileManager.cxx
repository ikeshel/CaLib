/*************************************************************************
 * Author: Dominik Werthmueller, Irakli Keshelashvili
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// iFileManager                                                         //
//                                                                      //
// Histogram building class.                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "iFileManager.hh"

ClassImp(iFileManager)


//______________________________________________________________________________
iFileManager::iFileManager(Int_t set, CalibData_t data)
{
    // Constructor.
    
    // init members
    fSet = set;
    fCalibData = data;
    fFiles = new TList();
    fFiles->SetOwner(kTRUE);

    // read input file pattern
    if (TString* f = iReadConfig::GetReader()->GetConfig("File.Input.Rootfiles"))
    {
        fInputFilePatt = *f;
        
        // check file pattern
        if (!fInputFilePatt.Contains("RUN"))
        {
            Error("iFileManager", "Error in file pattern configuration!");
            return;
        }
    }
    else
    {
        Error("iFileManager", "Could not load input file pattern from configuration!");
        return;
    }

    // build the list of files
    BuildFileList();
}

//______________________________________________________________________________
iFileManager::~iFileManager()
{
    // Destructor.

    if (fFiles) delete fFiles;
}

//______________________________________________________________________________
void iFileManager::BuildFileList()
{
    // Build the list of files belonging to the runset.
    
    // get the list of runs for this set
    Int_t nRun;
    Int_t* runs = iMySQLManager::GetManager()->GetRunsOfSet(fCalibData, fSet, &nRun);

    // loop over runs
    for (Int_t i = 0; i < nRun; i++)
    {
        // construct file name
        TString filename(fInputFilePatt);
        filename.ReplaceAll("RUN", TString::Format("%d", runs[i]));
        
        // open the file
        TFile* f = new TFile(filename.Data());
        
        // check nonexisting file
        if (!f) 
        {   
            Warning("BuildFileList", "Could not open file '%s'", filename.Data());
            continue;
        }

        // check bad file
        if (f->IsZombie())
        {
            Warning("BuildFileList", "Could not open file '%s'", filename.Data());
            continue;
        }

        // add good file to list
        fFiles->Add(f);

        // user information
        Info("BuildFileList", "%03d : added file '%s'", i, f->GetName());
    }

    // clean-up
    delete runs;
}

//______________________________________________________________________________
TH1* iFileManager::GetHistogram(const Char_t* name)
{
    // Get the summed-up histogram with name 'name'.
    // NOTE: the histogram has to be destroyed by the caller.

    TH1* hOut = 0;

    // check if there are some runs
    if (!fFiles->GetSize())
    {
        Error("GetHistogram", "ROOT file list is empty!");
        return 0;
    }
    
    // do not keep histograms in memory
    TH1::AddDirectory(kFALSE);

    // user information
    Info("GetHistogram", "Adding histogram - please wait");

    // loop over files
    TIter next(fFiles);
    TFile* f;
    Bool_t first = kTRUE;
    while ((f = (TFile*)next()))
    {
        // get histogram
        TH1* h = (TH1*) f->Get(name);

        // check if histogram is there
        if (h)
        {
            // correct destroying
            h->ResetBit(kMustCleanup);  

            // check if object is really a histogram
            if (h->InheritsFrom("TH1"))
            {
                // check if it is the first one
                if (first) hOut = (TH1*) h->Clone();      
                else hOut->Add(h);
            }
            else
            {
                Error("GetHistogram", "Object '%s' found in file '%s' is not a histogram!",
                                      name, f->GetName());
            }

            // clean-up
            delete h;
        }
        else
        {
            Warning("GetHistogram", "Histogram '%s' was not found in file '%s'",
                                    name, f->GetName());
        }
    } // loop over files

    return hOut;
}

