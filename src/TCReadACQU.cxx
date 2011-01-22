// SVN Info: $Id$

/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCReadACQU                                                           //
//                                                                      //
// Read raw ACQU MK1 files.                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TCReadACQU.h"

ClassImp(TCReadACQU)


//______________________________________________________________________________
TCReadACQU::TCReadACQU(const Char_t* path) 
{
    // Constructor using the path of the raw files 'path'.
    
    // init members
    fPath = new Char_t[256];
    fFiles = new TList();
    fFiles->SetOwner(kTRUE);
    
    // copy path
    strcpy(fPath, path);

    // read all files
    ReadFiles();
}

//______________________________________________________________________________
TCReadACQU::~TCReadACQU()
{
    // Destructor.
    
    if (fPath) delete [] fPath;
    if (fFiles) delete fFiles;
}

//______________________________________________________________________________
void TCReadACQU::ReadFiles()
{
    // Read all raw files.

    // user information
    Info("ReadFiles", "Looking for ACQU raw files in '%s'", fPath);
    
    // try to get directory content
    TSystemDirectory dir("rawdir", fPath);
    TList* list = dir.GetListOfFiles();
    if (!list)
    {
        Error("ReadFiles", "'%s' is not a directory!", fPath);
        return;
    }

    // loop over directory content
    TIter next(list);
    TSystemFile* f;
    while ((f = (TSystemFile*)next()))
    {
        // look for ACQU raw files
        TString str(f->GetName());
        if (str.BeginsWith("CB_") && str.EndsWith(".dat"))
        {
            // user information
            Info("ReadFiles", "Reading '%s/%s'", fPath, f->GetName());
            
            // create file object
            TCACQUFile* acqufile = new TCACQUFile();
            acqufile->ReadFile(fPath, f->GetName());
            
            // check file 
            if (!acqufile->IsGoodStartMarker())
            {
                Warning("ReadFiles", "Bad file header found in '%s/%s'", fPath, f->GetName());
            }
            
            // add file to list
            fFiles->Add(acqufile);
        }
    }

    // clean-up
    delete list;
}
