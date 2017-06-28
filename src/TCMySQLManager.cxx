/*************************************************************************
 * Author: Dominik Werthmueller, Irakli Keshelashvili, Thomas Strub
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCMySQLManager                                                       //
//                                                                      //
// This class handles all the communication with the MySQL server.      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <fstream>

#include "THashList.h"
#include "TError.h"
#include "TSystem.h"
#include "TSQLServer.h"
#include "TSQLRow.h"
#include "TSQLResult.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"

#include "TCMySQLManager.h"
#include "TCReadConfig.h"
#include "TCReadACQU.h"
#include "TCReadARCalib.h"
#include "TCACQUFile.h"
#include "TCCalibData.h"
#include "TCCalibType.h"
#include "TCBadScRElement.h"
#include "TCContainer.h"

ClassImp(TCMySQLManager)

// init static class members
TCMySQLManager* TCMySQLManager::fgMySQLManager = 0;

//______________________________________________________________________________
TCMySQLManager::TCMySQLManager()
{
    // Constructor.

    fDB = 0;
    fDBType = kNoType;
    fSilence = kFALSE;
    fData = new THashList();
    fData->SetOwner(kTRUE);
    fTypes = new THashList();
    fTypes->SetOwner(kTRUE);

    // read CaLib data
    if (!ReadCaLibData())
    {
        if (!fSilence) Error("TCMySQLManager", "Could not read the CaLib data definitions!");
        return;
    }

    // read CaLib types
    if (!ReadCaLibTypes())
    {
        if (!fSilence) Error("TCMySQLManager", "Could not read the CaLib type definitions!");
        return;
    }

    //
    // get database configuration
    //

    TString* strDBFile;

    // get database file (sqlite)
    if ((strDBFile = TCReadConfig::GetReader()->GetConfig("DB.File")))
    {
        // open connection to sqlite server
        Char_t szMySQL[200];
        Char_t* exp = gSystem->ExpandPathName(strDBFile->Data());
        sprintf(szMySQL, "sqlite://%s", exp);
        fDB = TSQLServer::Connect(szMySQL, "", "");

        // check DB connection
        if (!fDB)
        {
            if (!fSilence) Error("TCMySQLManager", "Cannot connect to the database '%s'!",
                                 exp);
            delete exp;
            return;
        }
        else if (fDB->IsZombie())
        {
            if (!fSilence) Error("TCMySQLManager", "Cannot connect to the database '%s'!",
                                 exp);
            delete exp;
            return;
        }
        else
        {
            if (!fSilence) Info("TCMySQLManager", "Connected to the database '%s' using CaLib %s",
                                strDBFile->Data(), TCConfig::kCaLibVersion);
            fDBType = kSQLite;
            delete exp;
        }
    }
    else
    {
        TString* strDBHost;
        TString* strDBName;
        TString* strDBUser;
        TString* strDBPass;

        // get database hostname
        if (!(strDBHost = TCReadConfig::GetReader()->GetConfig("DB.Host")))
        {
            if (!fSilence) Error("TCMySQLManager", "Database host not included in configuration file!");
            return;
        }

        // get database name
        if (!(strDBName = TCReadConfig::GetReader()->GetConfig("DB.Name")))
        {
            if (!fSilence) Error("TCMySQLManager", "Database name not included in configuration file!");
            return;
        }

        // get database user
        if (!(strDBUser = TCReadConfig::GetReader()->GetConfig("DB.User")))
        {
            if (!fSilence) Error("TCMySQLManager", "Database user not included in configuration file!");
            return;
        }

        // get database password
        if (!(strDBPass = TCReadConfig::GetReader()->GetConfig("DB.Pass")))
        {
            if (!fSilence) Error("TCMySQLManager", "Database password not included in configuration file!");
            return;
        }

        // open connection to MySQL server on localhost
        Char_t szMySQL[200];
        sprintf(szMySQL, "mysql://%s/%s", strDBHost->Data(), strDBName->Data());
        fDB = TSQLServer::Connect(szMySQL, strDBUser->Data(), strDBPass->Data());

        // check DB connection
        if (!fDB)
        {
            if (!fSilence) Error("TCMySQLManager", "Cannot connect to the database '%s' on '%s@%s'!",
                                 strDBName->Data(), strDBUser->Data(), strDBHost->Data());
            return;
        }
        else if (fDB->IsZombie())
        {
            if (!fSilence) Error("TCMySQLManager", "Cannot connect to the database '%s' on '%s@%s'!",
                                 strDBName->Data(), strDBUser->Data(), strDBHost->Data());
            return;
        }
        else
        {
            if (!fSilence) Info("TCMySQLManager", "Connected to the database '%s' on '%s@%s' using CaLib %s",
                                strDBName->Data(), strDBUser->Data(), strDBHost->Data(), TCConfig::kCaLibVersion);
            fDBType = kMySQL;
        }
    }
}

//______________________________________________________________________________
TCMySQLManager::~TCMySQLManager()
{
    // Destructor.

    // close DB
    if (fDB) delete fDB;
    if (fData) delete fData;
    if (fTypes) delete fTypes;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ReadCaLibData()
{
    // Read the CaLib data definitions from the configuration file.
    // Return kTRUE on success, kFALSE if an error occurred.

    // try to get the CaLib source path from the shell variable CALIB
    // otherwise use the current directory
    TString path;
    if (gSystem->Getenv("CALIB")) path = gSystem->Getenv("CALIB");
    else path = gSystem->pwd();
    path.Append("/data/data.def");

    // open the file
    std::ifstream infile;
    infile.open(path.Data());

    // check if file is open
    if (!infile.is_open())
    {
        if (!fSilence) Error("ReadCaLibData", "Could not open data definition file '%s'", path.Data());
        return kFALSE;
    }
    else
    {
        // read the file
        while (infile.good())
        {
            TString line;
            line.ReadLine(infile);

            // trim line
            line.Remove(TString::kBoth, ' ');

            // skip comments and empty lines
            if (line.BeginsWith("#") || line == "") continue;
            else
            {
                // tokenize the string
                TObjArray* token = line.Tokenize(",");

                // some variables for data extraction
                Int_t count = 0;
                TString name;
                TString title;
                TString table;
                Int_t size = 0;

                // iterate over tokens
                TIter next(token);
                TObjString* s;
                while ((s = (TObjString*)next()))
                {
                    // get the string and trim it
                    TString str(s->GetString());
                    str.Remove(TString::kBoth, ' ');

                    // set data
                    switch (count)
                    {
                        case 0:
                            name = str;
                            break;
                        case 1:
                            title = str;
                            break;
                        case 2:
                            table = str;
                            break;
                        case 3:
                            size = atoi(str.Data());
                            break;
                    }

                    count++;
                }

                // clean-up
                delete token;

                // check number of tokens
                if (count == 4)
                {
                    // create a new CaLib data class
                    TCCalibData* data = new TCCalibData(name, title, size);
                    data->SetTableName(table);

                    // add it to list
                    fData->Add(data);
                }
                else
                {
                    if (!fSilence) Error("ReadCaLibData", "Error in CaLib data definition '%s'!", line.Data());
                    return kFALSE;
                }
            }
        }
    }

    // close the file
    infile.close();

    // user information
    if (!fSilence) Info("ReadCaLibData", "Registered %d data definitions", fData->GetSize());

    return kTRUE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ReadCaLibTypes()
{
    // Read the CaLib type definitions from the configuration file.
    // Return kTRUE on success, kFALSE if an error occurred.

    // try to get the CaLib source path from the shell variable CALIB
    // otherwise use the current directory
    TString path;
    if (gSystem->Getenv("CALIB")) path = gSystem->Getenv("CALIB");
    else path = gSystem->pwd();
    path.Append("/data/types.def");

    // open the file
    std::ifstream infile;
    infile.open(path.Data());

    // check if file is open
    if (!infile.is_open())
    {
        if (!fSilence) Error("ReadCaLibTypes", "Could not open types definition file '%s'", path.Data());
        return kFALSE;
    }
    else
    {
        // read the file
        while (infile.good())
        {
            TString line;
            line.ReadLine(infile);

            // trim line
            line.Remove(TString::kBoth, ' ');

            // skip comments and empty lines
            if (line.BeginsWith("#") || line == "") continue;
            else
            {
                // tokenize the string
                TObjArray* token = line.Tokenize(",");

                // some variables for data extraction
                Int_t count = 0;
                TString name;
                TString title;
                TCCalibType* type = 0;

                // iterate over tokens
                TIter next(token);
                TObjString* s;
                while ((s = (TObjString*)next()))
                {
                    // get the string and trim it
                    TString str(s->GetString());
                    str.Remove(TString::kBoth, ' ');

                    // set name
                    if (count == 0) name = str;
                    else if (count == 1) title = str;
                    else
                    {
                        // create type object if necessary
                        if (!type) type = new TCCalibType(name, title);

                        // get calibration data
                        TCCalibData* d = GetCalibData(str.Data());

                        // add calibration data
                        if (d) type->AddData(d);
                        else
                        {
                            if (!fSilence) Error("ReadCaLibTypes", "CaLib data '%s' was not found!", str.Data());
                            return kFALSE;
                        }
                    }

                    count++;
                }

                // clean-up
                delete token;

                // check number of tokens
                if (count > 2)
                {
                    // add type
                    fTypes->Add(type);
                }
                else
                {
                    if (!fSilence) Error("ReadCaLibTypes", "Error in CaLib type definition '%s'!", line.Data());
                    return kFALSE;
                }
            }
        }
    }

    // close the file
    infile.close();

    // user information
    if (!fSilence) Info("ReadCaLibTypes", "Registered %d type definitions", fTypes->GetSize());

    return kTRUE;
}

//______________________________________________________________________________
TSQLResult* TCMySQLManager::SendQuery(const Char_t* query)
{
    // Send a query to the database and return the result.

    // check server connection
    if (!IsConnected())
    {
        if (!fSilence) Error("SendQuery", "No connection to the database!");
        return 0;
    }

    // execute query
    return fDB->Query(query);
}

//______________________________________________________________________________
const Char_t* TCMySQLManager::GetDBName() const
{
    // Return the database name.

    return fDB ? fDB->GetDB() : 0;
}

//______________________________________________________________________________
const Char_t* TCMySQLManager::GetDBHost() const
{
    // Return the database hostname.

    return fDB ? fDB->GetHost() : 0;
}

//______________________________________________________________________________
TCCalibType* TCMySQLManager::GetCalibType(const Char_t* type) const
{
    // Returns a pointer to the calibration type 'type'.
    // If no calibration type named 'type' is found 0 is returned.

    // get calibration type
    TCCalibType* t = (TCCalibType*) fTypes->FindObject(type);
    if (!t)
    {
        if (!fSilence) Error("GetCalibType", "Calibration type named '%s' is not defined!", type);
            return 0;
    }

    return t;
}

//______________________________________________________________________________
TCCalibData* TCMySQLManager::GetCalibData(const Char_t* data) const
{
    // Returns a pointer to the calibration data 'data'.
    // If no calibration data named 'data' is found 0 is returned.

    // get calibration data
    TCCalibData* d = (TCCalibData*) fData->FindObject(data);
    if (!d)
    {
        if (!fSilence) Error("GetCalibData", "Calibration data named '%s' is not defined!", data);
            return 0;
    }

    return d;
}

//______________________________________________________________________________
TCMySQLManager* TCMySQLManager::GetManager()
{
    // Return a pointer to the static instance of this class.

    if (!fgMySQLManager) fgMySQLManager = new TCMySQLManager();
    if (!fgMySQLManager->IsConnected())
    {
        Error("GetManager", "No connection to the database!");
        delete fgMySQLManager;
        return 0;
    }
    else return fgMySQLManager;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::SendExec(const Char_t* sql)
{
    // Send the SQL command 'sql' without returned result to the database.
    // Return kTRUE on success, otherwise kFALSE.

    // check server connection
    if (!IsConnected())
    {
        if (!fSilence) Error("SendExec", "No connection to the database!");
        return kFALSE;
    }

    // execute command
    return fDB->Exec(sql);
}

//______________________________________________________________________________
Bool_t TCMySQLManager::IsConnected()
{
    // Check if the connection to the database is open.

    if (!fDB)
    {
        if (!fSilence) Error("IsConnected", "Cannot access database!");
        return kFALSE;
    }
    else
        return kTRUE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::SearchTable(const Char_t* data, Char_t* outTableName)
{
    // Search the table name of the calibration quantity 'data' and write it
    // to 'outTableName'.
    // Return kTRUE when a true table was found, otherwise kFALSE.

    // get CaLib data
    TCCalibData* d = GetCalibData(data);

    // check data
    if (d)
    {
        // copy table name
        strcpy(outTableName, d->GetTableName());
        return kTRUE;
    }
    else return kFALSE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::SearchRunEntry(Int_t run, const Char_t* name, TString& outInfo)
{
    // Search the information 'name' for the run 'run' and write it to 'outInfo'.
    // Return kTRUE when the information was found, otherwise kFALSE.

    TString query;

    // create the query
    query.Form("SELECT %s FROM %s "
               "WHERE run = %d",
               name, TCConfig::kCalibMainTableName, run);

    // read from database
    TSQLResult* res = SendQuery(query.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("SearchRunEntry", "Could not find the information '%s' for run %d!",
                                               name, run);
        return kFALSE;
    }

    // get row
    TSQLRow* row = res->Next();

    // check row
    if (!row)
    {
        if (!fSilence) Error("SearchRunEntry", "Could not find the information '%s' for run %d!",
                                               name, run);
        delete res;
        return kFALSE;
    }

    // write the information
    const Char_t* field = row->GetField(0);
    if (!field) outInfo = "";
    outInfo = field;

    // clean-up
    delete row;
    delete res;

    return kTRUE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::SearchSetEntry(const Char_t* data, const Char_t* calibration, Int_t set,
                                      const Char_t* name, Char_t* outInfo)
{
    // Search the information 'name' for the calibration identifier 'calibration' and
    // the calibration data 'data' for the set number 'set' and write it to 'outInfo'.
    // Return kTRUE when the information was found, otherwise kFALSE.

    TString query;
    Char_t table[256];

    // check for data
    if (!GetCalibData(data)) return kFALSE;

    // get the data table
    if (!SearchTable(data, table))
    {
        if (!fSilence) Error("SearchSetEntry", "No data table for '%s' in calibration '%s' found!",
                                               data, calibration);
        return kFALSE;
    }

    // create the query
    query.Form("SELECT %s FROM %s WHERE "
               "calibration = '%s' "
               "ORDER BY first_run ASC LIMIT 1 OFFSET %d",
               name, table, calibration, set);

    // read from database
    TSQLResult* res = SendQuery(query.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("SearchSetEntry", "No runset %d found in table '%s' of '%s' in calibration '%s'!",
                                               set, table, data, calibration);
        return kFALSE;
    }

    // get data
    TSQLRow* row = res->Next();

    // check row
    if (!row)
    {
        if (!fSilence) Error("SearchSetEntry", "No runset %d found in table '%s' of '%s' in calibration '%s'!",
                                               set, table, data, calibration);
        return kFALSE;
    }

    // extract data
    const Char_t* d = row->GetField(0);
    if (!d) strcpy(outInfo, "");
    else strcpy(outInfo, d);

    // clean-up
    delete row;
    delete res;

    return kTRUE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeRunEntries(Int_t first_run, Int_t last_run,
                                        const Char_t* name, const Char_t* value)
{
    // Change the run entry 'name' for the runs 'first_run' to 'last_run' to 'value'.

    TString query;

    // check if first run is smaller than last run
    if (first_run > last_run)
    {
        if (!fSilence) Error("ChangeRunEntries", "First run of has to be smaller than last run!");
        return kFALSE;
    }

    // create the query
    query.Form("UPDATE %s SET %s = '%s' "
               "WHERE run >= %d AND"
               "      run <= %d",
               TCConfig::kCalibMainTableName, name, value, first_run, last_run);

    // read from database
    Bool_t res = SendExec(query.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("ChangeRunEntries", "Could not set the information '%s' for runs %d to %d to '%s'!",
                                                 name, first_run, last_run, value);
        return kFALSE;
    }

    return kTRUE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeSetEntry(const Char_t* data, const Char_t* calibration, Int_t set,
                                      const Char_t* name, const Char_t* value)
{
    // Change the information 'name' of the 'set'-th set of the calibration 'calibration'
    // for the calibration data 'data' to 'value'.

    TString query;
    Char_t table[256];

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return kFALSE;

    // get the data table
    if (!SearchTable(data, table))
    {
        if (!fSilence) Error("ChangeSetEntry", "No data table found!");
        return 0;
    }

    // get the first run of the set
    Int_t first_run = GetFirstRunOfSet(data, calibration, set);

    // create the query
    query.Form("UPDATE %s SET %s = '%s' "
               "WHERE calibration = '%s' AND "
               "first_run = %d",
               table, name, value, calibration, first_run);

    // read from database
    Bool_t res = SendExec(query.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("ChangeSetEntry", "Could not set the information '%s' for set %d of '%s' to '%s'!",
                             name, set, d->GetTitle(), value);
        return kFALSE;
    }

    return kTRUE;
}

//______________________________________________________________________________
Int_t TCMySQLManager::GetNsets(const Char_t* data, const Char_t* calibration)
{
    // Get the number of runsets for the calibration identifier 'calibration'
    // and the calibration data 'data'.

    TString query;
    Char_t table[256];

    // check for data
    if (!GetCalibData(data)) return 0;

    // get the data table
    if (!SearchTable(data, table))
    {
        if (!fSilence) Error("GetNsets", "No data table found!");
        return 0;
    }

    // create the query
    query.Form("SELECT first_run FROM %s WHERE "
               "calibration = '%s' "
               "ORDER BY first_run ASC",
               table, calibration);

    // read from database
    TSQLResult* res = SendQuery(query.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("GetNsets", "No runsets found in table '%s'!", table);
        return 0;
    }

    // count rows
    Int_t rows = 0;
    TSQLRow* r = res->Next();
    while (r)
    {
        rows++;
        delete r;
        r = res->Next();
    }
    delete res;

    return rows;
}

//______________________________________________________________________________
Int_t TCMySQLManager::GetFirstRunOfSet(const Char_t* data, const Char_t* calibration, Int_t set)
{
    // Get the first run of the runsets 'set' for the calibration identifier
    // 'calibration' and the calibration data 'data'.

    Char_t tmp[256];

    // check for data
    if (!GetCalibData(data)) return 0;

    // get the data
    if (SearchSetEntry(data, calibration, set, "first_run", tmp)) return atoi(tmp);
    else
    {
        if (!fSilence) Error("GetFirstRunOfSet", "Could not find first run of set!");
        return 0;
    }
}

//______________________________________________________________________________
Int_t TCMySQLManager::GetLastRunOfSet(const Char_t* data, const Char_t* calibration, Int_t set)
{
    // Get the last run of the runsets 'set' for the calibration identifier
    // 'calibration' and the calibration data 'data'.

    Char_t tmp[256];

    // check for data
    if (!GetCalibData(data)) return 0;

    // get the data
    if (SearchSetEntry(data, calibration, set, "last_run", tmp)) return atoi(tmp);
    else
    {
        if (!fSilence) Error("GetLastRunOfSet", "Could not find last run of set!");
        return 0;
    }
}

//______________________________________________________________________________
void TCMySQLManager::GetDescriptionOfSet(const Char_t* data, const Char_t* calibration,
                                         Int_t set, Char_t* outDesc)
{
    // Get the description of the runsets 'set' for the calibration identifier
    // 'calibration' and the calibration data 'data'.

    Char_t tmp[256];

    // check for data
    if (!GetCalibData(data)) return;

    // get the data
    if (SearchSetEntry(data, calibration, set, "description", tmp)) strcpy(outDesc, tmp);
    else
    {
        if (!fSilence) Error("GetDescriptionOfSet", "Could not find description of set!");
    }
}

//______________________________________________________________________________
void TCMySQLManager::GetChangeTimeOfSet(const Char_t* data, const Char_t* calibration,
                                        Int_t set, Char_t* outTime)
{
    // Get the change time of the runsets 'set' for the calibration identifier
    // 'calibration' and the calibration data 'data'.

    Char_t tmp[256];

    // check for data
    if (!GetCalibData(data)) return;

    // get the data
    if (SearchSetEntry(data, calibration, set, "changed", tmp)) strcpy(outTime, tmp);
    else
    {
        if (!fSilence) Error("GetChangeTimeOfSet", "Could not find change time of set!");
    }
}

//______________________________________________________________________________
Int_t* TCMySQLManager::GetRunsOfCalibration(const Char_t* calibration, Int_t* outNruns)
{
    // Returns a list of all runs of the calibration identifier 'calibration'.
    // If 'outNruns' is not zero the number of runs will be written to this variable.
    // NOTE: the run array must be destroyed by the caller.

    // default data
    const Char_t* default_data = "Data.Tagger.T0";

    // init result variables
    Int_t nruns = 0;
    Int_t* runs = 0;

    // get no. of sets
    Int_t nsets = GetNsets(default_data, calibration);

    // init helpers for sets
    Int_t* nruns_set = new Int_t[nsets];
    Int_t** runs_set = new Int_t*[nsets];

    // loop over sets
    for (Int_t i = 0; i < nsets; i++)
    {
        // get list of runs for set
        runs_set[i] = GetRunsOfSet(default_data, calibration, i, &nruns_set[i]);

        // update total number of runs
        nruns += nruns_set[i];
    }

    // create array of run numbers
    runs = new Int_t[nruns];

    // loop over sets
    Int_t index = 0;
    for (Int_t i = 0; i < nsets; i++)
    {
        // loop over runs of this set and fill runs
        for (Int_t j = 0; j < nruns_set[i]; j++)
            runs[index + j] = runs_set[i][j];

        // set new index
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

    // return number of runs
    if (outNruns) *outNruns = nruns;

    return runs;
}

//______________________________________________________________________________
Int_t* TCMySQLManager::GetRunsOfSet(const Char_t* data, const Char_t* calibration,
                                    Int_t set, Int_t* outNruns = 0)
{
    // Return the list of runs that are in the set 'set' of the calibration data 'data'
    // for the calibration identifier 'calibration'.
    // If 'outNruns' is not zero the number of runs will be written to this variable.
    // NOTE: the run array must be destroyed by the caller.

    TString query;
    // check for data
    if (!GetCalibData(data)) return 0;

    // get first and last run
    Int_t first_run = GetFirstRunOfSet(data, calibration, set);
    Int_t last_run = GetLastRunOfSet(data, calibration, set);

    // check first run
    if (!first_run)
    {
        if (!fSilence) Error("GetRunsOfSet", "Could not find runs of set %d!", set);
        return 0;
    }

    //
    // get all the runs that lie between first and last run
    //

    // create the query
    query.Form("SELECT run FROM %s "
               "WHERE run >= %d "
               "AND run <= %d "
               "ORDER by run,time",
               TCConfig::kCalibMainTableName, first_run, last_run);

    // read from database
    TSQLResult* res = SendQuery(query.Data());

    // create list for run numbers
    TList run_numbers;
    run_numbers.SetOwner(kTRUE);

    // read all rows/runs
    TSQLRow* r = res->Next();
    while (r)
    {
        run_numbers.Add(new TObjString(r->GetField(0)));
        delete r;
        r = res->Next();
    }

    // get number of runs
    Int_t nruns = run_numbers.GetSize();

    // create run array
    Int_t* runs = new Int_t[nruns];

    // read all runs
    for (Int_t i = 0; i < nruns; i++)
    {
        // get next run
        TObjString* rn = (TObjString*) run_numbers.At(i);

        // save run number
        runs[i] = atoi(rn->GetString().Data());
    }

    // clean-up
    delete res;

    // write number of runs
    if (outNruns) *outNruns = nruns;

    return runs;
}

//______________________________________________________________________________
Int_t TCMySQLManager::GetSetForRun(const Char_t* data, const Char_t* calibration, Int_t run)
{
    // Return the number of the set of the calibration data 'data' for the calibration
    // identifier 'calibration' the run 'run' belongs to.
    // Return -1 if there is no such set.

    // check for data
    if (!GetCalibData(data)) return -1;

    // get number of sets
    Int_t nSet = GetNsets(data, calibration);
    if (!nSet) return -1;

    // check if run exists
    TString tmp;
    if (!SearchRunEntry(run, "run", tmp))
    {
        if (!fSilence) Error("GetSetForRun", "Run has no valid run number!");
        return -1;
    }

    // loop over sets
    for (Int_t i = 0; i < nSet; i++)
    {
        // get first and last run
        Int_t first_run = GetFirstRunOfSet(data, calibration, i);
        Int_t last_run = GetLastRunOfSet(data, calibration, i);

        // check if run is in this set
        if (run >= first_run && run <= last_run) return i;
    }

    // no set found here
    return -1;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ReadParametersRun(const Char_t* data, const Char_t* calibration, Int_t run,
                                         Double_t* par, Int_t length)
{
    // Read 'length' parameters of the calibration data 'data' for the calibration identifier
    // 'calibration' valid for the run 'run' from the database to the value array 'par'.
    // Return kFALSE if an error occurred, otherwise kTRUE.

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return kFALSE;

    // get set
    Int_t set = GetSetForRun(data, calibration, run);

    // check set
    if (set == -1)
    {
        if (!fSilence) Error("ReadParametersRun", "No set of '%s' found for run %d",
                             d->GetTitle(), run);
        return kFALSE;
    }

    // call main parameter reading method
    return ReadParameters(data, calibration, set, par, length);
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ReadParameters(const Char_t* data, const Char_t* calibration, Int_t set,
                                      Double_t* par, Int_t length)
{
    // Read 'length' parameters of the 'set'-th set of the calibration data 'data'
    // for the calibration identifier 'calibration' from the database to the value array 'par'.
    // Return kFALSE if an error occurred, otherwise kTRUE.

    TString query;
    Char_t table[256];

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return kFALSE;

    // get the data table
    if (!SearchTable(data, table))
    {
        if (!fSilence) Error("ReadParameters", "No data table found!");
        return kFALSE;
    }

    // get the first run of the set
    Int_t first_run = GetFirstRunOfSet(data, calibration, set);

    // check first run
    if (!first_run)
    {
        if (!fSilence) Error("ReadParameters", "No calibration found for set %d of '%s'!",
                             set, d->GetTitle());
        return kFALSE;
    }

    // create the query
    query.Form("SELECT * FROM %s WHERE "
               "calibration = '%s' AND "
               "first_run = %d",
               table, calibration, first_run);

    // read from database
    TSQLResult* res = SendQuery(query.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("ReadParameters", "No calibration found for set %d of '%s'!",
                             set, d->GetTitle());
        return kFALSE;
    }

    // get data (parameters start at field 5)
    TSQLRow* row = res->Next();

    // check row
    if (!row)
    {
        if (!fSilence) Error("ReadParameters", "No calibration found for set %d of '%s'!",
                             set, d->GetTitle());
        delete res;
        return kFALSE;
    }

    for (Int_t i = 0; i < length; i++) par[i] = atof(row->GetField(i+5));

    // clean-up
    delete row;
    delete res;

    // user information
    if (!fSilence) Info("ReadParameters", "Read %d parameters of '%s' from the database",
                        length, d->GetTitle());

    return kTRUE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::WriteParameters(const Char_t* data, const Char_t* calibration, Int_t set,
                                       Double_t* par, Int_t length)
{
    // Write 'length' parameters of the 'set'-th set of the calibration data 'data'
    // for the calibration identifier 'calibration' from the value array 'par' to the database.
    // Return kFALSE if an error occurred, otherwise kTRUE.

    Char_t table[256];

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return kFALSE;

    // get the data table
    if (!SearchTable(data, table))
    {
        if (!fSilence) Error("WriteParameters", "No data table found!");
        return kFALSE;
    }

    // get the first run of the set
    Int_t first_run = GetFirstRunOfSet(data, calibration, set);

    // check first run
    if (!first_run)
    {
        if (!fSilence) Error("WriteParameters", "Could not write parameters of '%s'!",
                             d->GetTitle());
        return kFALSE;
    }

    // prepare the insert query
    TString query = TString::Format("UPDATE %s SET ", table);

    // read all parameters and write them to new query
    for (Int_t j = 0; j < length; j++)
    {
        // append parameter to query
        query.Append(TString::Format("par_%03d = %.17g", j, par[j]));
        if (j != length - 1) query.Append(",");
    }

    // finish query
    query.Append(TString::Format(" WHERE calibration = '%s' AND first_run = %d",
                                 calibration, first_run));

    // write data to database
    Bool_t res = SendExec(query.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("WriteParameters", "Could not write parameters of '%s'!",
                             d->GetTitle());
        return kFALSE;
    }
    else
    {
        if (!fSilence) Info("WriteParameters", "Wrote %d parameters of '%s' to the database",
                            length, d->GetTitle());
        return kTRUE;
    }
}

//______________________________________________________________________________
Bool_t TCMySQLManager::InitDatabase(Bool_t interact)
{
    // Init a new CaLib database on a MySQL server in interactive mode when
    // 'interact' is kTRUE.
    // Return kTRUE on success, otherwise kFALSE.

    // check server connection
    if (!IsConnected())
    {
        if (!fSilence) Error("InitDatabase", "No connection to the database!");
        return kFALSE;
    }

    // ask for user confirmation
    if (interact)
    {
        Char_t answer[256];
        if (fDBType == kSQLite)
        {
            printf("\nWARNING: You are about to initialize a new CaLib database.\n"
                   "         All existing tables in the database '%s'\n"
                   "         will be deleted!\n\n", fDB->GetDB());
        }
        else
        {
            printf("\nWARNING: You are about to initialize a new CaLib database.\n"
                   "         All existing tables in the database '%s' on '%s'\n"
                   "         will be deleted!\n\n", fDB->GetDB(), fDB->GetHost());
        }
        printf("Are you sure to continue? (yes/no) : ");
        Int_t ret = scanf("%s", answer);
        if (strcmp(answer, "yes"))
        {
            printf("Aborted.\n");
            return kFALSE;
        }
    }

    // create the main table
    CreateMainTable();

    // create the data tables
    TIter next(fData);
    TCCalibData* d;
    Bool_t err = kFALSE;
    while ((d = (TCCalibData*)next()))
    {
        // create the data table
        if (!CreateDataTable(d->GetName(), d->GetSize()))
            err = kTRUE;
    }

    // add MC run
    AddRunMC();

    // check for errors
    if (err)
    {
        if (!fSilence) Error("InitDatabase", "An error occurred during database initialization!");
        return kFALSE;
    }

    return kTRUE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::UpgradeDatabase(Int_t version)
{
    // Upgrade the CaLib database to version 'version'.
    // Return kTRUE on success, otherwise kFALSE.

    // check server connection
    if (!IsConnected())
    {
        if (!fSilence) Error("UpgradeDatabase", "No connection to the database!");
        return kFALSE;
    }

    // ask for user confirmation
    Char_t answer[256];
    if (fDBType == kSQLite)
    {
        printf("\nWARNING: You are about to update your existing CaLib database '%s'\n"
               "         It is highly recommended to create a complete backup\n"
               "         before proceeding!\n\n", fDB->GetDB());
    }
    else
    {
        printf("\nWARNING: You are about to update your existing CaLib database '%s' on '%s'\n"
               "         It is highly recommended to create a complete backup using\n"
               "         e.g. mysqldump before proceeding!\n\n", fDB->GetDB(), fDB->GetHost());
    }
    printf("Are you sure to continue? (yes/no) : ");
    Int_t ret = scanf("%s", answer);
    if (strcmp(answer, "yes"))
    {
        printf("Aborted.\n");
        return kFALSE;
    }

    // create queries (update only this part in the future)
    Int_t nQuery = 0;
    Char_t query[3][256];
    switch (version)
    {
        // version 3:
        // - add two columns in run_main for bad scaler reads
        // - remove livetime table
        case 3:
        {
            nQuery = 3;
            strcpy(query[0], "ALTER TABLE run_main ADD scr_n INT DEFAULT -1 AFTER size");
            strcpy(query[1], "ALTER TABLE run_main ADD scr_bad VARCHAR(512) AFTER scr_n");
            strcpy(query[2], "DROP TABLE IF EXISTS livetime");
            break;
        }
        // version 4:
        // - change type of table scr_bad to TEXT
        // - add Data.Tagger.Pol data table
        case 4:
        {
            nQuery = 1;
            strcpy(query[0], "ALTER TABLE run_main MODIFY scr_bad TEXT");

            // add Data.Tagger.Pol data table
            if (!AddNewDataTable("Data.Tagger.Pol"))
                Error("UpgradeDatabase", "Some errors occurred while adding data table for 'Data.Tagger.Pol'!");

            // add pizza detector tables
            if (!AddNewDataTable("Data.Pizza.Phi"))
                Error("UpgradeDatabase", "Some errors occurred while adding data table for 'Data.Pizza.Phi'!");
            if (!AddNewDataTable("Data.Pizza.T0"))
                Error("UpgradeDatabase", "Some errors occurred while adding data table for 'Data.Pizza.T0'!");
            if (!AddNewDataTable("Data.Pizza.E0"))
                Error("UpgradeDatabase", "Some errors occurred while adding data table for 'Data.Pizza.E0'!");
            if (!AddNewDataTable("Data.Pizza.E1"))
                Error("UpgradeDatabase", "Some errors occurred while adding data table for 'Data.Pizza.E1'!");

            break;
        }
        default:
        {
            Error("UpgradeDatabase", "Database upgrade to version %d not implemented!", version);
            return kFALSE;
        }
    }

    // loop over queries
    Int_t queryOk = 0;
    for (Int_t i = 0; i < nQuery; i++)
    {
        // send query
        Bool_t res = SendExec(query[i]);

        // check result
        if (!res) Error("UpgradeDatabase", "Could not execute query %d!", i+1);
        else
        {
            Info("UpgradeDatabase", "Executed query %d", i+1);
            queryOk++;
        }
    }

    // check final result
    if (queryOk == nQuery)
    {
        Info("UpgradeDatabase", "Performed upgrade of database");
        return kTRUE;
    }
    else
    {
        Error("UpgradeDatabase", "Some errors occurred during database upgrade - please check the consistency of your database!");
        return kFALSE;
    }
}

//______________________________________________________________________________
Bool_t TCMySQLManager::AddNewDataTable(const Char_t* data)
{
    // Add a new table for the calibration data 'data' to the database including
    // a default calibration set for all calibrations found in the database.

    Char_t tmp[256];

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return kFALSE;

    // check if table exists already
    if (SendExec(TString::Format("SELECT 1 from %s LIMIT 1", d->GetTableName())))
    {
        if (!fSilence) Error("AddNewDataTable", "Data table for '%s' exists already!", data);
        return kFALSE;
    }

    // add data table
    if (CreateDataTable(d->GetName(), d->GetSize()))
    {
        if (!fSilence) Info("AddNewDataTable", "Added data table for '%s'", data);
    }
    else
    {
        if (!fSilence) Error("AddNewDataTable", "Could not add data table for '%s'!", data);
        return kFALSE;
    }

    //
    // create default set for all calibrations
    //

    // get all calibrations
    TList* list = GetAllCalibrations();
    if (!list)
    {
        if (!fSilence) Warning("AddNewDataTable", "No default set added because no calibration was found!");
        return kTRUE;
    }

    // loop over all calibrations
    Bool_t ret = kTRUE;
    TIter next(list);
    TObjString* s;
    while ((s = (TObjString*)next()))
    {
        // get set configuration from the Data.Tagger.T0 data table
        GetDescriptionOfSet("Data.Tagger.T0", s->GetString().Data(), 0, tmp);
        Int_t first_run = GetFirstRunOfSet("Data.Tagger.T0", s->GetString().Data(), 0);
        Int_t last_run = GetLastRunOfSet("Data.Tagger.T0", s->GetString().Data(), GetNsets("Data.Tagger.T0", s->GetString().Data())-1);
        if (!fSilence) Info("AddNewDataTable", "Adding set [%d,%d] for calibration '%s'", first_run, last_run, s->GetString().Data());

        // add set
        if (!AddDataSet(data, s->GetString().Data(), tmp, first_run, last_run, 0.0))
        {
            if (!fSilence) Error("AddNewDataTable", "Could not add set to data table '%s'!", data);
            ret = kFALSE;
        }
    }

    // clean-up
    delete list;

    return ret;
}

//______________________________________________________________________________
void TCMySQLManager::AddRunFiles(const Char_t* path, const Char_t* target,
                                 const Char_t* runPrefix)
{
    // Look for raw ACQU files in 'path' and add all runs with the prefix 'runPrefix'
    // to the database using the target specifier 'target'.

    struct tm tm;
    Char_t time[256];

    // read the raw files
    TCReadACQU r(path, runPrefix);
    Int_t nRun = r.GetNFiles();

    // ask for user confirmation
    Char_t answer[256];
    if (fDBType == kSQLite)
    {
        printf("\n%d runs were found in '%s'\n"
               "They will be added to the database '%s'\n",
               nRun, path, fDB->GetDB());
    }
    else
    {
        printf("\n%d runs were found in '%s'\n"
               "They will be added to the database '%s' on '%s'\n",
               nRun, path, fDB->GetDB(), fDB->GetHost());
    }
    printf("Are you sure to continue? (yes/no) : ");
    Int_t ret = scanf("%s", answer);
    if (strcmp(answer, "yes"))
    {
        printf("Aborted.\n");
        return;
    }

    // loop over runs
    Int_t nRunAdded = 0;
    for (Int_t i = 0; i < nRun; i++)
    {
        TCACQUFile* f = r.GetFile(i);

        // convert the time string
        strptime(f->GetTime(), "%a %b %d %H:%M:%S %Y", &tm);
        strftime(time, sizeof(time), "%Y-%m-%d %H:%M:%S", &tm);

        // prepare the insert query
        TString ins_query = TString::Format("INSERT INTO %s (run, path, filename, time, description, run_note, size, target) "
                                            "VALUES ( "
                                            "%d, "
                                            "\"%s\", "
                                            "\"%s\", "
                                            "\"%s\", "
                                            "\"%s\", "
                                            "\"%s\", "
                                            "%lld, "
                                            "\"%s\" )",
                                            TCConfig::kCalibMainTableName,
                                            f->GetRun(),
                                            path,
                                            f->GetFileName(),
                                            time,
                                            f->GetDescription(),
                                            f->GetRunNote(),
                                            f->GetSize(),
                                            target);

        // try to write data to database
        Bool_t res = SendExec(ins_query.Data());
        if (!res)
        {
            Warning("AddRunFiles", "Run %d of file '%s/%s' could not be added to the database!",
                    f->GetRun(), path, f->GetFileName());
        }
        else
        {
            nRunAdded++;
        }
    }

    // user information
    if (!fSilence) Info("AddRunFiles", "Added %d runs to the database", nRunAdded);
}

//______________________________________________________________________________
void TCMySQLManager::AddRun(Int_t run, const Char_t* target, const Char_t* desc)
{
    // Add a run to the run database using the run number 'run', the target
    // identifier 'target' and the description 'desc'.

    // string arrays
    Char_t t[256] = "";
    Char_t d[256] = "";
    if (target) strcpy(t, target);
    if (desc) strcpy(d, desc);

    // prepare the insert query
    TString ins_query = TString::Format("INSERT INTO %s (run, description, target) "
                                        "VALUES ( "
                                        "%d, "
                                        "'%s', "
                                        "'%s' )",
                                        TCConfig::kCalibMainTableName,
                                        run,
                                        d,
                                        t);

    // try to write data to database
    Bool_t res = SendExec(ins_query.Data());
    if (!res)
    {
        if (!fSilence) Warning("AddRun", "Run %d could not be added to the database!", run);
    }
    else
    {
        if (!fSilence) Info("AddRun", "Run %d ('%s') for target '%s' was added to the database",
                            run, d, t);
    }
}

//______________________________________________________________________________
void TCMySQLManager::AddRunMC(Int_t run, const Char_t* target, const Char_t* desc)
{
    // Add an MC run to the run database using optionally the run number 'run', the target
    // identifier 'target' and the description 'desc'.

    // add MC run
    AddRun(run, target, desc);
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeRunPath(Int_t first_run, Int_t last_run, const Char_t* path)
{
    // Change the path of all runs between 'first_run' and 'last_run' to 'path'.

    // change path
    return ChangeRunEntries(first_run, last_run, "path", path);
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeRunTarget(Int_t first_run, Int_t last_run, const Char_t* target)
{
    // Change the target of all runs between 'first_run' and 'last_run' to 'target'.

    // change target
    return ChangeRunEntries(first_run, last_run, "target", target);
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeRunTargetPol(Int_t first_run, Int_t last_run, const Char_t* target_pol)
{
    // Change the target polarization of all runs between 'first_run' and
    // 'last_run' to 'target_pol'.

    // change target polarization
    return ChangeRunEntries(first_run, last_run, "target_pol", target_pol);
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeRunTargetPolDeg(Int_t first_run, Int_t last_run, Double_t target_pol_deg)
{
    // Change the degree of target polarization of all runs between 'first_run' and
    // 'last_run' to 'target_pol_deg'.

    // create string
    Char_t tmp[256];
    sprintf(tmp, "%lf", target_pol_deg);

    // change degree of target polarization
    return ChangeRunEntries(first_run, last_run, "target_pol_deg", tmp);
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeRunBeamPol(Int_t first_run, Int_t last_run, const Char_t* beam_pol)
{
    // Change the beam polarization of all runs between 'first_run' and
    // 'last_run' to 'beam_pol'.

    // change beam polarization
    return ChangeRunEntries(first_run, last_run, "beam_pol", beam_pol);
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeRunBeamPolDeg(Int_t first_run, Int_t last_run, Double_t beam_pol_deg)
{
    // Change the degree of beam polarization of all runs between 'first_run' and
    // 'last_run' to 'beam_pol_deg'.

    // create string
    Char_t tmp[256];
    sprintf(tmp, "%lf", beam_pol_deg);

    // change degree of beam polarization
    return ChangeRunEntries(first_run, last_run, "beam_pol_deg", tmp);
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeRunNScR(const Int_t run, const Int_t nscr)
{
    // Change the number of scaler reads for run 'run' to 'nscr'.  Returns kTRUE
    // on success, kFALSE otherwise.

    // create string
    Char_t tmp[16];
    sprintf(tmp, "%d", nscr);

    // change number of scaler reads
    return ChangeRunEntries(run, run, "scr_n", tmp);
}

//______________________________________________________________________________
Int_t TCMySQLManager::GetRunNScR(const Int_t run)
{
   // Reads the number of scaler reads for run 'run' from the database. Returns
   // the read number of scaler read on success, or -2;

   // create string
   TString tmp;

   // get number of scaler reads
   if (!SearchRunEntry(run, "scr_n", tmp))
   {
       //
       return -2;
   }
   else
   {
       //
       return atoi(tmp.Data());
   }
}

/*
//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeRunBadScR(const Int_t run, const Int_t nbadscr, const Int_t* const badscr)
{
    // Change the bad scaler reads list string for run 'run' to a string of
    // comma and colon separated (ASCII) integers describing the 'nbadscr' bad
    // scaler reads in 'badscr'. Returns kTRUE on success, kFALSE otherwise.
    //
    // e.g.: badscr = {1,3,4,5,7,8,10} --> s = "1,3:5,7,8,10"
    //
    // With a 512 byte string and number of scaler reads <= 1000:
    //     Supported number of bad scaler reads: 128 (guaranteed) or more
    //     Supported number of scaler reads:     218 (guaranteed) or more

    // init query string
    TString s = "";

    // loop over input bad scaler reads list
    for (Int_t i = 0; i < nbadscr; i++)
    {
        Char_t tmp[16];
        Int_t j = 0;

        // loop over subsequent values in order to detect series, i.e., badscr[k]+1 == badscr[k+1]
        for (j = 0; i + j < nbadscr - 1; j++)
            if (badscr[i+j] + 1 != badscr[i+j+1]) break;

        // check for series with more than 2 elements
        if (j > 1)
        {
            // separate first and last value of series with ':'
            sprintf(tmp, "%d:%d,", badscr[i], badscr[i+j]);
            i += j;
        }
        else
        {
            // separate with ','
            sprintf(tmp, "%d,", badscr[i]);
        }

        // append to query string
        s.Append(tmp);
    }

    // remove tailing comma
    s.Chop();

    // check length
    if (s.Length() >= 512)
    {
        if (!fSilence) Error("ChangeRunBadScR", "Too many bad scaler reads!");
        return kFALSE;
    }

    // change list of bad scaler reads
    return ChangeRunEntries(run, run, "scr_bad", s.Data());
}
*/

//______________________________________________________________________________
Bool_t TCMySQLManager::ContainsCalibration(const Char_t* calibration)
{
    // Check if the calibration 'calibration' exists in the database.

    // loop over calibration data
    TIter next(fData);
    TCCalibData* d;
    while ((d = (TCCalibData*)next()))
    {
        // get list of calibrations
        TList* list = GetAllCalibrations(d->GetName());

        // no calibrations
        if (!list) continue;

        // check if calibration exists
        TIter next(list);
        Bool_t found = kFALSE;
        TObjString* s;
        while ((s = (TObjString*)next()))
        {
            if (!strcmp(calibration, s->GetString().Data()))
            {
                found = kTRUE;
                break;
            }
        }

        // clean-up
        delete list;

        // check if calibration was found
        if (found) return kTRUE;
    }

    return kFALSE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeCalibrationName(const Char_t* calibration, const Char_t* newCalibration)
{
    // Change the calibration identifer 'calibration' in all calibration sets
    // to 'newCalibration'.

    TString query;

    // check if calibration was not found
    if (!ContainsCalibration(calibration))
    {
        if (!fSilence) Error("ChangeCalibrationName", "Calibration '%s' was not found in database!",
                             calibration);
        return kFALSE;
    }

    // loop over calibration data
    TIter next(fData);
    TCCalibData* d;
    Bool_t succ = kTRUE;
    while ((d = (TCCalibData*)next()))
    {
        // create the query
        query.Form("UPDATE %s SET calibration = '%s' "
                   "WHERE calibration = '%s'",
                   d->GetTableName(), newCalibration, calibration);

        // read from database
        Bool_t res = SendExec(query.Data());

        // check result
        if (!res)
        {
            if (!fSilence) Error("ChangeCalibrationName", "Could not rename table '%s' of calibration '%s'!",
                                 d->GetTableName(), calibration);
            succ = kFALSE;
        }
    }

    if (!fSilence)
    {
        if (succ)
            Info("ChangeCalibrationName", "Renamed calibration '%s' to '%s'", calibration, newCalibration);
        else
            Error("ChangeCalibrationName", "Some problem(s) occurred - check your database for inconsistencies!");
    }

    return succ;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeCalibrationDescription(const Char_t* calibration, const Char_t* newDesc)
{
    // Change the calibration description of the 'calibration' in all calibration sets
    // to 'newDesc'.

    TString query;

    // check if calibration was not found
    if (!ContainsCalibration(calibration))
    {
        if (!fSilence) Error("ChangeCalibrationDescription", "Calibration '%s' was not found in database!",
                             calibration);
        return kFALSE;
    }

    // loop over calibration data
    TIter next(fData);
    TCCalibData* d;
    Bool_t succ = kTRUE;
    while ((d = (TCCalibData*)next()))
    {
        // create the query
        query.Form("UPDATE %s SET description = '%s' "
                   "WHERE calibration = '%s'",
                   d->GetTableName(), newDesc, calibration);

        // read from database
        Bool_t res = SendExec(query.Data());

        // check result
        if (!res)
        {
            if (!fSilence) Error("ChangeCalibrationDescription", "Could not change description in table '%s' of calibration '%s'!",
                                 d->GetTableName(), calibration);
            succ = kFALSE;
        }
    }

    if (!fSilence)
    {
        if (succ)
            Info("ChangeCalibrationDescription", "Changed description of calibration '%s' to '%s'", calibration, newDesc);
        else
            Error("ChangeCalibrationDescription", "Some problem(s) occurred - check your database for inconsistencies!");
    }

    return succ;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeCalibrationRunRange(const Char_t* calibration, const UInt_t firstRun,
                                                 const UInt_t lastRun)
{
    // Change the run range in all sets of the calibration 'calibration' to start from 'firstRun'
    // and to end at 'lastRun'. The runs must be present in the run database.
    // If 'firstRun'/'lastRun' is zero, the current start/stop run are kept.

    TString query;
    TString tmp;

    // check if calibration was not found
    if (!ContainsCalibration(calibration))
    {
        if (!fSilence) Error("ChangeCalibrationRunRange", "Calibration '%s' was not found in database!",
                             calibration);
        return kFALSE;
    }

    // quit if both firstRun and lastRun are zero
    if (!firstRun && !lastRun)
    {
        if (!fSilence) Error("ChangeCalibrationRunRange", "First and/or last run has to be non-zero!");
        return kFALSE;
    }

    // check if first run exists
    if (firstRun && !SearchRunEntry(firstRun, "run", tmp))
    {
        if (!fSilence) Error("ChangeCalibrationRunRange", "First run %d has no valid run number!",
                             firstRun);
        return kFALSE;
    }

    // check if last run exists
    if (lastRun && !SearchRunEntry(lastRun, "run", tmp))
    {
        if (!fSilence) Error("ChangeCalibrationRunRange", "Last run %d has no valid run number!",
                             lastRun);
        return kFALSE;
    }

    // check if first runs is smaller than last run
    if (firstRun && lastRun && firstRun > lastRun)
    {
        if (!fSilence) Error("ChangeCalibrationRunRange", "First run %d has to be smaller than last run %d",
                             firstRun, lastRun);
        return kFALSE;
    }

    // loop over calibration data
    TIter next(fData);
    TCCalibData* d;
    Bool_t succ = kTRUE;
    while ((d = (TCCalibData*)next()))
    {
        //
        // first run
        //

        if (firstRun)
        {
            // get the first run of the first set
            Int_t oldFirstRun = GetFirstRunOfSet(d->GetName(), calibration, 0);

            // execute the query
            query.Form("UPDATE %s SET first_run = %d WHERE calibration = '%s' and first_run = %d",
                       d->GetTableName(), firstRun, calibration, oldFirstRun);
            Bool_t res = SendExec(query.Data());

            // check result
            if (!res)
            {
                if (!fSilence) Error("ChangeCalibrationRunRange", "Could not change first run in table '%s' of calibration '%s'!",
                                     d->GetTableName(), calibration);
                succ = kFALSE;
            }
        }

        //
        // last run
        //

        if (lastRun)
        {
           // get the last run of the last set
            Int_t oldLastRun = GetLastRunOfSet(d->GetName(), calibration, GetNsets(d->GetName(), calibration)-1);

            // execute the query
            query.Form("UPDATE %s SET last_run = %d WHERE calibration = '%s' and last_run = %d",
                       d->GetTableName(), lastRun, calibration, oldLastRun);
            Bool_t res = SendExec(query.Data());

            // check result
            if (!res)
            {
                if (!fSilence) Error("ChangeCalibrationRunRange", "Could not change last run in table '%s' of calibration '%s'!",
                                     d->GetTableName(), calibration);
                succ = kFALSE;
            }
        }
    }

    if (!fSilence)
    {
        if (succ)
            Info("ChangeCalibrationRunRange", "Changed run range of calibration '%s' to [%d,%d]", calibration, firstRun, lastRun);
        else
            Error("ChangeCalibrationRunRange", "Some problem(s) occurred - check your database for inconsistencies!");
    }

    return succ;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::RemoveCalibration(const Char_t* calibration, const Char_t* data)
{
    // Remove the calibration data 'data' of the calibration 'calibration'.

    TString query;

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return kFALSE;

    // check if calibration was not found
    if (!ContainsCalibration(calibration))
    {
        if (!fSilence) Error("RemoveCalibration", "Calibration '%s' of calibration '%s' was not found in database!",
                             d->GetTitle(), calibration);
        return kFALSE;
    }

    // create the query
    query.Form("DELETE from %s "
               "WHERE calibration = '%s'",
               d->GetTableName(), calibration);

    // read from database
    Bool_t res = SendExec(query.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("RemoveCalibration", "Could not remove calibration '%s' of calibration '%s'!",
                             d->GetTitle(), calibration);
        return kFALSE;
    }

    if (!fSilence) Info("RemoveCalibration", "Removed calibration '%s' of calibration '%s'",
                        d->GetTitle(), calibration);

    return kTRUE;
}

//______________________________________________________________________________
Int_t TCMySQLManager::RemoveAllCalibrations(const Char_t* calibration)
{
    // Remove all calibrations with the calibration identifer 'calibration'.

    Int_t nCalib = 0;

    // loop over calibration data
    TIter next(fData);
    TCCalibData* d;
    while ((d = (TCCalibData*)next()))
    {
        // remove calibration
        if (RemoveCalibration(calibration, d->GetName())) nCalib++;
    }

    return nCalib;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::RemoveAllRuns()
{
    // Remove all runs from the database.

    TString query;

    // create the query
    query.Form("DELETE FROM %s", TCConfig::kCalibMainTableName);

    // read from database
    Bool_t res = SendExec(query.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("RemoveAllRuns", "Could not remove all runs!");
        return kFALSE;
    }
    else
    {
        if (!fSilence) Info("RemoveAllRuns", "Removed all runs");
        return kTRUE;
    }
}

//______________________________________________________________________________
void TCMySQLManager::AddCalibAR(CalibDetector_t det, const Char_t* calibFileAR,
                                const Char_t* calib, const Char_t* desc,
                                Int_t first_run, Int_t last_run)
{
    // Read the calibration for the detector 'det' from the AcquRoot calibration
    // file 'calibFileAR' and create calibration sets for the runs 'first_run'
    // to 'last_run' using the calibration name 'calib' and the description
    // 'desc'.

    // read the calibration file
    TCReadARCalib r(calibFileAR, kFALSE);

    // get the number of elements
    Int_t nDet = r.GetNelements();

    // user information
    if (nDet)
    {
        if (!fSilence) Info("AddCalibAR", "Found calibrations for %d detector elements", nDet);
    }
    else
    {
        if (!fSilence) Error("AddCalibAR", "No detector elements found in calibration file!");
        return;
    }

    // create generic parameter arrays
    Double_t eL[nDet];
    Double_t e0[nDet];
    Double_t e1[nDet];
    Double_t t0[nDet];
    Double_t t1[nDet];

    // read generic parameters
    for (Int_t i = 0; i < nDet; i++)
    {
        eL[i] = r.GetElement(i)->GetEnergyLow();
        e0[i] = r.GetElement(i)->GetPedestal();
        e1[i] = r.GetElement(i)->GetADCGain();
        t0[i] = r.GetElement(i)->GetOffset();
        t1[i] = r.GetElement(i)->GetTDCGain();
    }

    // read detector specific calibration values
    // and write the data to the database (depends also
    // on the detector)
    switch (det)
    {
        // tagger
        case kDETECTOR_TAGG:
        {
            // write to database
            AddDataSet("Data.Tagger.T0", calib, desc, first_run, last_run, t0, nDet);

            break;
        }
        // CB
        case kDETECTOR_CB:
        {
            // write to database
            AddDataSet("Data.CB.E1", calib, desc, first_run, last_run, e1, nDet);
            AddDataSet("Data.CB.T0", calib, desc, first_run, last_run, t0, nDet);

            break;
        }
        // TAPS
        case kDETECTOR_TAPS:
        {
            // read the calibration file again for the SG stuff
            TCReadARCalib rSG(calibFileAR, kFALSE, "TAPSSG:");

            // get the number of elements
            Int_t nDetSG = rSG.GetNelements();

            // user information
            if (!nDetSG)
            {
                if (!fSilence) Error("AddCalibAR", "No TAPS SG detector elements found in calibration file!");
                return;
            }

            // create SG parameter arrays
            Double_t e0SG[nDetSG];
            Double_t e1SG[nDetSG];

            // read SG parameters
            for (Int_t i = 0; i < nDetSG; i++)
            {
                e0SG[i] = rSG.GetElement(i)->GetPedestal();
                e1SG[i] = rSG.GetElement(i)->GetADCGain();
            }

            // write to database
            AddDataSet("Data.TAPS.CFD", calib, desc, first_run, last_run, eL, nDet);
            AddDataSet("Data.TAPS.LG.E0", calib, desc, first_run, last_run, e0, nDet);
            AddDataSet("Data.TAPS.LG.E1", calib, desc, first_run, last_run, e1, nDet);
            AddDataSet("Data.TAPS.SG.E0", calib, desc, first_run, last_run, e0SG, nDetSG);
            AddDataSet("Data.TAPS.SG.E1", calib, desc, first_run, last_run, e1SG, nDetSG);
            AddDataSet("Data.TAPS.T0", calib, desc, first_run, last_run, t0, nDet);
            AddDataSet("Data.TAPS.T1", calib, desc, first_run, last_run, t1, nDet);

            break;
        }
        // PID
        case kDETECTOR_PID:
        {
            // create special parameter arrays
            Double_t phi[nDet];

            // read special parameters
            for (Int_t i = 0; i < nDet; i++)
            {
                phi[i] = r.GetElement(i)->GetZ();
            }

            // write to database
            AddDataSet("Data.PID.Phi", calib, desc, first_run, last_run, phi, nDet);
            AddDataSet("Data.PID.E0", calib, desc, first_run, last_run, e0, nDet);
            AddDataSet("Data.PID.E1", calib, desc, first_run, last_run, e1, nDet);
            AddDataSet("Data.PID.T0", calib, desc, first_run, last_run, t0, nDet);

            break;
        }
        // Veto
        case kDETECTOR_VETO:
        {
            // write to database
            AddDataSet("Data.Veto.LED", calib, desc, first_run, last_run, eL, nDet);
            AddDataSet("Data.Veto.E0", calib, desc, first_run, last_run, e0, nDet);
            AddDataSet("Data.Veto.E1", calib, desc, first_run, last_run, e1, nDet);
            AddDataSet("Data.Veto.T0", calib, desc, first_run, last_run, t0, nDet);
            AddDataSet("Data.Veto.T1", calib, desc, first_run, last_run, t1, nDet);

            break;
        }
        // Pizza
        case kDETECTOR_PIZZA:
        {
            // create special parameter arrays
            Double_t phi[nDet];

            // read special parameters
            for (Int_t i = 0; i < nDet; i++)
            {
                phi[i] = r.GetElement(i)->GetZ();
            }

            // write to database
            AddDataSet("Data.Pizza.Phi", calib, desc, first_run, last_run, phi, nDet);
            AddDataSet("Data.Pizza.E0", calib, desc, first_run, last_run, e0, nDet);
            AddDataSet("Data.Pizza.E1", calib, desc, first_run, last_run, e1, nDet);
            AddDataSet("Data.Pizza.T0", calib, desc, first_run, last_run, t0, nDet);

            break;
        }
        // no detector
        case kDETECTOR_NODET:
        {
            // do nothing
            break;
        }
    }
}

//______________________________________________________________________________
void TCMySQLManager::CreateMainTable()
{
    // Create the main table for CaLib.

    // user information
    if (!fSilence) Info("CreateMainTable", "Creating main CaLib table");

    // delete the old table if it exists
    SendExec(TString::Format("DROP TABLE IF EXISTS %s", TCConfig::kCalibMainTableName).Data());

    // create the table
    SendExec(TString::Format("CREATE TABLE %s ( %s )",
                             TCConfig::kCalibMainTableName, TCConfig::kCalibMainTableFormat).Data());

    // add timestamp update mechanism
    if (fDBType == kMySQL)
    {
        SendExec(TString::Format("ALTER TABLE %s DROP changed", TCConfig::kCalibMainTableName).Data());
        SendExec(TString::Format("ALTER TABLE %s ADD changed TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP "
                                 "ON UPDATE CURRENT_TIMESTAMP", TCConfig::kCalibMainTableName).Data());
    }
    else if (fDBType == kSQLite)
    {
        // delete the old timestamp update trigger if it exists
        SendExec(TString::Format("DROP TRIGGER IF EXISTS timestamp_update_%s", TCConfig::kCalibMainTableName).Data());

        // create the timestamp update trigger
        SendExec(TString::Format("CREATE TRIGGER after_%s_update "
                                 "AFTER UPDATE "
                                 "ON %s "
                                 "FOR EACH ROW "
                                 "WHEN NEW.changed <= OLD.changed "
                                 "BEGIN "
                                 "UPDATE %s SET changed = CURRENT_TIMESTAMP WHERE run = OLD.run; "
                                 "END",
                                 TCConfig::kCalibMainTableName, TCConfig::kCalibMainTableName, TCConfig::kCalibMainTableName).Data());
    }
}

//______________________________________________________________________________
Bool_t TCMySQLManager::CreateDataTable(const Char_t* data, Int_t nElem)
{
    // Create the table for the calibration data 'data' for 'nElem' elements.

    Char_t table[256];

    // check for data
    if (!GetCalibData(data)) return kFALSE;

    // get the table name
    if (!SearchTable(data, table))
    {
        if (!fSilence) Error("CreateDataTable", "No data table found!");
        return kFALSE;
    }

    if (!fSilence) Info("CreateDataTable", "Adding data table '%s' for %d elements", table, nElem);

    // delete the old table if it exists
    SendExec(TString::Format("DROP TABLE IF EXISTS %s", table));

    // prepare CREATE TABLE query
    TString query;
    query.Append(TString::Format("CREATE TABLE %s ( %s ", table, TCConfig::kCalibDataTableHeader));

    // loop over elements
    for (Int_t j = 0; j < nElem; j++)
    {
        query.Append(TString::Format("par_%03d DOUBLE DEFAULT 0", j));
        if (j != nElem - 1) query.Append(", ");
    }

    // finish preparing the query
    query.Append(TCConfig::kCalibDataTableSettings);
    query.Append(" )");

    // submit the query
    if (!SendExec(query.Data()))
    {
        if (!fSilence) Error("CreateDataTable", "An error occurred during data table creation for '%s'!", data);
        return kFALSE;
    }

    // add timestamp update mechanism
    if (fDBType == kMySQL)
    {
        SendExec(TString::Format("ALTER TABLE %s DROP changed", table).Data());
        SendExec(TString::Format("ALTER TABLE %s ADD changed TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP "
                                 "ON UPDATE CURRENT_TIMESTAMP AFTER last_run", table).Data());
    }
    else if (fDBType == kSQLite)
    {
        // delete the old timestamp update trigger if it exists
        SendExec(TString::Format("DROP TRIGGER IF EXISTS timestamp_update_%s", table).Data());

        // create the timestamp update trigger
        SendExec(TString::Format("CREATE TRIGGER after_%s_update "
                                 "AFTER UPDATE "
                                 "ON %s "
                                 "FOR EACH ROW "
                                 "WHEN NEW.changed <= OLD.changed "
                                 "BEGIN "
                                 "UPDATE %s SET changed = CURRENT_TIMESTAMP WHERE calibration = OLD.calibration AND first_run = old.first_run ; "
                                 "END",
                                 table, table, table).Data());
    }

    return kTRUE;
}

//______________________________________________________________________________
TList* TCMySQLManager::SearchDistinctEntries(const Char_t* field, const Char_t* table)
{
    // Return a list of TStrings containing all distinct entries of the type
    // 'field' in the table 'table' of the CaLib database.
    // If nothing is found 0 is returned.
    // NOTE: The list must be destroyed by the caller.

    TString query;

    // get all entries
    query.Form("SELECT DISTINCT %s from %s", field, table);
    TSQLResult* res = SendQuery(query.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("SearchDistinctEntries", "Could find field '%s' in table '%s'!",
                             field, table);
        return 0;
    }

    // create list for rows
    TList rows;
    rows.SetOwner(kTRUE);

    // read all rows
    TSQLRow* r = res->Next();
    while (r)
    {
        rows.Add(new TObjString(r->GetField(0)));
        delete r;
        r = res->Next();
    }

    // get number of entries
    Int_t n = rows.GetSize();

    // count rows
    if (!n)
    {
        delete res;
        return 0;
    }

    // create list
    TList* list = new TList();
    list->SetOwner(kTRUE);

    // read all entries and add them to the list
    for (Int_t i = 0; i < n; i++)
    {
        // get next entry
        TObjString* entry = (TObjString*) rows.At(i);

        // save entry
        list->Add(new TObjString(*entry));
    }

    // clean-up
    delete res;

    // return the list
    return list;
}

//______________________________________________________________________________
TList* TCMySQLManager::GetAllTargets()
{
    // Return a list of TStrings containing all targets in the database.
    // If no targets were found 0 is returned.
    // NOTE: The list must be destroyed by the caller.

    return SearchDistinctEntries("target", TCConfig::kCalibMainTableName);
}

//______________________________________________________________________________
TList* TCMySQLManager::GetAllCalibrations(const Char_t* data)
{
    // Return a list of TStrings containing all calibration identifiers in the database
    // for the calibration data 'data'.
    // If no calibrations were found 0 is returned.
    // NOTE: The list must be destroyed by the caller.

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return 0;

    return SearchDistinctEntries("calibration", d->GetTableName());
}

//______________________________________________________________________________
Bool_t TCMySQLManager::AddDataSet(const Char_t* data, const Char_t* calibration, const Char_t* desc,
                                  Int_t first_run, Int_t last_run, Double_t* par, Int_t length,
                                  Bool_t skipChecks)
{
    // Create a new set of the calibration data 'data' with the calibration identifier
    // 'calibration' for the runs 'first_run' to 'last_run'. Use 'desc' as a
    // description. Read the 'length' parameters from 'par'.
    // Check run numbers if 'skipChecks' is kFALSE (default) and abort on errors.
    // Return kFALSE when an error occurred, otherwise kTRUE.

    TString table;

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return kFALSE;

    // do some checks concerning the run numbers
    if (!skipChecks)
    {
        //
        // check if first and last run exist
        //

        // check first run
        if (!SearchRunEntry(first_run, "run", table))
        {
            if (!fSilence) Error("AddDataSet", "First run has no valid run number!");
            return kFALSE;
        }

        // check last run
        if (!SearchRunEntry(last_run, "run", table))
        {
            if (!fSilence) Error("AddDataSet", "Last run has no valid run number!");
            return kFALSE;
        }

        //
        // check if the run range is ok
        //

        // check if first run is smaller than last run
        if (first_run > last_run)
        {
            if (!fSilence) Error("AddDataSet", "First run of set has to be smaller than last run!");
            return kFALSE;
        }

        // get number of existing sets
        Int_t nSet = GetNsets(data, calibration);

        // loop over sets
        for (Int_t i = 0; i < nSet; i++)
        {
            // get first and last runs of set
            Int_t setLow = GetFirstRunOfSet(data, calibration, i);
            Int_t setHigh = GetLastRunOfSet(data, calibration, i);

            // check if first run is member of this set
            if (first_run >= setLow && first_run <= setHigh)
            {
                if (!fSilence) Error("AddDataSet", "First run is already member of set %d", i);
                return kFALSE;
            }

            // check if last run is member of this set
            if (last_run >= setLow && last_run <= setHigh)
            {
                if (!fSilence) Error("AddDataSet", "Last run is already member of set %d", i);
                return kFALSE;
            }

            // check if sets are not overlapping
            if ((setLow >= first_run && setLow <= last_run) ||
                (setHigh >= first_run && setHigh <= last_run))
            {
                if (!fSilence) Error("AddDataSet", "Run overlap with set %d", i);
                return kFALSE;
            }
        }
    }

    //
    // create the set
    //

    // get the data table
    Char_t ctable[256];
    if (!SearchTable(data, ctable))
    {
        if (!fSilence) Error("AddDataSet", "No data table found!");
        return kFALSE;
    }

    // prepare the insert query
    TString ins_query_1 = TString::Format("INSERT INTO %s (calibration, description, first_run, last_run",
                                          ctable);
    TString ins_query_2 = TString::Format(" VALUES ( '%s', '%s', %d, %d",
                                          calibration, desc, first_run, last_run);

    // read all parameters and update the partial queries
    for (Int_t j = 0; j < length; j++)
    {
        // append parameter to query
        ins_query_1.Append(TString::Format(", par_%03d", j));
        ins_query_2.Append(TString::Format(", %.17g", par[j]));
    }

    // finalize query
    ins_query_1.Append(")");
    ins_query_2.Append(")");
    ins_query_1.Append(ins_query_2.Data());

    // write data to database
    Bool_t res = SendExec(ins_query_1.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("AddDataSet", "Could not add the set of '%s' for runs %d to %d!",
                             d->GetTitle(), first_run, last_run);
        return kFALSE;
    }
    else
    {
        if (!fSilence) Info("AddDataSet", "Added set of '%s' for runs %d to %d",
                            d->GetTitle(), first_run, last_run);
        return kTRUE;
    }
}

//______________________________________________________________________________
Bool_t TCMySQLManager::AddDataSet(const Char_t* data, const Char_t* calibration, const Char_t* desc,
                                  Int_t first_run, Int_t last_run, Double_t par)
{
    // Create a new set of the calibration data 'data' with the calibration identifier
    // 'calibration' for the runs 'first_run' to 'last_run'. Use 'desc' as a
    // description. Set all parameters to the value 'par'.
    // Return kFALSE when an error occurred, otherwise kTRUE.

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return kFALSE;

    // get maximum number of parameters
    Int_t length = d->GetSize();

    // create and fill parameter array
    Double_t par_array[length];
    for (Int_t i = 0; i < length; i++) par_array[i] = par;

    // set parameters
    return AddDataSet(data, calibration, desc, first_run, last_run, par_array, length);
}

//______________________________________________________________________________
Bool_t TCMySQLManager::AddSet(const Char_t* type, const Char_t* calibration, const Char_t* desc,
                              Int_t first_run, Int_t last_run, Double_t par)
{
    // Create new sets for the calibration type 'type' with the calibration identifier
    // 'calibration' for the runs 'first_run' to 'last_run'. Use 'desc' as a
    // description. Set all parameters to the value 'par'.
    // Return kFALSE when an error occurred, otherwise kTRUE.

    Bool_t ret = kTRUE;

    // create and fill parameter array
    Double_t par_array[TCConfig::kMaxCrystal];
    for (Int_t i = 0; i < TCConfig::kMaxCrystal; i++) par_array[i] = par;

    // get calibration type
    TCCalibType* t = GetCalibType(type);
    if (!t) return kFALSE;

    // get data list
    TList* data = t->GetData();

    // loop over calibration data of this calibration type
    TIter next(data);
    TCCalibData* d;
    while ((d = (TCCalibData*)next()))
    {
        // add set
        if (AddDataSet(d->GetName(), calibration, desc, first_run, last_run, par_array,
                       d->GetSize())) ret = kFALSE;
    }

    return ret;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::RemoveDataSet(const Char_t* data, const Char_t* calibration, Int_t set)
{
    // Remove the set 'set' from the calibration 'calibration' of the calibration data
    // 'data'.

    TString query;
    Char_t table[256];

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return kFALSE;

    // get the data table
    if (!SearchTable(data, table))
    {
        if (!fSilence) Error("RemoveDataSet", "No data table found!");
        return kFALSE;
    }

    // get the first run of the set
    Int_t first_run = GetFirstRunOfSet(data, calibration, set);

    // check first run
    if (!first_run)
    {
        if (!fSilence) Error("RemoveDataSet", "Could not delete set %d in '%s' of calibration '%s'!",
                             set, d->GetTitle(), calibration);
        return kFALSE;
    }

    // create the query
    query.Form("DELETE FROM %s WHERE "
               "calibration = '%s' AND "
               "first_run = %d",
               table, calibration, first_run);

    // read from database
    Bool_t res = SendExec(query.Data());

    // check result
    if (!res)
    {
        if (!fSilence) Error("RemoveDataSet", "Could not delete set %d in '%s' of calibration '%s'!",
                             set, d->GetTitle(), calibration);
        return kFALSE;
    }
    else
    {
        if (!fSilence) Info("RemoveDataSet", "Deleted set %d in '%s' of calibration '%s'",
                            set, d->GetTitle(), calibration);
        return kTRUE;
    }
}

//______________________________________________________________________________
Bool_t TCMySQLManager::RemoveSet(const Char_t* type, const Char_t* calibration, Int_t set)
{
    // Remove all sets 'set' from the calibration 'calibration' that are needed by the
    // calibration type 'type'.

    Bool_t ret = kTRUE;

    // get calibration type
    TCCalibType* t = GetCalibType(type);
    if (!t) return kFALSE;

    // get data list
    TList* data = t->GetData();

    // loop over calibration data of this calibration type
    TIter next(data);
    TCCalibData* d;
    while ((d = (TCCalibData*)next()))
    {
        // remove set of calibration data
        if (!RemoveDataSet(d->GetName(), calibration, set))
            ret = kFALSE;
    }

    return ret;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::SplitDataSet(const Char_t* data, const Char_t* calibration, Int_t set,
                                    Int_t lastRunFirstSet)
{
    // Split the set 'set' of the calibration data 'data' and the calibration identifier
    // 'calibration' into two sets. 'lastRunFirstSet' will be the last run of the first
    // set. The first run of the second set will be the next found run in the database.

    TString query;
    TString tmp;

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return kFALSE;

    // check if splitting run exists
    if (!SearchRunEntry(lastRunFirstSet, "run", tmp))
    {
        if (!fSilence) Error("SplitDataSet", "Splitting run has no valid run number!");
        return kFALSE;
    }

    // check if splitting run is in set
    Int_t first_run = GetFirstRunOfSet(data, calibration, set);
    Int_t last_run = GetLastRunOfSet(data, calibration, set);
    if (lastRunFirstSet < first_run || lastRunFirstSet > last_run)
    {
        if (!fSilence) Error("SplitDataSet", "Splitting run has to be in set!");
        return kFALSE;
    }

    //
    // get the first and last run of the second set
    //

    // create the query
    query.Form("SELECT run FROM %s "
               "WHERE run > %d "
               "ORDER by run LIMIT 1",
               TCConfig::kCalibMainTableName, lastRunFirstSet);

    // read from database
    TSQLResult* res = SendQuery(query.Data());
    if (!res)
    {
        if (!fSilence) Error("SplitDataSet", "Cannot find first run of second set!");
        return kFALSE;
    }

    // get row
    TSQLRow* row = res->Next();

    // check row
    if (!row)
    {
        if (!fSilence) Error("SplitDataSet", "Cannot find first run of second set!");
        delete res;
        return kFALSE;
    }

    // get the first run of the second set
    Int_t firstRunSecondSet = 0;
    const Char_t* field = row->GetField(0);
    if (!field)
    {
        if (!fSilence) Error("SplitDataSet", "Cannot find first run of second set!");
        return kFALSE;
    }
    else firstRunSecondSet = atoi(field);

    // clean-up
    delete row;
    delete res;

    // get the last run of the second set
    Int_t lastRunSecondSet = last_run;

    //
    // backup values
    //

    // backup the description
    Char_t desc[256];
    if (!SearchSetEntry(data, calibration, set, "description", desc))
    {
        if (!fSilence) Error("SplitDataSet", "Cannot read description of set %d in '%s' of calibration '%s'",
                             set, d->GetTitle(), calibration);
        return kFALSE;
    }

    // get number of parameters
    Int_t nPar = d->GetSize();

    // backup parameters
    Double_t par[nPar];
    if (!ReadParameters(data, calibration, set, par, nPar))
    {
        if (!fSilence) Error("SplitDataSet", "Cannot backup parameters of set %d in '%s' of calibration '%s'",
                             set, d->GetTitle(), calibration);
    }

    //
    // modify first set
    //

    // change the last run of the first set
    tmp.Form("%d", lastRunFirstSet);
    if (!ChangeSetEntry(data, calibration, set, "last_run", tmp.Data()))
    {
        if (!fSilence) Error("SplitDataSet", "Cannot change last run of set %d in '%s' of calibration '%s'",
                             set, d->GetTitle(), calibration);
        return kFALSE;
    }
    else
    {
        if (!fSilence) Info("SplitDataSet", "Changed last run of set %d in '%s' of calibration '%s' to %d",
                            set, d->GetTitle(), calibration, lastRunFirstSet);
    }

    //
    // add second set
    //

    if (!AddDataSet(data, calibration, desc, firstRunSecondSet, lastRunSecondSet, par, nPar))
    {
        if (!fSilence) Error("SplitDataSet", "Cannot split set %d in '%s' of calibration '%s'",
                             set, d->GetTitle(), calibration);
        return kFALSE;
    }

    return kTRUE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::SplitSet(const Char_t* type, const Char_t* calibration, Int_t set,
                                Int_t lastRunFirstSet)
{
    // Split all sets 'set' of the calibration type 'type' and the calibration identifier
    // 'calibration' into two sets. 'lastRunFirstSet' will be the last run of the first
    // set. The first run of the second set will be the next found run in the database.

    Bool_t ret = kTRUE;

    // get calibration type
    TCCalibType* t = GetCalibType(type);
    if (!t) return kFALSE;

    // get data list
    TList* data = t->GetData();

    // loop over calibration data of this calibration type
    TIter next(data);
    TCCalibData* d;
    while ((d = (TCCalibData*)next()))
    {
        // split set of calibration data
        if (!SplitDataSet(d->GetName(), calibration, set, lastRunFirstSet))
            ret = kFALSE;
    }

    return ret;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::MergeDataSets(const Char_t* data, const Char_t* calibration,
                                     Int_t set1, Int_t set2)
{
    // Merge the two adjacent sets 'set1' and 'set2' of the calibration data 'data' and the
    // calibration identifier 'calibration' into one set. Description and parameters
    // of 'set1' will be used in the merged set.
    // get the first run of the set

    Char_t tmp[256];

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return kFALSE;

    // check if the two sets are adjacent
    if (TMath::Abs(set2 - set1) != 1)
    {
        if (!fSilence) Error("MergeDataSets", "Only adjacent sets can be merged!");
        return kFALSE;
    }

    // get information of set 1
    Int_t firstSet1 = GetFirstRunOfSet(data, calibration, set1);
    Int_t lastSet1 = GetLastRunOfSet(data, calibration, set1);

    // get information of set 2
    Int_t firstSet2 = GetFirstRunOfSet(data, calibration, set2);
    Int_t lastSet2 = GetLastRunOfSet(data, calibration, set2);

    // check set 1
    if (!firstSet1 || !lastSet1)
    {
        if (!fSilence) Error("MergeDataSets", "Could not find set %d in '%s' of calibration '%s'!",
                             set1, d->GetTitle(), calibration);
        return kFALSE;
    }

    // check set 2
    if (!firstSet2 || !lastSet2)
    {
        if (!fSilence) Error("MergeDataSets", "Could not find set %d in '%s' of calibration '%s'!",
                             set2, d->GetTitle(), calibration);
        return kFALSE;
    }

    // delete set 2
    if (!RemoveDataSet(data, calibration, set2))
    {
        if (!fSilence) Error("MergeDataSets", "Could not remove set %d!", set2);
        return kFALSE;
    }

    // adjust run interval
    // check if set 2 was in front of set 1
    if (firstSet2 < firstSet1)
    {
        // adjust first run of remaining set
        sprintf(tmp, "%d", firstSet2);
        if (!ChangeSetEntry(data, calibration, set1-1, "first_run", tmp))
        {
            if (!fSilence) Error("MergeDataSets", "Could not adjust run interval of set %d!", set1-1);
            return kFALSE;
        }
    }
    else
    {
        // adjust last run of remaining set
        sprintf(tmp, "%d", lastSet2);
        if (!ChangeSetEntry(data, calibration, set1, "last_run", tmp))
        {
            if (!fSilence) Error("MergeDataSets", "Could not adjust run interval of set %d!", set1);
            return kFALSE;
        }
    }

    return kTRUE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ChangeRunBadScR(Int_t run, Int_t nbadscr, const Int_t* badscr, const Char_t* data)
{
    // ...

    // get short name (only last part of calib data, e.g. 'Data.Run.BadScR.NaI' --> 'NaI')
    data = strrchr(data, '.') + 1;

    // read old values from database -------------------------------------------

    // init return values
    TCBadScRElement** badscr_data = 0;
    Int_t ndata = 0;

    // read all data
    if (!ReadAllBadScR(run, badscr_data, ndata)) return kFALSE;


    // change scaler reads in list --------------------------------------------

    // init flag
    Bool_t isAdded = kFALSE;

    // loop over data
    for (Int_t d = 0; d < ndata; d++)
    {
        // check whether names match
        if (strcmp(badscr_data[d]->GetCalibData(), data) == 0)
        {
            // remove bad scaler reads
            badscr_data[d]->RemBad();

            // add new bad scaler reads
            badscr_data[d]->AddBad(nbadscr, badscr);

            // set flag
            isAdded = kTRUE;
        }
    }

    // check whether it was added already (if not append a new element)
    if (!isAdded)
    {
        // create new array
        TCBadScRElement** badscr_data_tmp = new TCBadScRElement*[ndata+1];

        // copy old pointers
        for (Int_t i = 0; i < ndata; i++)
            badscr_data_tmp[i] = badscr_data[i];

        // delete old array
        delete [] badscr_data;

        // set pointer
        badscr_data = badscr_data_tmp;

        // append the new data (do not need to set totnscr)
        badscr_data[ndata] = new TCBadScRElement(run, nbadscr, badscr);

        // set name
        badscr_data[ndata]->SetCalibData(data);

        // increment number of data
        ndata++;
    }


    // create query string -------------------------------------------------------

    // init query string
    TString s = "";

    // loop over data
    for (Int_t d = 0; d < ndata; d++)
    {
        // append name
        s.Append(badscr_data[d]->GetCalibData());
        s.Append(":");

        // get bad scaler read values
        Int_t nbadscr = badscr_data[d]->GetNBad();
        const Int_t* badscr = badscr_data[d]->GetBad();

        // loop over input bad scaler reads list
        for (Int_t i = 0; i < nbadscr; i++)
        {
            Char_t tmp[16];
            Int_t j = 0;

            // loop over subsequent values in order to detect series, i.e., badscr[k]+1 == badscr[k+1]
            for (j = 0; i + j < nbadscr - 1; j++)
                if (badscr[i+j] + 1 != badscr[i+j+1]) break;

            // check for series with more than 2 elements
            if (j > 1)
            {
                // separate first and last value of series with ':'
                sprintf(tmp, "%d-%d,", badscr[i], badscr[i+j]);
                i += j;
            }
            else
            {
                // separate with ','
                sprintf(tmp, "%d,", badscr[i]);
            }

            // append to query string
            s.Append(tmp);
        }

        // remove tailing comma
        s.Chop();

        // append data delimiter
        s.Append(";");
    }


    // finish ------------------------------------------------------------------

    // clean up
    for (Int_t i = 0; i < ndata; i++)
        delete badscr_data[i];
    delete [] badscr_data;

    // change list of bad scaler reads
    return ChangeRunEntries(run, run, "scr_bad", s.Data());
}

//______________________________________________________________________________
Bool_t TCMySQLManager::GetRunBadScR(Int_t run, Int_t& nbadscr, Int_t*& badscr, const Char_t* data)
{
    // Returns a list bad scaler reads for 'data' of the run 'run'. If data is
    // ommited all bad scaler read of the run are returned. The number of bad
    // scaler reads is stored to 'nbadscr' and the array of bad scaler reads is
    // stored to 'badscr'.
    // If there is no list of bad scaler reads for 'data' the NULL pointer for
    // 'badscr' is returnd.
    // NOTE: The array has to be destroyed by the caller.

    // get short name (only last part of calib data, e.g. 'Data.Run.BadScR.NaI' --> 'NaI')
    data = strrchr(data, '.') + 1;

    // declare temporary bad scr element
    TCBadScRElement* badscr_tmp = 0;

    // init return variables
    nbadscr = 0;
    badscr = 0;

    // init variables for database read call
    TCBadScRElement** badscr_data = 0;
    Int_t ndata = 0;

    // read bad scr for all data
    if (!ReadAllBadScR(run, badscr_data, ndata)) return kFALSE;

    // loop over data
    for (Int_t d = 0; d < ndata; d++)
    {
        // check whether names match
        if (!data || strcmp(badscr_data[d]->GetCalibData(), data) == 0)
        {
            // create bad scr element
            if (!badscr_tmp) badscr_tmp = new TCBadScRElement();

            // add bad scr
            badscr_tmp->AddBad(badscr_data[d]->GetNBad(), badscr_data[d]->GetBad());
        }

        // clean up
        delete badscr_data[d];
    }

    // clean up
    delete [] badscr_data;

    // check whether data was found
    if (badscr_tmp)
    {
        // set number of bad scaler reads
        nbadscr = badscr_tmp->GetNBad();

        // create array of bad scaler reads
        if (nbadscr) badscr = new Int_t[nbadscr];

        // copy values
        for (Int_t i = 0; i < nbadscr; i++)
            badscr[i] = badscr_tmp->GetBad()[i];

        // clean up
        delete badscr_tmp;
    }

    return kTRUE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ReadAllBadScR(Int_t run, TCBadScRElement**& badscr_data, Int_t& ndata)
{
    // Reads the all bad scaler read of the run 'run'. The number of data found
    // is stored to 'ndata'. The bad scaler reads of each data is stored to a
    // bad scaler read element
    // an array of pointers to bad scaler read elements
    // If the NULL
    // Returns kTRUE if the database readout was successful, kFALSE otherwise.
    // NOTE: The array has to be destroyed by the caller.

    // init result variables
    ndata = 0;
    badscr_data = 0;

    // init query return string
    TString tmp;

    // get run entry
    if (!SearchRunEntry(run, "scr_bad", tmp)) return kFALSE;

    // pointers to data tokens
    Char_t** data = 0;

    // get token for first data (e.g., "NaI:1,2,5-9\0")
    Char_t* datatok = strtok((Char_t*)tmp.Data(), ";");

    // loop over data tokens
    while (datatok)
    {
        // create new pointer array
        Char_t** data_tmp = new Char_t*[ndata + 1];

        // copy old pointers
        for (Int_t i = 0; i < ndata; i++)
            data_tmp[i] = data[i];

        // delete old pointer array
        delete [] data;

        // set pointer array to new pointer array
        data = data_tmp;

        // add new pointer
        data[ndata] = datatok;

        // increment data counter
        ndata++;

        // get next data token
        datatok = strtok(0, ";");
    }

    // create bad scaler read element array
    if (ndata) badscr_data = new TCBadScRElement*[ndata];

    // loop over data
    for (Int_t i = 0; i < ndata; i++)
    {
        // create new element
        badscr_data[i] = new TCBadScRElement();

        // get data name (e.g., "NaI\0")
        Char_t* name = strtok(data[i], ":");

        // set name
        badscr_data[i]->SetCalibData(name);

        // get first bad scr token (e.g., "1\0")
        Char_t* bad = strtok(0, ",");

        // loop over bad scr
        while (bad)
        {
            // detect series ('-' separated, e.g. "5-9\0")
            Char_t* loc = strchr(bad, (Int_t) '-');

            // check for series
            if (!loc)
            {
                // add single bad scr
                badscr_data[i]->AddBad(atoi(bad));
            }
            else
            {
                // get start value of series
                Int_t min = atoi(bad);

                // get stop value of series
                Int_t max = atoi(&loc[1]);

                // add series of bad scr
                for (Int_t j = min; j <= max; j++)
                    badscr_data[i]->AddBad(j);
            }

            // get next bad scr
            bad = strtok(0, ",");
        }
    }

    return kTRUE;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::MergeSets(const Char_t* type, const Char_t* calibration,
                                 Int_t set1, Int_t set2)
{
    // Merge all adjacent pairs of the sets 'set1' and 'set2' of the calibration type 'type' and the
    // calibration identifier 'calibration' into one set. Description and parameters
    // of 'set1' will be used in the merged set.
    // get the first run of the set

    Bool_t ret = kTRUE;

    // get calibration type
    TCCalibType* t = GetCalibType(type);
    if (!t) return kFALSE;

    // get data list
    TList* data = t->GetData();

    // loop over calibration data of this calibration type
    TIter next(data);
    TCCalibData* d;
    while ((d = (TCCalibData*)next()))
    {
        // merge sets of calibration data
        if (!MergeDataSets(d->GetName(), calibration, set1, set2))
            ret = kFALSE;
    }

    return ret;
}

//______________________________________________________________________________
Int_t TCMySQLManager::DumpRuns(TCContainer* container, Int_t first_run, Int_t last_run)
{
    // Dump the run information from run 'first_run' to run 'last_run' to
    // the CaLib container 'container'.
    // If first_run and last_run is zero all available runs will be dumped.
    // Return the number of dumped runs.

    TString query;
    TString tmp;

    // create the query
    if (!first_run && !last_run)
    {
        query.Form("SELECT run FROM %s "
                   "ORDER by run",
                   TCConfig::kCalibMainTableName);
    }
    else
    {
        query.Form("SELECT run FROM %s "
                   "WHERE run >= %d "
                   "AND run <= %d "
                   "ORDER by run",
                   TCConfig::kCalibMainTableName, first_run, last_run);
    }

    // read from database
    TSQLResult* res = SendQuery(query.Data());

    // create list for run numbers
    TList run_numbers;
    run_numbers.SetOwner(kTRUE);

    // read all rows/runs
    TSQLRow* r = res->Next();
    while (r)
    {
        run_numbers.Add(new TObjString(r->GetField(0)));
        delete r;
        r = res->Next();
    }

    // get number of runs
    Int_t nruns = run_numbers.GetSize();

    // loop over runs
    for (Int_t i = 0; i < nruns; i++)
    {
        // get next run
        TObjString* rn = (TObjString*) run_numbers.At(i);

        // get run number
        Int_t run_number = atoi(rn->GetString().Data());

        // add new run
        TCRun* run = container->AddRun(run_number);

        // set path
        if (SearchRunEntry(run_number, "path", tmp)) run->SetPath(tmp.Data());

        // set filename
        if (SearchRunEntry(run_number, "filename", tmp)) run->SetFileName(tmp.Data());

        // set time
        if (SearchRunEntry(run_number, "time", tmp)) run->SetTime(tmp.Data());

        // set description
        if (SearchRunEntry(run_number, "description", tmp)) run->SetDescription(tmp.Data());

        // set run_note
        if (SearchRunEntry(run_number, "run_note", tmp)) run->SetRunNote(tmp.Data());

        // set size
        if (SearchRunEntry(run_number, "size", tmp))
        {
            Long64_t size;
            sscanf(tmp.Data(), "%lld", &size);
            run->SetSize(size);
        }

        // set scaler reads
        if (SearchRunEntry(run_number, "scr_n", tmp)) run->SetNScalerReads(atoi(tmp.Data()));

        // set bad scaler reads
        if (SearchRunEntry(run_number, "scr_bad", tmp)) run->SetBadScalerReads(tmp.Data());

        // set target
        if (SearchRunEntry(run_number, "target", tmp)) run->SetTarget(tmp.Data());

        // set target polarization
        if (SearchRunEntry(run_number, "target_pol", tmp)) run->SetTargetPol(tmp.Data());

        // set target polarization degree
        if (SearchRunEntry(run_number, "target_pol_deg", tmp)) run->SetTargetPolDeg(atof(tmp.Data()));

        // set beam polarization
        if (SearchRunEntry(run_number, "beam_pol", tmp)) run->SetBeamPol(tmp.Data());

        // set beam polarization degree
        if (SearchRunEntry(run_number, "beam_pol_deg", tmp)) run->SetBeamPolDeg(atof(tmp.Data()));

        // user information
        if (!fSilence) Info("DumpRuns", "Dumped run %d", run_number);
    }

    // clean-up
    delete res;

    return nruns;
}

//______________________________________________________________________________
Int_t TCMySQLManager::DumpCalibrations(TCContainer* container, const Char_t* calibration,
                                       const Char_t* data)
{
    // Dump calibrations of the calibration data 'data' with the calibration
    // identifier 'calibration' to the CaLib container 'container'.
    // Return the number of dumped calibrations.

    Char_t tmp[256];

    // get data
    TCCalibData* d = GetCalibData(data);
    if (!d) return 0;

    // get number of parameters
    Int_t nPar = d->GetSize();

    // create the parameter array
    Double_t par[nPar];

    // get the number of sets
    Int_t nSet = GetNsets(data, calibration);

    // check calibration
    if (!nSet)
    {
        if (!fSilence) Error("DumpCalibrations", "No sets of '%s' of the calibration '%s' found!",
                             d->GetTitle(), calibration);
        return 0;
    }

    // loop over sets
    for (Int_t i = 0; i < nSet; i++)
    {
        // read parameters
        ReadParameters(data, calibration, i, par, nPar);

        // add the calibration
        TCCalibration* c = container->AddCalibration(calibration);

        // set calibration data
        c->SetCalibData(data);

        // set description
        GetDescriptionOfSet(data, calibration, i, tmp);
        c->SetDescription(tmp);

        // set first and last run
        c->SetFirstRun(GetFirstRunOfSet(data, calibration, i));
        c->SetLastRun(GetLastRunOfSet(data, calibration, i));

        // set fill time
        GetChangeTimeOfSet(data, calibration, i, tmp);
        c->SetChangeTime(tmp);

        // set parameters
        c->SetParameters(nPar, par);
    }

    // user information
    if (!fSilence) Info("DumpCalibrations", "Dumped %d sets of '%s' of the calibration '%s'",
                        nSet, d->GetTitle(), calibration);

    return nSet;
}

//______________________________________________________________________________
Int_t TCMySQLManager::DumpAllCalibrations(TCContainer* container, const Char_t* calibration)
{
    // Dump all calibrations with the calibration identifier 'calibration' to
    // the CaLib container 'container'.
    // Return the number of dumped calibrations.

    Int_t nDump = 0;

    // loop over calibration data
    TIter next(fData);
    TCCalibData* d;
    while ((d = (TCCalibData*)next()))
    {
        // dump calibrations
        if (Int_t nd = DumpCalibrations(container, calibration, d->GetName())) nDump += nd;
    }

    return nDump;
}

//______________________________________________________________________________
Int_t TCMySQLManager::ImportRuns(TCContainer* container)
{
    // Import all runs from the CaLib container 'container' to the database.
    // Return the number of imported runs.

    // get number of runs
    Int_t nRun = container->GetNRuns();

    // loop over runs
    Int_t nRunAdded = 0;
    for (Int_t i = 0; i < nRun; i++)
    {
        // get the run
        TCRun* r = container->GetRun(i);

        // prepare the insert query
        TString ins_query = TString::Format("INSERT INTO %s (run, path, filename, time, description, run_note, size, scr_n, scr_bad, "
                                            "target, target_pol, target_pol_deg, beam_pol, beam_pol_deg) "
                                            "VALUES ( "
                                            "%d, "
                                            "'%s', "
                                            "'%s', "
                                            "'%s', "
                                            "'%s', "
                                            "'%s', "
                                            "%lld, "
                                            "%d, "
                                            "'%s', "
                                            "'%s', "
                                            "'%s', "
                                            "%lf, "
                                            "'%s', "
                                            "%lf )",
                                            TCConfig::kCalibMainTableName,
                                            r->GetRun(),
                                            r->GetPath(),
                                            r->GetFileName(),
                                            r->GetTime(),
                                            r->GetDescription(),
                                            r->GetRunNote(),
                                            r->GetSize(),
                                            r->GetNScalerReads(),
                                            r->GetBadScalerReads(),
                                            r->GetTarget(),
                                            r->GetTargetPol(),
                                            r->GetTargetPolDeg(),
                                            r->GetBeamPol(),
                                            r->GetBeamPolDeg());

        // try to write data to database
        Bool_t res = SendExec(ins_query.Data());
        if (res == kFALSE)
        {
            Warning("ImportRuns", "Run %d could not be added to the database!",
                    r->GetRun());
        }
        else
        {
            if (!fSilence) Info("ImportRuns", "Added run %d to the database", r->GetRun());
            nRunAdded++;
        }
    }

    // user information
    if (!fSilence) Info("ImportRuns", "Added %d runs to the database", nRunAdded);

    return nRunAdded;
}

//______________________________________________________________________________
Int_t TCMySQLManager::ImportCalibrations(TCContainer* container, const Char_t* newCalibName,
                                         const Char_t* data)
{
    // Import all calibrations from the CaLib container 'container' to the database.
    // If 'newCalibName' is non-zero rename the calibration to 'newCalibName'
    // If 'data' is not 0 import only calibrations of the data 'data'.
    // Return the number of imported calibrations.

    // check data
    if (data && !GetCalibData(data)) return 0;

    // get number of calibrations
    Int_t nCalib = container->GetNCalibrations();

    // loop over calibrations
    Int_t nCalibAdded = 0;
    for (Int_t i = 0; i < nCalib; i++)
    {
        // get the calibration
        TCCalibration* c = container->GetCalibration(i);

        // skip unwanted calibration data
        if (data != 0 && strcmp(c->GetCalibData(), data)) continue;

        // add the set with new calibration identifer or the same
        const Char_t* calibration;
        if (newCalibName) calibration = newCalibName;
        else calibration = c->GetCalibration();

        TCCalibData* d = GetCalibData(c->GetCalibData());
        if (!d) continue;

        // add the set
        if (AddDataSet(c->GetCalibData(), calibration, c->GetDescription(),
                       c->GetFirstRun(), c->GetLastRun(), c->GetParameters(), c->GetNParameters(), kTRUE))
        {
            if (!fSilence) Info("ImportCalibrations", "Added calibration '%s' of '%s' to the database",
                                calibration, d->GetTitle());
            nCalibAdded++;
        }
        else
        {
            if (!fSilence) Error("ImportCalibrations", "Calibration '%s' of '%s' could not be added to the database!",
                                 calibration, d->GetTitle());
        }
    }

    // user information
    if (!fSilence) Info("ImportCalibrations", "Added %d calibrations to the database", nCalibAdded);

    return nCalibAdded;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::CloneCalibration(const Char_t* calibration, const Char_t* newCalibrationName,
                                        const Char_t* newDesc, Int_t new_first_run, Int_t new_last_run)
{
    // Create a clone of the calibration 'calibration' using 'newCalibrationName'
    // as the new calibration name and 'newDesc' as the calibration description.
    // If 'new_first_run' and 'new_last_run' are both equal to zero for each
    // calibration data all sets are copied.
    // Otherwise, if 'new_first_run' and 'new_last_run' are not equal to zero,
    // for each calibration data one set from 'new_first_run' to 'new_last_run'
    // is created with the values of the last set of the original calibration.
    // Return kFALSE if an error occurred, otherwise kTRUE.

    // check if original calibration exists
    if (!ContainsCalibration(calibration))
    {
        if (!fSilence) Error("CloneCalibration", "Original calibration '%s' does not exist!", calibration);
        return kFALSE;
    }

    // set "true clone" flag
    Bool_t true_clone = (new_first_run == 0 && new_last_run == 0) ? kTRUE : kFALSE;

    // loop over calibration data
    TIter next(fData);
    TCCalibData* d;
    Bool_t error = kFALSE;
    while ((d = (TCCalibData*)next()))
    {
        // get number of original sets
        Int_t nSets = GetNsets(d->GetName(), calibration);

        // check for true clone
        if (true_clone)
        {
            // loop over sets
            for (Int_t i = 0; i < nSets; i++)
            {
                // get run range of set
                Int_t first_run = GetFirstRunOfSet(d->GetName(), calibration, i);
                Int_t last_run = GetLastRunOfSet(d->GetName(), calibration, i);

                // read parameters
                Double_t par[d->GetSize()];
                if (ReadParameters(d->GetName(), calibration, i, par, d->GetSize()))
                {
                    // add new set
                    if (!AddDataSet(d->GetName(), newCalibrationName, newDesc, first_run, last_run, par, d->GetSize()))
                    {
                        if (!fSilence) Error("CloneCalibration", "Could not clone calibration data '%s'!", d->GetName());
                        error = kTRUE;
                    }
                }
                else
                {
                    if (!fSilence) Error("CloneCalibration", "Could not read original data '%s'!", d->GetName());
                    error = kTRUE;
                }
            }
        }
        else
        {
            // read parameters
            Double_t par[d->GetSize()];
            if (ReadParameters(d->GetName(), calibration, nSets-1, par, d->GetSize()))
            {
                // add new set
                if (!AddDataSet(d->GetName(), newCalibrationName, newDesc, new_first_run, new_last_run, par, d->GetSize()))
                {
                    if (!fSilence) Error("CloneCalibration", "Could not clone calibration data '%s'!", d->GetName());
                    error = kTRUE;
                }
            }
            else
            {
                if (!fSilence) Error("CloneCalibration", "Could not read original data '%s'!", d->GetName());
                error = kTRUE;
            }
        }
    }

    // check if an error occurred
    if (error)
    {
        // remove partially cloned calibrations when an error occurred
        RemoveAllCalibrations(newCalibrationName);
        return kFALSE;
    }
    else
    {
        return kTRUE;
    }
}

//______________________________________________________________________________
void TCMySQLManager::Export(const Char_t* filename, Int_t first_run, Int_t last_run,
                            const Char_t* calibration)
{
    // Export run and/or calibration data to the ROOT file 'filename'
    //
    // If 'first_run' is non-zero and 'last_run' is non-zero run information from run
    // 'first_run' to run 'last_run' is exported.
    // If 'first_run' is zero and 'last_run' is zero all run information is exported.
    // If 'first_run' is -1 or 'last_run' is -1 no run information is exported.
    //
    // If 'calibration' is non-zero the calibration with the identifier 'calibration'
    // is exported.

    // create new container
    TCContainer* container = new TCContainer(TCConfig::kCaLibDumpName);

    // dump runs to container
    if (first_run != -1 && last_run != -1)
    {
        DumpRuns(container, first_run, last_run);
        if (!fSilence)
        {
            if (container->GetNRuns())
                Info("Export", "Dumped %d runs to '%s'", container->GetNRuns(), filename);
            else
                Error("Export", "No runs were dumped to '%s'!", filename);
        }
    }

    // dump runs to container
    if (calibration)
    {
        DumpAllCalibrations(container, calibration);
        if (!fSilence)
        {
            if (container->GetNCalibrations())
                Info("Export", "Dumped %d calibrations to '%s'", container->GetNCalibrations(), filename);
            else
                Error("Export", "No calibrations were dumped to '%s'!", filename);
        }
    }

    // save container to ROOT file
    container->Save(filename, fSilence);

    // clean-up
    delete container;
}

//______________________________________________________________________________
TCContainer* TCMySQLManager::LoadContainer(const Char_t* filename)
{
    // Load the CaLib container from the ROOT file 'filename'.
    // Return the container if found, otherwise 0.
    // NOTE: The container must be destroyed by the caller.

    // try to open the ROOT file
    TFile* f = new TFile(filename);
    if (!f)
    {
        if (!fSilence) Error("LoadContainer", "Could not open the ROOT file '%s'!", filename);
        return 0;
    }
    if (f->IsZombie())
    {
        if (!fSilence) Error("LoadContainer", "Could not open the ROOT file '%s'!", filename);
        return 0;
    }

    // to load the CaLib container
    TCContainer* c_orig = (TCContainer*) f->Get(TCConfig::kCaLibDumpName);
    if (!c_orig)
    {
        if (!fSilence) Error("LoadContainer", "No CaLib container found in ROOT file '%s'!", filename);
        delete f;
        return 0;
    }

    // clone the container
    TCContainer* c = (TCContainer*) c_orig->Clone();

    // check container format
    if (c->GetVersion() != TCConfig::kContainerFormatVersion)
    {
        if (!fSilence) Error("LoadContainer", "Cannot load CaLib container format version %d - "
                                              "use the corresponding CaLib version instead!", c->GetVersion());
        return 0;
    }

    // clean-up
    delete f;

    return c;
}

//______________________________________________________________________________
void TCMySQLManager::Import(const Char_t* filename, Bool_t runs, Bool_t calibrations,
                            const Char_t* newCalibName)
{
    // Import run and/or calibration data from the ROOT file 'filename'
    //
    // If 'runs' is kTRUE all run information is imported.
    // If 'calibrations' is kTRUE all calibration information is imported.
    // If 'newCalibName' is non-zero rename the calibration to 'newCalibName'

    // try to load the container
    TCContainer* c = LoadContainer(filename);
    if (!c)
    {
        if (!fSilence) Error("Import", "CaLib container could not be loaded from '%s'", filename);
        return;
    }

    // import runs
    if (runs)
    {
        // get number of runs
        Int_t nRun = c->GetNRuns();

        // check if some runs were found
        if (nRun)
        {
            // ask for user confirmation
            Char_t answer[256];
            if (fDBType == kSQLite)
            {
                printf("\n%d runs were found in the ROOT file '%s'\n"
                       "They will be added to the database '%s'\n",
                       nRun, filename, fDB->GetDB());
            }
            else
            {
                printf("\n%d runs were found in the ROOT file '%s'\n"
                       "They will be added to the database '%s' on '%s'\n",
                       nRun, filename, fDB->GetDB(), fDB->GetHost());
            }
            printf("Are you sure to continue? (yes/no) : ");
            Int_t ret = scanf("%s", answer);
            if (strcmp(answer, "yes"))
            {
                printf("Aborted.\n");
            }
            else
            {
                // import all runs
                ImportRuns(c);
            }
        }
        else
        {
            if (!fSilence) Error("Import", "No runs were found in ROOT file '%s'!", filename);
        }
    }

    // import calibrations
    if (calibrations)
    {
        // get number of calibrations
        Int_t nCalib = c->GetNCalibrations();

        // check if some calibrations were found
        if (nCalib)
        {
            // get name of calibrations
            const Char_t* calibName = c->GetCalibration(0)->GetCalibration();

            // ask for user confirmation
            Char_t answer[256];
            if (fDBType == kSQLite)
            {
                printf("\n%d calibrations named '%s' were found in the ROOT file '%s'\n"
                       "They will be added to the database '%s'\n",
                       nCalib, calibName, filename, fDB->GetDB());
            }
            else
            {
                printf("\n%d calibrations named '%s' were found in the ROOT file '%s'\n"
                       "They will be added to the database '%s' on '%s'\n",
                       nCalib, calibName, filename, fDB->GetDB(), fDB->GetHost());
            }
            if (newCalibName) printf("The calibrations will be renamed to '%s'\n", newCalibName);
            printf("Are you sure to continue? (yes/no) : ");
            Int_t ret = scanf("%s", answer);
            if (strcmp(answer, "yes"))
            {
                printf("Aborted.\n");
            }
            else
            {
                // import all runs
                ImportCalibrations(c, newCalibName);
            }
        }
        else
        {
            if (!fSilence) Error("Import", "No calibrations were found in ROOT file '%s'!", filename);
        }
    }

    // clean-up
    delete c;
}

//______________________________________________________________________________
Bool_t TCMySQLManager::ExportDatabase(const Char_t* filename)
{
    // Export the complete database to the SQLite database 'filename'.

    Char_t tmp[256];
    Char_t fn[256];

    // expand filename
    Char_t* fnt = gSystem->ExpandPathName(filename);
    strcpy(fn, fnt);
    delete fnt;

    // check if file exists
    if (!gSystem->AccessPathName(fn))
    {
        if (!fSilence) Error("ExportDatabase", "File '%s' exists already!", fn);
        return kFALSE;
    }

    // open SQLite database
    sprintf(tmp, "sqlite://%s", fn);
    TSQLServer* db = TSQLServer::Connect(tmp, "", "");

    // check DB connection
    if (!db || db->IsZombie())
    {
        if (!fSilence) Error("ExportDatabase", "Cannot connect to the database '%s'!", fn);
        return kFALSE;
    }

    // create new container
    TCContainer* container = new TCContainer(TCConfig::kCaLibDumpName);

    // dump runs to container
    Int_t nRunExp = TCMySQLManager::GetManager()->DumpRuns(container);
    if (!fSilence) Info("ExportDatabase", "Exported %d runs to '%s'", nRunExp, fn);

    // get all calibrations
    TList* c = TCMySQLManager::GetManager()->GetAllCalibrations();

    // dump all calibrations to container
    Int_t nCalibExp = 0;
    Int_t nCalib = c->GetSize();
    if (!fSilence) Info("ExportDatabase", "Exporting %d calibrations to '%s'", nCalib, fn);
    for (Int_t i = 0; i < nCalib; i++)
    {
        TObjString* s = (TObjString*) c->At(i);
        nCalibExp += TCMySQLManager::GetManager()->DumpAllCalibrations(container, s->GetString().Data());
        if (!fSilence) Info("ExportDatabase", "Exported calibration '%s'", s->GetString().Data());
    }

    // backup original database config
    TSQLServer* db_orig = fDB;
    ServerType_t type_orig = fDBType;

    // configure db connection to SQLite database
    fDB = db;
    fDBType = kSQLite;

    // init the database
    InitDatabase(kFALSE);

    // import runs
    Int_t nRunImp = TCMySQLManager::GetManager()->ImportRuns(container);

    // import calibrations
    Int_t nCalibImp = TCMySQLManager::GetManager()->ImportCalibrations(container);

    // restore original db connection
    fDB = db_orig;
    fDBType = type_orig;

    // clean-up
    delete c;
    delete container;
    delete db;

    // make some checks
    if (nRunImp != nRunExp)
    {
        if (!fSilence) Error("ExportDatabase", "Could not export all runs!");
        return kFALSE;
    }
    if (nCalibImp != nCalibExp)
    {
        if (!fSilence) Error("ExportDatabase", "Could not export all calibrations!");
        return kFALSE;
    }

    return kTRUE;
}

