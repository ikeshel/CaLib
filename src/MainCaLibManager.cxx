/*************************************************************************
 * Author: Dominik Werthmueller, Thomas Strub
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// CaLibManager                                                         //
//                                                                      //
// Manage a CaLib database.                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <ncurses.h>
#include <signal.h>

#include "TList.h"
#include "TObjString.h"
#include "THashList.h"

#include "TCMySQLManager.h"
#include "TCCalibType.h"
#include "TCContainer.h"
#include "TCCalibData.h"

#define KEY_ENTER_MINE 13
#define KEY_ESC 27
#define MAX_SCR_BAD_STRING 30000

// locally used enum for run entries
enum ERunEntry
{
    kPATH,
    kTARGET,
    kTARGET_POL,
    kTARGET_POL_DEG,
    kBEAM_POL,
    kBEAM_POL_DEG
};
typedef ERunEntry RunEntry_t;

// global variables
Int_t gNrow;
Int_t gNcol;
Char_t gCalibration[256];
TCCalibType* gCalibrationType;
TCCalibData* gCalibrationData;
Char_t gFinishMessage[256];

// function prototypes
void MainMenu();
void RunEditor();
void CalibEditor();
void PrintStatusMessage(const Char_t* message);
Int_t SelectCalibration(const Char_t* mMsg = 0, Int_t active = 0);
Int_t SelectCalibrationData();
Int_t SelectCalibrationType();
void Administration();
void WriteTableEntry(WINDOW* win, const Char_t* str, Int_t colLength, Int_t att);
Int_t ShowMenu(const Char_t* title, const Char_t* message,
               Int_t nEntries, Char_t* entries[], Int_t active = 0);

//______________________________________________________________________________
void Finish(Int_t sig)
{
    // Exit the program.

    fclose(stderr);
    endwin();
    if (strcmp(gFinishMessage, "")) printf("%s\n", gFinishMessage);
    exit(0);
}

//______________________________________________________________________________
void DrawHeader(const Char_t* title = "CaLib Manager")
{
    // Draw a header.

    Char_t tmp[256];

    // turn on bold
    attron(A_BOLD);

    // draw lines
    for (Int_t i = 0; i < gNcol; i++)
    {
        mvprintw(0, i, "#");
        mvprintw(2, i, "#");
    }

    // draw lines
    for (Int_t i = 1; i < 3; i++)
    {
        mvprintw(i, 0, "#");
        mvprintw(i, gNcol-1, "#");
    }

    // draw title
    if (TCMySQLManager::GetManager()->GetDBType() == kSQLite)
    {
        sprintf(tmp, "%s using %s with CaLib %s",
                title,
                TCMySQLManager::GetManager()->GetDBName(),
                TCConfig::kCaLibVersion);
    }
    else
    {
        sprintf(tmp, "%s using %s@%s with CaLib %s",
                title,
                TCMySQLManager::GetManager()->GetDBName(),
                TCMySQLManager::GetManager()->GetDBHost(),
                TCConfig::kCaLibVersion);
    }
    mvprintw(1, (gNcol-strlen(tmp)) / 2, "%s", tmp);

    // turn off bold
    attroff(A_BOLD);
}

//______________________________________________________________________________
void DrawMenu(const Char_t* title, const Char_t* message,
               Int_t nEntries, Char_t* entries[], Int_t active)
{
    // Draw the main menu.

    // clear the screen
    clear();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "%s", title);
    attroff(A_UNDERLINE);

    // draw message
    mvprintw(6, 2, "%s:", message);

    // calculate the entries window width
    Int_t e_window_w = 0;
    for (Int_t i = 0; i < nEntries; i++)
        e_window_w = TMath::Max(e_window_w, (Int_t)strlen(entries[i]));

    // create the entries window
    WINDOW* entries_window = newpad(nEntries, e_window_w);

    // draw entries
    for (Int_t i = 0; i < nEntries; i++)
    {
        // move cursor
        wmove(entries_window, i, 0);

        // print entry
        if (i == active) wattron(entries_window, A_REVERSE);
        wprintw(entries_window, "%s", entries[i]);
        if (i == active) wattroff(entries_window, A_REVERSE);
    }

    // user information
    PrintStatusMessage("Use UP and DOWN keys to select - hit ENTER or RIGHT key to confirm - use LEFT key to go back");

    // calculate mininum pad row
    Int_t minp = active+1 > gNrow-3-8 ? active - (gNrow-3-8) : 0;

    // refresh window
    refresh();
    prefresh(entries_window, minp, 0, 8, 2, gNrow-3, gNcol-3);
    move(gNrow-1, gNcol-1);

}

//______________________________________________________________________________
Int_t ShowMenu(const Char_t* title, const Char_t* message,
               Int_t nEntries, Char_t* entries[], Int_t active)
{
    // Navigate the main menu.

    // draw
    DrawMenu(title, message, nEntries, entries, active);

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        //
        // decide what to do
        //

        // go to previous menue
        if (c == KEY_LEFT)
        {
            return -1;
        }
        // go to previous entry
        if (c == KEY_UP)
        {
            if (active > 0)
            {
                DrawMenu(title, message, nEntries, entries, --active);
            }
        }
        // go to next entry
        else if (c == KEY_DOWN)
        {
            if (active < nEntries-1)
            {
                DrawMenu(title, message, nEntries, entries, ++active);
            }
        }
        // choose entry
        else if (c == KEY_ENTER_MINE || c == KEY_ENTER || c == KEY_RIGHT)
        {
            return active;
        }
    }
}

//______________________________________________________________________________
void WriteTableEntry(WINDOW* win, const Char_t* str, Int_t colLength,
                     Int_t att = A_NORMAL)
{
    // Write the table entry 'str' to the window 'win' using the column length 'colLength'
    // at the current cursor position.

    // set attribute
    wattron(win, att);

    // check if string is empty
    if (strlen(str) == 0)
    {
        // write '-' with spaces
        wprintw(win, "-");
        for (Int_t i = 1; i < colLength+4; i++) wprintw(win, " ");
    }
    else
    {
        // write entry
        wprintw(win, "%s", str);

        // fill with spaces
        for (Int_t i = strlen(str); i < colLength+4; i++) wprintw(win, " ");
    }

    // unset attribute
    wattroff(win, att);
}

//______________________________________________________________________________
WINDOW* FormatRunTable(TCContainer& c, Int_t* outColLengthTot, WINDOW** outHeader)
{
    // Format the run table using the runs in 'c' and return the created window.
    // Save the maximum column length to 'outColLengthTot'.
    // Save the header window in 'outHeader'.

    // get number of runs
    Int_t nRuns = c.GetNRuns();

    // define col headers
    const Char_t* colHead[] = { "Run", "Path", "File name", "Time", "Description",
                                "Run note", "File size", "Scaler reads", "Bad scaler reads", "Target", "Target pol.",
                                "Target pol. deg.", "Beam pol.", "Beam pol. deg." };

    // determine maximum col lengths
    Char_t tmp_str[65536];
    const Char_t* tmp = 0;
    UInt_t colLength[14];

    // loop over headers
    for (Int_t i = 0; i < 14; i++) colLength[i] = strlen(colHead[i]);

    // loop over runs
    for (Int_t i = 0; i < nRuns; i++)
    {
        // run number
        sprintf(tmp_str, "%d", c.GetRun(i)->GetRun());
        if (strlen(tmp_str) > colLength[0]) colLength[0] = strlen(tmp_str);

        // path
        tmp = c.GetRun(i)->GetPath();
        if (strlen(tmp) > colLength[1]) colLength[1] = strlen(tmp);

        // file name
        tmp = c.GetRun(i)->GetFileName();
        if (strlen(tmp) > colLength[2]) colLength[2] = strlen(tmp);

        // time
        tmp = c.GetRun(i)->GetTime();
        if (strlen(tmp) > colLength[3]) colLength[3] = strlen(tmp);

        // description
        tmp = c.GetRun(i)->GetDescription();
        if (strlen(tmp) > colLength[4]) colLength[4] = strlen(tmp);

        // run note
        tmp = c.GetRun(i)->GetRunNote();
        if (strlen(tmp) > colLength[5]) colLength[5] = strlen(tmp);

        // size
        sprintf(tmp_str, "%lld", c.GetRun(i)->GetSize());
        if (strlen(tmp_str) > colLength[6]) colLength[6] = strlen(tmp_str);

        // scaler reads
        sprintf(tmp_str, "%d", c.GetRun(i)->GetNScalerReads());
        if (strlen(tmp_str) > colLength[7]) colLength[7] = strlen(tmp_str);

        // bad scaler reads
        tmp = c.GetRun(i)->GetBadScalerReads();
        if (strlen(tmp) > colLength[8]) colLength[8] = strlen(tmp);
        if (colLength[8] > MAX_SCR_BAD_STRING) colLength[8] = MAX_SCR_BAD_STRING;

        // target
        tmp = c.GetRun(i)->GetTarget();
        if (strlen(tmp) > colLength[9]) colLength[9] = strlen(tmp);

        // target polarization
        tmp = c.GetRun(i)->GetTargetPol();
        if (strlen(tmp) > colLength[10]) colLength[10] = strlen(tmp);

        // target polarization degree
        sprintf(tmp_str, "%f", c.GetRun(i)->GetTargetPolDeg());
        if (strlen(tmp_str) > colLength[11]) colLength[11] = strlen(tmp_str);

        // beam polarization
        tmp = c.GetRun(i)->GetBeamPol();
        if (strlen(tmp) > colLength[12]) colLength[12] = strlen(tmp);

        // beam polarization degree
        sprintf(tmp_str, "%f", c.GetRun(i)->GetBeamPolDeg());
        if (strlen(tmp_str) > colLength[13]) colLength[13] = strlen(tmp_str);
    }

    // calculate the maximum col length (14*4 spaces)
    Int_t colLengthTot = 14*4;
    for (Int_t i = 0; i < 14; i++) colLengthTot += colLength[i];
    *outColLengthTot = colLengthTot;

    // create the table and the header window
    WINDOW* table = newpad(nRuns, colLengthTot);
    WINDOW* header = newpad(1, colLengthTot);

    // add header content
    wmove(header, 0, 0);
    for (Int_t i = 0; i < 14; i++)
        WriteTableEntry(header, colHead[i], colLength[i], A_BOLD);

    // add table content
    for (Int_t i = 0; i < nRuns; i++)
    {
        // move cursor
        wmove(table, i, 0);

        // run number
        sprintf(tmp_str, "%d", c.GetRun(i)->GetRun());
        WriteTableEntry(table, tmp_str, colLength[0]);

        // path
        WriteTableEntry(table, c.GetRun(i)->GetPath(), colLength[1]);

        // file name
        WriteTableEntry(table, c.GetRun(i)->GetFileName(), colLength[2]);

        // time
        WriteTableEntry(table, c.GetRun(i)->GetTime(), colLength[3]);

        // description
        WriteTableEntry(table, c.GetRun(i)->GetDescription(), colLength[4]);

        // run note
        WriteTableEntry(table, c.GetRun(i)->GetRunNote(), colLength[5]);

        // size
        sprintf(tmp_str, "%lld", c.GetRun(i)->GetSize());
        WriteTableEntry(table, tmp_str, colLength[6]);

        // scaler reads
        sprintf(tmp_str, "%d", c.GetRun(i)->GetNScalerReads());
        WriteTableEntry(table, tmp_str, colLength[7]);

        // bad scaler reads
        strncpy(tmp_str, c.GetRun(i)->GetBadScalerReads(), MAX_SCR_BAD_STRING);
        tmp_str[MAX_SCR_BAD_STRING] = '\0';
        WriteTableEntry(table, tmp_str, colLength[8]);

        // target
        WriteTableEntry(table, c.GetRun(i)->GetTarget(), colLength[9]);

        // target polarization
        WriteTableEntry(table, c.GetRun(i)->GetTargetPol(), colLength[10]);

        // target polarization degree
        sprintf(tmp_str, "%f", c.GetRun(i)->GetTargetPolDeg());
        WriteTableEntry(table, tmp_str, colLength[11]);

        // beam polarization
        WriteTableEntry(table, c.GetRun(i)->GetBeamPol(), colLength[12]);

        // beam polarization degree
        sprintf(tmp_str, "%f", c.GetRun(i)->GetBeamPolDeg());
        WriteTableEntry(table, tmp_str, colLength[13]);
    }

    *outHeader = header;

    return table;
}

//______________________________________________________________________________
WINDOW* FormatCalibTable(TCContainer& c, Int_t* outColLengthTot, WINDOW** outHeader,
                         Bool_t withPar)
{
    // Format the calibration table using the calibrations in 'c' and return
    // the created window.
    // Save the maximum column length to 'outColLengthTot'.
    // Save the header window in 'outHeader'.
    // If 'withPar' is kTRUE include also the parameters (gCalibrationData has to be set).

    // get number of calibrations
    Int_t nCalib = c.GetNCalibrations();

    // get number of paramters
    Int_t nPar = 0;
    if (withPar) nPar = gCalibrationData->GetSize();

    // define col headers
    Int_t nCol = 5+nPar;
    Char_t* colHead[nCol];
    for (Int_t i = 0; i < nCol; i++) colHead[i] = new Char_t[16];
    strcpy(colHead[0], "Set");
    strcpy(colHead[1], "First run");
    strcpy(colHead[2], "Last run");
    strcpy(colHead[3], "Change time");
    strcpy(colHead[4], "Description");
    for (Int_t i = 0; i < nPar; i++) sprintf(colHead[5+i], "Par. %03d", i);

    // determine maximum col lengths
    Char_t tmp_str[256];
    const Char_t* tmp = 0;
    UInt_t colLength[nCol];

    // loop over headers
    for (Int_t i = 0; i < nCol; i++) colLength[i] = strlen(colHead[i]);

    // loop over calibrations
    for (Int_t i = 0; i < nCalib; i++)
    {
        // set
        sprintf(tmp_str, "%d", i);
        if (strlen(tmp_str) > colLength[0]) colLength[0] = strlen(tmp_str);

        // first run
        sprintf(tmp_str, "%d", c.GetCalibration(i)->GetFirstRun());
        if (strlen(tmp_str) > colLength[1]) colLength[1] = strlen(tmp_str);

        // last run
        sprintf(tmp_str, "%d", c.GetCalibration(i)->GetLastRun());
        if (strlen(tmp_str) > colLength[2]) colLength[2] = strlen(tmp_str);

        // change time
        tmp = c.GetCalibration(i)->GetChangeTime();
        if (strlen(tmp) > colLength[3]) colLength[3] = strlen(tmp);

        // description
        tmp = c.GetCalibration(i)->GetDescription();
        if (strlen(tmp) > colLength[4]) colLength[4] = strlen(tmp);

        // parameters
        Double_t* par = c.GetCalibration(i)->GetParameters();
        for (Int_t j = 0; j < nPar; j++)
        {
            sprintf(tmp_str, "%lf", par[j]);
            if (strlen(tmp_str) > colLength[5+j]) colLength[5+j] = strlen(tmp_str);
        }
    }

    // calculate the maximum col length (nCol*4 spaces)
    Int_t colLengthTot = nCol*4;
    for (Int_t i = 0; i < nCol; i++) colLengthTot += colLength[i];
    *outColLengthTot = colLengthTot;

    // create the table and the header window
    WINDOW* table = newpad(nCalib, colLengthTot);
    WINDOW* header = newpad(1, colLengthTot);

    // add header content
    wmove(header, 0, 0);
    for (Int_t i = 0; i < nCol; i++)
        WriteTableEntry(header, colHead[i], colLength[i], A_BOLD);

    // add table content
    for (Int_t i = 0; i < nCalib; i++)
    {
        // move cursor
        wmove(table, i, 0);

        // set
        sprintf(tmp_str, "%d", i);
        WriteTableEntry(table, tmp_str, colLength[0]);

        // first run
        sprintf(tmp_str, "%d", c.GetCalibration(i)->GetFirstRun());
        WriteTableEntry(table, tmp_str, colLength[1]);

        // last run
        sprintf(tmp_str, "%d", c.GetCalibration(i)->GetLastRun());
        WriteTableEntry(table, tmp_str, colLength[2]);

        // change time
        WriteTableEntry(table, c.GetCalibration(i)->GetChangeTime(), colLength[3]);

        // description
        WriteTableEntry(table, c.GetCalibration(i)->GetDescription(), colLength[4]);

        // parameters
        Double_t* par = c.GetCalibration(i)->GetParameters();
        for (Int_t j = 0; j < nPar; j++)
        {
            sprintf(tmp_str, "%lf", par[j]);
            WriteTableEntry(table, tmp_str, colLength[5+j]);
        }
    }

    // clean-up
    for (Int_t i = 0; i < nCol; i++) delete [] colHead[i];

    *outHeader = header;

    return table;
}

//______________________________________________________________________________
void PrintStatusMessage(const Char_t* message)
{
    // Print a status message.

    // print message
    mvprintw(gNrow-1, 0, message);

    // delete rest of line with spaces
    for (Int_t i = strlen(message); i < gNcol; i++) printw(" ");
}

//______________________________________________________________________________
void SplitSet()
{
    // Split a calibration set into two.
    // (gCalibration and gCalibrationType have to be set)

    Char_t answer[16];
    Int_t set;
    Int_t last_run;

    // echo input
    echo();

    // ask set
    mvprintw(gNrow-10, 2, "Enter number of set to split               : ");
    scanw((Char_t*)"%d", &set);

    // ask last run
    mvprintw(gNrow-9, 2, "Enter last run of first set                : ");
    scanw((Char_t*)"%d", &last_run);

    // ask confirmation
    mvprintw(gNrow-7, 2, "Splitting set %d after run %d", set, last_run);
    mvprintw(gNrow-5, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(gNrow-3, 2, "Aborted.");
    }
    else
    {
        // split set
        Bool_t ret = TCMySQLManager::GetManager()->SplitSet(gCalibrationType->GetName(), gCalibration,
                                                            set, last_run);

        // check return value
        if (ret)
            mvprintw(gNrow-3, 2, "Split set %d successfully into two sets", set);
        else
            mvprintw(gNrow-3, 2, "There was an error during splitting set %d!", set);
    }

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to the calibration editor)
    return;
}

//______________________________________________________________________________
void MergeSets()
{
    // Merge two calibration sets into one.
    // (gCalibration and gCalibrationType have to be set)

    Char_t answer[16];
    Int_t set1;
    Int_t set2;

    // echo input
    echo();

    // ask set 1
    mvprintw(gNrow-10, 2, "Enter number of first set (use its paramters) : ");
    scanw((Char_t*)"%d", &set1);

    // ask set 2
    mvprintw(gNrow-9, 2, "Enter number of second set                    : ");
    scanw((Char_t*)"%d", &set2);

    // ask confirmation
    mvprintw(gNrow-7, 2, "Merging sets %d and %d using the paramters of set %d", set1, set2, set1);
    mvprintw(gNrow-5, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(gNrow-3, 2, "Aborted.");
    }
    else
    {
        // merge sets
        Bool_t ret = TCMySQLManager::GetManager()->MergeSets(gCalibrationType->GetName(), gCalibration,
                                                            set1, set2);

        // check return value
        if (ret)
            mvprintw(gNrow-3, 2, "Merged sets %d and %d successfully into one set", set1, set2);
        else
            mvprintw(gNrow-3, 2, "There was an error during merging sets %d and %d!", set1, set2);
    }

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to the calibration editor)
    return;
}

//______________________________________________________________________________
void RunBrowser()
{
    // Show the run browser.

    // create a CaLib container
    TCContainer c("container");

    // dump all runs
    TCMySQLManager::GetManager()->DumpRuns(&c);

    // get number of runs
    Int_t nRuns = c.GetNRuns();

    // clear the screen
    clear();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "RUN BROWSER");
    attroff(A_UNDERLINE);

    // build the windows
    Int_t colLengthTot;
    WINDOW* header = 0;
    WINDOW* table = FormatRunTable(c, &colLengthTot, &header);

    // user information
    Char_t tmp[256];
    sprintf(tmp, "%d runs found. Use UP/DOWN keys to scroll "
                 "(PAGE-UP or 'p' / PAGE-DOWN or 'n' for fast mode) - hit ESC or 'q' to exit", nRuns);
    PrintStatusMessage(tmp);

    // refresh windows
    refresh();
    prefresh(header, 0, 0, 6, 2, 7, gNcol-3);
    prefresh(table, 0, 0, 7, 2, gNrow-3, gNcol-3);

    Int_t first_row = 0;
    Int_t first_col = 0;
    Int_t winHeight = gNrow-3-7;
    Int_t winWidth = gNcol;

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        //
        // decide what to do
        //

        // go up one entry
        if (c == KEY_UP)
        {
            if (first_row > 0) first_row--;
        }
        // go down one entry
        else if (c == KEY_DOWN)
        {
            if (first_row < nRuns-winHeight-1) first_row++;
        }
        // go up one page
        else if (c == KEY_PPAGE || c == 'p')
        {
            if (first_row > winHeight-1) first_row -= winHeight+1;
            else first_row = 0;
        }
        // go down one page
        else if (c == KEY_NPAGE || c == 'n')
        {
            if (first_row < nRuns-winHeight-winHeight) first_row += winHeight+1;
            else first_row = nRuns-winHeight-1;
        }
        // go right
        else if (c == KEY_RIGHT)
        {
            if (first_col < colLengthTot-winWidth) first_col += 10;
            else first_col = colLengthTot-winWidth;
        }
        // go left
        else if (c == KEY_LEFT)
        {
            if (first_col > 0) first_col -= 10;
            else continue;
        }
        // exit
        else if (c == KEY_ESC || c == 'q') break;

        // update window
        prefresh(header, 0, first_col, 6, 2, 7, gNcol-3);
        prefresh(table, first_row, first_col, 7, 2, gNrow-3, gNcol-3);
    }

    // clean-up
    delwin(table);
    delwin(header);

    // go back (to run editor)
    return;
}

//______________________________________________________________________________
void CalibBrowser(Bool_t browseTypes)
{
    // Show the calibration browser.
    // Browse calibration types if 'browseTypes' is kTRUE, otherwise browser
    // calibration data.

    // init choice
    Int_t choice = 0;

    while (kTRUE)
    {
        // select a calibration
        if ((choice = SelectCalibration(0, choice)) == -1) return;

        // select calibration data or type
        if (browseTypes)
        {
            if (SelectCalibrationType() >= 0) break;
        }
        else
        {
            if (SelectCalibrationData() >= 0) break;
        }
    }

    Bool_t redo = kTRUE;
    while (redo)
    {
        // create a CaLib container
        TCContainer c("container");

        // dump calibrations
        if (browseTypes)
            TCMySQLManager::GetManager()->DumpCalibrations(&c, gCalibration, gCalibrationType->GetData(0)->GetName());
        else
            TCMySQLManager::GetManager()->DumpCalibrations(&c, gCalibration, gCalibrationData->GetName());

        // get number of calibrations
        Int_t nCalib = c.GetNCalibrations();

        // clear the screen
        clear();

        // draw header
        DrawHeader();

        // draw title
        attron(A_UNDERLINE);
        mvprintw(4, 2, "CALIBRATION BROWSER");
        attroff(A_UNDERLINE);

        // draw calibration data or type
        if (browseTypes) mvprintw(6, 2, "Calibration type: %s", gCalibrationType->GetTitle());
        else mvprintw(6, 2, "Calibration data: %s", gCalibrationData->GetTitle());

        // build the windows
        Int_t colLengthTot;
        WINDOW* header = 0;
        WINDOW* table;
        if (browseTypes) table = FormatCalibTable(c, &colLengthTot, &header, kFALSE);
        else table = FormatCalibTable(c, &colLengthTot, &header, kTRUE);

        // user information
        Char_t tmp[256];
        sprintf(tmp, "%d sets found. Use UP/DOWN keys to scroll "
                     "(PAGE-UP or 'p' / PAGE-DOWN or 'n' for fast mode) - hit ESC or 'q' to exit", nCalib);
        PrintStatusMessage(tmp);

        // user interface geometry
        Int_t rOffset;
        if (browseTypes) rOffset = 15;
        else rOffset = 3;
        Int_t first_row = 0;
        Int_t first_col = 0;
        Int_t winHeight = gNrow-rOffset-9;
        Int_t winWidth = gNcol;

        // set operations
        if (browseTypes)
        {
            mvprintw(gNrow-13, 2, "Set operations:");
            mvprintw(gNrow-12, 2, "[s] split set    [m] merge sets");
        }

        // refresh windows
        refresh();
        prefresh(header, 0, 0, 8, 2, 9, gNcol-3);
        prefresh(table, 0, 0, 9, 2, gNrow-rOffset, gNcol-3);

        // wait for input
        for (;;)
        {
            // get key
            Int_t c = getch();

            //
            // decide what to do
            //

            // go up one entry
            if (c == KEY_UP)
            {
                if (first_row > 0) first_row--;
            }
            // go down one entry
            else if (c == KEY_DOWN)
            {
                if (first_row < nCalib-winHeight-1) first_row++;
            }
            // go up one page
            else if (c == KEY_PPAGE || c == 'p')
            {
                if (first_row > winHeight-1) first_row -= winHeight+1;
                else first_row = 0;
            }
            // go down one page
            else if (c == KEY_NPAGE || c == 'n')
            {
                if (first_row < nCalib-winHeight-winHeight) first_row += winHeight+1;
                else first_row = nCalib-winHeight-1;
            }
            // go right
            else if (c == KEY_RIGHT)
            {
                if (first_col < colLengthTot-winWidth) first_col += winWidth/2;
                else first_col = colLengthTot-winWidth;
            }
            // go left
            else if (c == KEY_LEFT)
            {
                if (first_col > 0) first_col -= (winWidth/2);
                else continue;
            }
            // split set
            else if (browseTypes && c == 's')
            {
                SplitSet();
                break;
            }
            // merge set
            else if (browseTypes && c == 'm')
            {
                MergeSets();
                break;
            }
            // exit
            else if (c == KEY_ESC || c == 'q')
            {
                redo = kFALSE;
                break;
            }

            // update window
            prefresh(header, 0, first_col, 8, 2, 9, gNcol-3);
            prefresh(table, first_row, first_col, 9, 2, gNrow-rOffset, gNcol-3);
        }

        // clean-up
        delwin(table);
        delwin(header);
    }

    // go back (to calibration editor)
    return;
}

//______________________________________________________________________________
void ChangeRunEntry(const Char_t* title, const Char_t* name, RunEntry_t entry)
{
    // Change the run entry 'entry' with name for a certain range of runs.
    // Use 'title' as the window title.

    Int_t first_run;
    Int_t last_run;
    Char_t new_value[256];
    Char_t answer[16];

    // clear the screen
    clear();

    // echo input
    echo();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, title);
    attroff(A_UNDERLINE);

    // ask first run
    mvprintw(6, 2, "Enter first run                            : ");
    scanw((Char_t*)"%d", &first_run);

    // ask last run
    mvprintw(7, 2, "Enter last run                             : ");
    scanw((Char_t*)"%d", &last_run);

    // ask new value
    mvprintw(8, 2, "Enter new %-32s : ", name);
    scanw((Char_t*)"%s", new_value);

    // ask confirmation
    mvprintw(10, 2, "Changing %s for runs %d to %d to '%s'", name, first_run, last_run, new_value);
    mvprintw(12, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(14, 2, "Aborted.");
    }
    else
    {
        Bool_t ret = kFALSE;

        // check what to change
        if (entry == kPATH)
            ret = TCMySQLManager::GetManager()->ChangeRunPath(first_run, last_run, new_value);
        else if (entry == kTARGET)
            ret = TCMySQLManager::GetManager()->ChangeRunTarget(first_run, last_run, new_value);
        else if (entry == kTARGET_POL)
            ret = TCMySQLManager::GetManager()->ChangeRunTargetPol(first_run, last_run, new_value);
        else if (entry == kTARGET_POL_DEG)
            ret = TCMySQLManager::GetManager()->ChangeRunTargetPolDeg(first_run, last_run, atof(new_value));
        else if (entry == kBEAM_POL)
            ret = TCMySQLManager::GetManager()->ChangeRunBeamPol(first_run, last_run, new_value);
        else if (entry == kBEAM_POL_DEG)
            ret = TCMySQLManager::GetManager()->ChangeRunBeamPolDeg(first_run, last_run, atof(new_value));

        // check return value
        if (ret)
        {
            mvprintw(14, 2, "%s for runs %d to %d was successfully changed to '%s'",
                            name, first_run, last_run, new_value);
        }
        else
            mvprintw(14, 2, "There was an error during %s changing for runs %d to %d!",
                            name, first_run, last_run);
    }

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to run editor)
    return;
}

//______________________________________________________________________________
Int_t SelectCalibration(const Char_t* mMsg, Int_t active)
{
    // Show the calibration selection.

    // get all calibrations
    TList* c = TCMySQLManager::GetManager()->GetAllCalibrations();

    // check if there are some calibrations
    if (!c)
    {
        strcpy(gFinishMessage, "No calibrations found!");
        Finish(-1);
    }

    // get number of calibrations
    Int_t nCalib = c->GetSize();

    // menu configuration
    const Char_t mTitle[] = "CALIBRATION SELECTION";
    const Char_t mMsgStd[] = "Select a calibration";
    const Int_t mN = nCalib;
    Char_t* mEntries[mN];
    for (Int_t i = 0; i < mN; i++) mEntries[i] = new Char_t[256];

    // fill calibrations
    TIter next(c);
    TObjString* s;
    Int_t d = 0;
    while ((s = (TObjString*)next()))
        strcpy(mEntries[d++], s->GetString().Data());

    // clean-up
    delete c;

    // show menu
    Int_t choice = 0;
    if (mMsg) choice = ShowMenu(mTitle, mMsg, mN, (Char_t**)mEntries, active);
    else choice = ShowMenu(mTitle, mMsgStd, mN, (Char_t**)mEntries, active);

    // save selected calibration
    if (choice >= 0 )
        strcpy(gCalibration, mEntries[choice]);

    // clean-up
    for (Int_t i = 0; i < mN; i++) delete [] mEntries[i];

    // return
    return choice;
}

//______________________________________________________________________________
Int_t SelectCalibrationType()
{
    // Show the calibration type selection.

    // menu configuration
    const Char_t mTitle[] = "CALIBRATION TYPE SELECTION";
    const Char_t mMsg[] = "Select a calibration type";

    // get type table
    THashList* table = TCMySQLManager::GetManager()->GetTypeTable();

    // number of types
    Int_t mN = table->GetSize();

    // create string array
    Char_t** mEntries = new Char_t*[mN];
    for (Int_t i = 0; i < mN; i++) mEntries[i] = new Char_t[256];

    // loop over types
    TIter next(table);
    TCCalibType* e;
    TCCalibType* list[mN];
    Int_t n = 0;
    while ((e = (TCCalibType*)next()))
    {
        strcpy(mEntries[n], e->GetTitle());
        list[n++] = e;
    }

    // show menu
    Int_t choice = ShowMenu(mTitle, mMsg, mN, mEntries, 0);

    // clean-up
    for (Int_t i = 0; i < mN; i++) delete [] mEntries[i];
    delete [] mEntries;

    // save selected calibration type
    if (choice >= 0)
        gCalibrationType = list[choice];

    return choice;
}

//______________________________________________________________________________
Int_t SelectCalibrationData()
{
    // Show the calibration data selection.

    // menu configuration
    const Char_t mTitle[] = "CALIBRATION DATA SELECTION";
    const Char_t mMsg[] = "Select a calibration data";

    // get type table
    THashList* table = TCMySQLManager::GetManager()->GetDataTable();

    // number of data
    Int_t mN = table->GetSize();

    // create string array
    Char_t** mEntries = new Char_t*[mN];
    for (Int_t i = 0; i < mN; i++) mEntries[i] = new Char_t[256];

    // loop over types
    TIter next(table);
    TCCalibData* e;
    TCCalibData* list[mN];
    Int_t n = 0;
    while ((e = (TCCalibData*)next()))
    {
        strcpy(mEntries[n], e->GetTitle());
        list[n++] = e;
    }

    // show menu
    Int_t choice = ShowMenu(mTitle, mMsg, mN, mEntries, 0);

    // clean-up
    for (Int_t i = 0; i < mN; i++) delete [] mEntries[i];
    delete [] mEntries;

    // save selected calibration type
    if (choice >= 0)
        gCalibrationData = list[choice];

    return choice;
}

//______________________________________________________________________________
void ExportRuns()
{
    // Export runs.

    Int_t first_run = 0;
    Int_t last_run = 0;
    Char_t filename[256];
    Char_t answer[16];

    // clear the screen
    clear();

    // echo input
    echo();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "EXPORT RUNS");
    attroff(A_UNDERLINE);

    // ask first run, last run and file name
    mvprintw(6, 2, "First run (enter 0 to select all)          : ");
    scanw((Char_t*)"%d", &first_run);
    mvprintw(7, 2, "Last run  (enter 0 to select all)          : ");
    scanw((Char_t*)"%d", &last_run);
    mvprintw(8, 2, "Name of output file                        : ");
    scanw((Char_t*)"%s", filename);

    // ask confirmation
    if (!first_run && !last_run)
        mvprintw(10, 2, "Saving all runs to '%s'", filename);
    else
        mvprintw(10, 2, "Saving runs %d to %d to '%s'", first_run, last_run, filename);
    mvprintw(12, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(14, 2, "Aborted.");
    }
    else
    {
        // create new container
        TCContainer* container = new TCContainer(TCConfig::kCaLibDumpName);

        // dump runs to container
        Int_t nRun = TCMySQLManager::GetManager()->DumpRuns(container, first_run, last_run);

        // save container to ROOT file
        Bool_t ret = container->Save(filename, kTRUE);

        // clean-up
        delete container;

        // check return value
        if (ret)
            mvprintw(14, 2, "Saved %d runs to '%s'", nRun, filename);
        else
            mvprintw(14, 2, "Could not save %d runs to '%s'!", nRun, filename);
    }

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to admin menue)
    return;
}

//______________________________________________________________________________
void ExportCalibration()
{
    // Export calibration.

    Char_t filename[256];
    Char_t answer[16];

    // select a calibration
    if (SelectCalibration() == -1) return;

    // clear the screen
    clear();

    // echo input
    echo();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "EXPORT CALIBRATION");
    attroff(A_UNDERLINE);

    // ask calibration and output file
    mvprintw(6, 2, "Calibration to export                      : %s", gCalibration);
    mvprintw(7, 2, "Name of output file                        : ");
    scanw((Char_t*)"%s", filename);

    // ask confirmation
    mvprintw(9, 2, "Saving the calibration '%s' to '%s'", gCalibration, filename);
    mvprintw(11, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(13, 2, "Aborted.");
    }
    else
    {
        // create new container
        TCContainer* container = new TCContainer(TCConfig::kCaLibDumpName);

        // dump calibration to container
        Int_t nCalib = TCMySQLManager::GetManager()->DumpAllCalibrations(container, gCalibration);

        // save container to ROOT file
        Bool_t ret = container->Save(filename, kTRUE);

        // clean-up
        delete container;

        // check return value
        if (ret)
            mvprintw(13, 2, "Saved %d calibrations to '%s'", nCalib, filename);
        else
            mvprintw(13, 2, "Could not save %d calibrations to '%s'!", nCalib, filename);
    }

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to admin menue)
    return;
}

//______________________________________________________________________________
void ImportRuns()
{
    // Import runs.

    Char_t filename[256];
    Char_t answer[16];

    // clear the screen
    clear();

    // echo input
    echo();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "IMPORT RUNS");
    attroff(A_UNDERLINE);

    // ask file name
    mvprintw(6, 2, "Name of input file                         : ");
    scanw((Char_t*)"%s", filename);

    // ask confirmation
    mvprintw(8, 2, "Importing all runs from '%s'", filename);
    mvprintw(10, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(12, 2, "Aborted.");
    }
    else
    {
        // try to load the container
        TCContainer* c = TCMySQLManager::GetManager()->LoadContainer(filename);
        if (!c)
        {
            mvprintw(12, 2, "No CaLib container found in file '%s'!", filename);
            goto error_import_run;
        }

        // import runs
        Int_t nRun = TCMySQLManager::GetManager()->ImportRuns(c);

        // check return value
        if (nRun)
            mvprintw(12, 2, "Imported %d runs from '%s'", nRun, filename);
        else
            mvprintw(12, 2, "No run was imported from '%s'!", filename);

        // clean-up
        delete c;
    }

    // ugly goto label
    error_import_run:

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to admin menue)
    return;
}

//______________________________________________________________________________
void ImportCalibration()
{
    // Import calibration.

    Char_t tmp[256];
    Char_t filename[256];
    Char_t answer[16];
    Bool_t importAll = kFALSE;
    Bool_t calibExists = kFALSE;
    Int_t nCalib;
    const Char_t* calibName;

    // menu configuration
    const Char_t mTitle[] = "IMPORT CALIBRATION";
    const Char_t mMsg[] = "Configure importing";
    const Int_t mN = 3;
    const Char_t* mEntries[] = { "Import complete calibration",
                                 "Import single calibration",
                                 "Go back" };

    // show menu
    Int_t choice = ShowMenu(mTitle, mMsg, mN, (Char_t**)mEntries, 0);

    // decide what do to
    switch (choice)
    {
        case -1:
            return;
        case 0:
            importAll = kTRUE;
            break;
        case 1:
            if (SelectCalibrationType() == -1) return;
            sprintf(tmp, "Select the calibration the imported '%s' is added to",
                    gCalibrationType->GetTitle());
            if (SelectCalibration(tmp) == -1) return;
            break;
        case 2:
            return;
    }

    // clear the screen
    clear();

    // echo input
    echo();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "IMPORT CALIBRATION");
    attroff(A_UNDERLINE);

    // ask file name
    mvprintw(6, 2, "Name of input file                         : ");
    scanw((Char_t*)"%s", filename);

    // try to load the container
    TCContainer* c = TCMySQLManager::GetManager()->LoadContainer(filename);
    if (!c)
    {
        mvprintw(8, 2, "No CaLib container found in file '%s'!", filename);
        goto error_import_calib;
    }

    // get number of calibrations
    nCalib = c->GetNCalibrations();
    if (!nCalib)
    {
        mvprintw(8, 2, "No calibrations found in file '%s'!", filename);
        goto error_import_calib;
    }

    // check if calibration exists
    calibName = c->GetCalibration(0)->GetCalibration();
    if (TCMySQLManager::GetManager()->ContainsCalibration(calibName))
        calibExists = kTRUE;

    // ask confirmation
    if (importAll)
    {
        // calibration exits or not
        if (calibExists)
        {
            mvprintw(8, 2, "Calibration '%s' exists already in database - it will be overwritten!",
                     calibName);
        }
        else
        {
            mvprintw(8, 2, "Importing %d calibrations of '%s' from '%s'",
                     nCalib, calibName, filename);
        }
    }
    else
    {
        mvprintw(8, 2, "'%s' from '%s' will be added to calibration '%s' replacing the existing data!",
                 gCalibrationType->GetTitle(), filename, gCalibration);
    }
    mvprintw(10, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(12, 2, "Aborted.");
    }
    else
    {
        // all calibrations
        if (importAll)
        {
            // delete old calibration before importing
            if (calibExists) TCMySQLManager::GetManager()->RemoveAllCalibrations(calibName);

            // import new calibration
            Int_t nAdded = TCMySQLManager::GetManager()->ImportCalibrations(c);

            // check return value
            if (nAdded)
                mvprintw(12, 2, "Imported %d calibrations of '%s' from '%s'", nAdded, calibName, filename);
            else
                mvprintw(12, 2, "No calibration was imported from '%s'!", filename);
        }
        else
        {
            Int_t nAdded = 0;

            // loop over calibration data for this calibration type
            TIter next(gCalibrationType->GetData());
            TCCalibData* d;
            while ((d = (TCCalibData*)next()))
            {
                // delete old calibration before importing
                TCMySQLManager::GetManager()->RemoveCalibration(gCalibration, d->GetName());

                // import new calibration
                nAdded += TCMySQLManager::GetManager()->ImportCalibrations(c, gCalibration, d->GetName());
            }

            // check return value
            if (nAdded)
                mvprintw(12, 2, "Imported %d calibrations of '%s' from '%s' and added them to '%s'",
                         nAdded, gCalibrationType->GetTitle(), filename, gCalibration);
            else
                mvprintw(12, 2, "No calibration was imported from '%s'!", filename);
        }

        // clean-up
        delete c;
    }

    // ugly goto label
    error_import_calib:

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to admin menue)
    return;
}

//______________________________________________________________________________
void CloneCalibration()
{
    // Clone calibration.

    Char_t calib[256];
    Char_t desc[256];
    Char_t answer[16];
    Int_t firstRun, lastRun;

    // select a calibration
    if (SelectCalibration() == -1) return;

    // clear the screen
    clear();

    // echo input
    echo();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "CLONE CALIBRATION");
    attroff(A_UNDERLINE);

    // ask new calibration, description, first and last runs
    mvprintw(6, 2, "Calibration to export                      : %s", gCalibration);
    mvprintw(7, 2, "Name of cloned calibration                 : ");
    scanw((Char_t*)"%s", calib);
    mvprintw(8, 2, "Description of cloned calibration          : ");
    scanw((Char_t*)"%s", desc);
    mvprintw(9, 2, "First run of cloned calibration            : ");
    scanw((Char_t*)"%d", &firstRun);
    mvprintw(10, 2, "Last run of cloned calibration             : ");
    scanw((Char_t*)"%d", &lastRun);

    // ask confirmation
    mvprintw(12, 2, "Cloning the calibration '%s' to '%s' for runs [%d,%d]",
             gCalibration, calib, firstRun, lastRun);
    mvprintw(14, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(16, 2, "Aborted.");
    }
    else
    {
        // clone the calibration
        Bool_t ret = TCMySQLManager::GetManager()->CloneCalibration(gCalibration, calib, desc,
                                                                    firstRun, lastRun);

        // check return value
        if (ret)
            mvprintw(16, 2, "Cloned the calibration '%s' to '%s'", gCalibration, calib);
        else
            mvprintw(16, 2, "Could not clone the calibration '%s' to '%s!", gCalibration, calib);
    }

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to admin menue)
    return;
}

//______________________________________________________________________________
void ExportDatabase()
{
    // Export the complete database.

    Char_t filename[256];
    Char_t answer[16];

    // clear the screen
    clear();

    // echo input
    echo();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "EXPORT DATABASE");
    attroff(A_UNDERLINE);

    // ask file name
    mvprintw(6, 2, "Name of output file                         : ");
    scanw((Char_t*)"%s", filename);

    // ask confirmation
    mvprintw(8, 2, "Exporting the complete database to '%s'", filename);
    mvprintw(10, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(12, 2, "Aborted.");
    }
    else
    {
        // try to export
        if (TCMySQLManager::GetManager()->ExportDatabase(filename))
            mvprintw(12, 2, "Exported the complete database to '%s'", filename);
        else
            mvprintw(12, 2, "Could not export the complete database to '%s'!", filename);
    }

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to admin menue)
    return;
}

//______________________________________________________________________________
void RenameCalibration()
{
    // Rename a calibration.

    Char_t newName[256];
    Char_t answer[16];

    // select a calibration
    if (SelectCalibration() == -1) return;

    // clear the screen
    clear();

    // echo input
    echo();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "RENAME CALIBRATION");
    attroff(A_UNDERLINE);

    // ask new name
    mvprintw(6, 2, "Old calibration identifier                 : %s", gCalibration);
    mvprintw(7, 2, "Enter new calibration identifier           : ");
    scanw((Char_t*)"%s", newName);

    // ask confirmation
    mvprintw(9, 2, "Renaming calibration '%s' to '%s'", gCalibration, newName);
    mvprintw(11, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(13, 2, "Aborted.");
    }
    else
    {
        // rename calibration
        Bool_t ret = TCMySQLManager::GetManager()->ChangeCalibrationName(gCalibration, newName);

        // check return value
        if (ret)
        {
            mvprintw(13, 2, "Renamed calibration '%s' to '%s'",
                            gCalibration, newName);
        }
        else
            mvprintw(13, 2, "There was an error during renaming the calibration '%s' to '%s'!",
                            gCalibration, newName);
    }

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to calib editor)
    return;
}

//______________________________________________________________________________
void ChangeRunRange()
{
    // Change the run range of a calibration.

    Char_t answer[16];
    Int_t newFirstRun;
    Int_t newLastRun;

    // select a calibration
    if (SelectCalibration() == -1) return;

    // clear the screen
    clear();

    // echo input
    echo();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "CHANGE RUN RANGE OF CALIBRATION");
    attroff(A_UNDERLINE);

    // ask new name
    mvprintw(6, 2, "Selected calibration                                       : %s", gCalibration);
    mvprintw(7, 2, "New first run (run has to be in DB, 0 to keep current one) : ");
    scanw((Char_t*)"%d", &newFirstRun);
    mvprintw(8, 2, "New last run (run has to be in DB, 0 to keep current one)  : ");
    scanw((Char_t*)"%d", &newLastRun);

    // ask confirmation
    mvprintw(10, 2, "Changing run range to [%d,%d]", newFirstRun, newLastRun);
    mvprintw(12, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(14, 2, "Aborted.");
    }
    else
    {
        // rename calibration
        Bool_t ret = TCMySQLManager::GetManager()->ChangeCalibrationRunRange(gCalibration, newFirstRun, newLastRun);

        // check return value
        if (ret)
        {
            mvprintw(14, 2, "Changed run range of calibration '%s' to [%d,%d]",
                            gCalibration, newFirstRun, newLastRun);
        }
        else
            mvprintw(14, 2, "There was an error during changing the run range!");
    }

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to calib editor)
    return;
}

//______________________________________________________________________________
void ChangeDescription()
{
    // Change the description of a calibration.

    Char_t newDesc[256];
    Char_t answer[16];

    // select a calibration
    if (SelectCalibration() == -1) return;

    // clear the screen
    clear();

    // echo input
    echo();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "CHANGE DESCRIPTION OF CALIBRATION");
    attroff(A_UNDERLINE);

    // ask new description
    mvprintw(6, 2, "Selected calibration                       : %s", gCalibration);
    mvprintw(7, 2, "Enter new description                      : ");
    scanw((Char_t*)"%[^\n]s", newDesc);

    // ask confirmation
    mvprintw(9, 2, "Changing description of calibration '%s' to '%s'", gCalibration, newDesc);
    mvprintw(11, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(13, 2, "Aborted.");
    }
    else
    {
        // change description
        Bool_t ret = TCMySQLManager::GetManager()->ChangeCalibrationDescription(gCalibration, newDesc);

        // check return value
        if (ret)
        {
            mvprintw(13, 2, "Changed description of calibration '%s' to '%s'",
                            gCalibration, newDesc);
        }
        else
            mvprintw(13, 2, "There was an error during changing the description of the calibration '%s' to '%s'!",
                            gCalibration, newDesc);
    }

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to calib editor)
    return;
}

//______________________________________________________________________________
void DeleteCalibration()
{
    // Delete a calibration.

    Char_t answer[16];

    // select a calibration
    if (SelectCalibration() == -1) return;

    // clear the screen
    clear();

    // echo input
    echo();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "DELETE CALIBRATION");
    attroff(A_UNDERLINE);

    // ask confirmation
    mvprintw(6, 2, "Deleting calibration '%s'", gCalibration);
    mvprintw(8, 6, "Are you sure to continue? (yes/no) : ");
    scanw((Char_t*)"%s", answer);
    if (strcmp(answer, "yes"))
    {
        mvprintw(10, 2, "Aborted.");
    }
    else
    {
        // delete calibration
        Bool_t ret = TCMySQLManager::GetManager()->RemoveAllCalibrations(gCalibration);

        // check return value
        if (ret)
            mvprintw(10, 2, "Deleted calibration '%s'", gCalibration);
        else
            mvprintw(10, 2, "There was an error during deleting the calibration '%s'!", gCalibration);
    }

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // don't echo input
    noecho();

    // go back (to calib editor)
    return;
}

//______________________________________________________________________________
void RunEditor()
{
    // Show the run editor.

    // menu configuration
    const Char_t mTitle[] = "RUN EDITOR";
    const Char_t mMsg[] = "Select a run operation";
    const Int_t mN = 8;
    const Char_t* mEntries[] = { "Browse runs",
                                 "Change path",
                                 "Change target",
                                 "Change target polarization",
                                 "Change degree of target polarization",
                                 "Change beam polarization",
                                 "Change degree of beam polarization",
                                 "Go back" };
    // menue index
    Int_t choice = 0;

    // menue loop
    while (kTRUE)
    {
        // show menu
        choice = ShowMenu(mTitle, mMsg, mN, (Char_t**)mEntries, choice);

        // decide what do to
        switch (choice)
        {
            case -1: return;
            case  0: RunBrowser();
                     break;
            case  1: ChangeRunEntry("CHANGE PATH", "path", kPATH);
                     break;
            case  2: ChangeRunEntry("CHANGE TARGET", "target", kTARGET);
                     break;
            case  3: ChangeRunEntry("CHANGE TARGET POLARIZATION",
                                    "target polarization", kTARGET_POL);
                     break;
            case  4: ChangeRunEntry("CHANGE DEGREE OF TARGET POLARIZATION",
                                    "degree of target polarization", kTARGET_POL_DEG);
                     break;
            case  5: ChangeRunEntry("CHANGE BEAM POLARIZATION",
                                    "beam polarization", kBEAM_POL);
                     break;
            case  6: ChangeRunEntry("CHANGE DEGREE OF BEAM POLARIZATION",
                                    "degree of beam polarization", kBEAM_POL_DEG);
                     break;
            case  7: return;
        }
    }
}

//______________________________________________________________________________
void CalibEditor()
{
    // Show the calibration editor.

    // menu configuration
    const Char_t mTitle[] = "CALIBRATION EDITOR";
    const Char_t mMsg[] = "Select a calibration operation";
    const Int_t mN = 7;
    const Char_t* mEntries[] = { "Browse calibration data",
                                 "Manipulate calibration sets",
                                 "Change run range",
                                 "Change description",
                                 "Rename calibration",
                                 "Delete calibration",
                                 "Go back" };

    // menue index
    Int_t choice = 0;

    // menue loop
    while (kTRUE)
    {
        // show menu
        choice = ShowMenu(mTitle, mMsg, mN, (Char_t**)mEntries, choice);

        // decide what do to
        switch (choice)
        {
            case -1: return;
            case  0: CalibBrowser(kFALSE);
                     break;
            case  1: CalibBrowser(kTRUE);
                     break;
            case  2: ChangeRunRange();
                     break;
            case  3: ChangeDescription();
                     break;
            case  4: RenameCalibration();
                     break;
            case  5: DeleteCalibration();
                     break;
            case  6: return;
        }
    }
}

//______________________________________________________________________________
void Administration()
{
    // Show the administration menu.

    // menu configuration
    const Char_t mTitle[] = "ADMINISTRATION";
    const Char_t mMsg[] = "Select an administration operation";
    const Int_t mN = 7;
    const Char_t* mEntries[] = { "Export runs",
                                 "Export calibration",
                                 "Import runs",
                                 "Import calibration",
                                 "Clone calibration",
                                 "Export complete database",
                                 "Go back" };

    // menue index
    Int_t choice = 0;

    // menue loop
    while (kTRUE)
    {
        // show menu
        choice = ShowMenu(mTitle, mMsg, mN, (Char_t**)mEntries, choice);

        // decide what do to
        switch (choice)
        {
            case -1: return;
            case  0: ExportRuns();
                     break;
            case  1: ExportCalibration();
                     break;
            case  2: ImportRuns();
                     break;
            case  3: ImportCalibration();
                     break;
            case  4: CloneCalibration();
                     break;
            case  5: ExportDatabase();
                     break;
            case  6: return;
        }
    }
}

//______________________________________________________________________________
void About()
{
    // Show the about screen.

    // clear the screen
    clear();

    // draw header
    DrawHeader();

    // draw title
    attron(A_UNDERLINE);
    mvprintw(4, 2, "ABOUT");
    attroff(A_UNDERLINE);

    // print some information
    mvprintw(6,  2, "CaLib Manager");
    mvprintw(7,  2, "an ncurses-based CaLib administration software");
    mvprintw(8,  2, "(c) 2011-2017 by Dominik Werthmueller, Thomas Strub");
    mvprintw(10, 2, "CaLib - A2 calibration system");
    mvprintw(11, 2, "Version %s", TCConfig::kCaLibVersion);
    mvprintw(12, 2, "(c) 2010-2017 by Dominik Werthmueller,");
    mvprintw(13, 2, "                 Irakli Keshelashvili,");
    mvprintw(14, 2, "                 Thomas Strub,");
    mvprintw(15, 2, "                 University of Basel");

    // user information
    PrintStatusMessage("Hit ESC or 'q' to exit");

    // wait for input
    for (;;)
    {
        // get key
        Int_t c = getch();

        // leave loop
        if (c == KEY_ESC || c == 'q') break;
    }

    // go back (to the main menu)
    return;
}

//______________________________________________________________________________
void MainMenu()
{
    // Show the main menu.

    // menu configuration
    const Char_t mTitle[] = "MAIN MENU";
    const Char_t mMsg[] = "Select an operation";
    const Int_t mN = 5;
    const Char_t* mEntries[] = { "Run editor",
                                 "Calibration editor",
                                 "Administration",
                                 "About",
                                 "Exit" };

    // menue index
    Int_t choice = 0;

    // main menue loop
    while (kTRUE)
    {
        // show menu
        choice = ShowMenu(mTitle, mMsg, mN, (Char_t**)mEntries, choice);

        // decide what do to
        switch (choice)
        {
            case -1: choice = 0;
                     break;
            case  0: RunEditor();
                     break;
            case  1: CalibEditor();
                     break;
            case  2: Administration();
                     break;
            case  3: About();
                     break;
            case  4: return;
        }
    }
}

//______________________________________________________________________________
Int_t main(Int_t argc, Char_t* argv[])
{
    // Main method.

    // init finish message
    gFinishMessage[0] = '\0';

    // set-up signal for CTRL-C
    signal(SIGINT, Finish);

    // redirect standard error
    FILE* f = freopen("calib_manager.log", "w", stderr);
    if (!f)
    {
        strcpy(gFinishMessage, "Cannot write log file to current directory!");
        Finish(-1);
    }

    // init ncurses
    initscr();

    // enable keyboard mapping (non-BSD feature)
    keypad(stdscr, TRUE);

    // tell curses not to do NL->CR/NL on output
    nonl();

    // take input chars one at a time, don't wait for <\\>n
    cbreak();

    // don't echo input
    noecho();

    // get number of rows and columns
    getmaxyx(stdscr, gNrow, gNcol);

    // check connection to database
    if (!TCMySQLManager::GetManager())
    {
        strcpy(gFinishMessage, "No connection to CaLib database!");
        Finish(-1);
    }

    // check dimensions
    if (gNrow < 40 || gNcol < 120)
    {
        sprintf(gFinishMessage, "Cannot run in a terminal smaller than 40 rows and 120 columns!");
        Finish(-1);
    }

    // set MySQL manager to silence mode
    //TCMySQLManager::GetManager()->SetSilenceMode(kTRUE);

    // show main menu
    MainMenu();

    // exit program
    Finish(0);
}

