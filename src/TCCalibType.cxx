/*************************************************************************
 * Author: Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibType                                                          //
//                                                                      //
// Calibration type class.                                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TList.h"

#include "TCCalibType.h"
#include "TCCalibData.h"

ClassImp(TCCalibType)

//______________________________________________________________________________
TCCalibType::TCCalibType(const Char_t* name, const Char_t* title)
    : TNamed(name, title)
{
    // Constructor.

    fData = new TList();
}

//______________________________________________________________________________
TCCalibType::~TCCalibType()
{
    // Destructor.

    //if (fData) delete fData; // crashes ROOT 6.12/04
}

//______________________________________________________________________________
TCCalibData* TCCalibType::GetData(Int_t n)
{
    // Return calibration data at index 'n'.

    return fData ? (TCCalibData*)fData->At(n) : 0;
}

//______________________________________________________________________________
void TCCalibType::AddData(TCCalibData* data)
{
    // Add calibration data 'data'.

    if (fData) fData->Add(data);
}

//______________________________________________________________________________
void TCCalibType::Print(Option_t* option) const
{
    // Print the content of this class.

    printf("CaLib Type Information\n");
    printf("Name           : %s\n", GetName());
    printf("Title          : %s\n", GetTitle());
    printf("CaLib data     : ");
    TIter next(fData);
    TCCalibData* d;
    while ((d = (TCCalibData*)next()))
    {
        printf("%s,", d->GetName());
    }
    printf("\r\n");
}

