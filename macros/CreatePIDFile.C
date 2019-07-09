#include "TOLoader.h"

void CreatePIDFile()
{
    // configuration
    const Int_t n          = 3;
    const Char_t* files[n] = { "/scratch/werthm/ARHistograms_single_part_electron_0001.root",
                               "/scratch/werthm/ARHistograms_single_part_pip_0001.root",
                               "/scratch/werthm/ARHistograms_single_part_proton_0001.root" };
    const Double_t mix[n]  = { 1, 1, 1 };

    // process data
    TH3* h3[n];
    TH1* h2[n];
    TH1* hSum = 0;
    for (Int_t i = 0; i < n; i++)
    {
        // load histogram
        TOLoader::LoadObject(files[i], "CaLib_PID_dE_E_000", &h3[i]);

        // create projection
        h2[i] = h3[i]->Project3D(TString::Format("Proj_%d_yx", i).Data());;

        // scale
        h2[i]->Scale(mix[i]);

        // add projections
        if (hSum)
            hSum->Add(h2[i]);
        else
            hSum = h2[i];
    }

    hSum->SetName("PID_MC");
    hSum->Draw("colz");
    hSum->SaveAs("PID_MC.root");
}

