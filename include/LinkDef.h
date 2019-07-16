/*************************************************************************
 * Author: Irakli Keshelashvili, Dominik Werthmueller, Thomas Strub
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// LinkDef.h                                                            //
//                                                                      //
// CaLib dictionary header file.                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifdef __CINT__

// turn everything off
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all typedef;

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedef;

// enums
#pragma link C++ enum ECalibDetector;
#pragma link C++ enum ERawFileType;
#pragma link C++ enum ERawFileFormat;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  #pragma link C++ enum EWalkCorrType;
#endif

// typedefs
#pragma link C++ typedef CalibDetector_t;
#pragma link C++ typedef RawFileType_t;
#pragma link C++ typedef RawFileFormat_t;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  #pragma link C++ typedef WalkCorrType_t;
#endif

// common classes
#pragma link C++ namespace TCConfig;
#pragma link C++ namespace TCUtils;
#pragma link C++ namespace TCFitUtils;
#pragma link C++ class TCFileManager+;
#pragma link C++ class TCReadConfig+;
#pragma link C++ class TCConfigElement+;
#pragma link C++ class TCReadARCalib+;
#pragma link C++ class TCWriteARCalib+;
#pragma link C++ class TCARElement+;
#pragma link C++ class TCLine+;
#pragma link C++ class TCARTimeWalk+;
#pragma link C++ class TCARNeighbours+;
#pragma link C++ class TCReadACQU+;
#pragma link C++ class TCACQUFile+;
#pragma link C++ class TCMySQLManager+;
#pragma link C++ class TCContainer+;
#pragma link C++ class TCRun+;
#pragma link C++ class TCCalibration+;
#pragma link C++ class TCCalibData+;
#pragma link C++ class TCCalibType+;
#pragma link C++ class TCCalib+;
#pragma link C++ class TCCalibPed+;
#pragma link C++ class TCCalibDiscrThr+;
#pragma link C++ class TCCalibTime+;
#pragma link C++ class TCCalibEnergy+;
#pragma link C++ class TCCalibQuadEnergy+;
#pragma link C++ class TCCalibPhi+;
#pragma link C++ class TCCalibDeltaETrad+;
#pragma link C++ class TCCalibDroop+;
#pragma link C++ class TCCalibPeakFit+;
#pragma link C++ class TCCalibPeakFitCB+;
#pragma link C++ class TCCalibPeakFitTAPS+;

// misc calibration classes
#pragma link C++ class TCCalibTargetPosition+;

// Tagger calibration classes
#pragma link C++ class TCCalibTaggerTime+;

// CB calibration classes
#pragma link C++ class TCCalibCBEnergy+;
#pragma link C++ class TCCalibCBQuadEnergy+;
#pragma link C++ class TCCalibCBTime+;
#pragma link C++ class TCCalibCBRiseTime+;
#pragma link C++ class TCCalibCBTimeWalk+;
#pragma link C++ class TCCalibCBLED+;

// TAPS calibration classes
#pragma link C++ class TCCalibTAPSEnergyLG+;
#pragma link C++ class TCCalibTAPSEnergySG+;
#pragma link C++ class TCCalibTAPSPedLG+;
#pragma link C++ class TCCalibTAPSPedSG+;
#pragma link C++ class TCCalibTAPSPedVeto+;
#pragma link C++ class TCCalibTAPSQuadEnergy+;
#pragma link C++ class TCCalibTAPSTime+;
#pragma link C++ class TCCalibTAPSLED1+;
#pragma link C++ class TCCalibTAPSLED2+;
#pragma link C++ class TCCalibTAPSCFD+;
#pragma link C++ class TCCalibTAPSPSA+;

// PID calibration classes
#pragma link C++ class TCCalibPIDPhi+;
#pragma link C++ class TCCalibPIDDroop+;
#pragma link C++ class TCCalibPIDEnergy+;
#pragma link C++ class TCCalibPIDEnergyTrad+;
#pragma link C++ class TCCalibPIDTime+;

// Veto calibration classes
#pragma link C++ class TCCalibVetoCorr+;
#pragma link C++ class TCCalibVetoEnergy+;
#pragma link C++ class TCCalibVetoEnergyTrad+;
#pragma link C++ class TCCalibVetoTime+;
#pragma link C++ class TCCalibVetoLED+;

// Pizza detector calibration classes
#pragma link C++ class TCCalibPizzaPhi+;
#pragma link C++ class TCCalibPizzaEnergyTrad+;
#pragma link C++ class TCCalibPizzaDroop+;
#pragma link C++ class TCCalibPizzaTime+;

// Run calibration classes
#pragma link C++ class TCARFileLoader+;
#pragma link C++ class TCARHistoLoader+;
#pragma link C++ class TCBadElement+;
#pragma link C++ class TCBadScRElement+;
#pragma link C++ class TCCalibRun+;
#pragma link C++ class TCCalibRunBadScR+;
#pragma link C++ class TCCalibRunBadScR_NaI+;
#pragma link C++ class TCCalibRunBadScR_PID+;
#pragma link C++ class TCCalibRunBadScR_MWPC+;
#pragma link C++ class TCCalibRunBadScR_BaF2PWO+;
#pragma link C++ class TCCalibRunBadScR_BaF2+;
#pragma link C++ class TCCalibRunBadScR_PWO+;
#pragma link C++ class TCCalibRunBadScR_Veto+;
#pragma link C++ class TCCalibRunBadScR_Pizza+;
#pragma link C++ class TCCalibRunBadScR_Ladder+;
#pragma link C++ class TCCalibRunBadScR_LadderScalers+;
#pragma link C++ class TCCalibRunBadScR_TimeShift+;
#pragma link C++ class TCCalibRunBadScR_TimeShiftTagger+;

#endif

