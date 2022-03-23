/* Additional info:
 *
 * One needs to set up the 2020 experiments s444 and s467, unpackers are:
 *
 * at $UCESB_DIR/../upexps/202002_s444 and $UCESB_DIR/../upexps/202002_s467
 *
 * Before executing the macro, set up the paths to lmd files, upexps and
 * califa mapping file "califamapfilename". Also set up the name of output
 * file "outputCalFile"
 *
 * @since January 18th, 2020
 */

typedef struct EXT_STR_h101_t
{
    EXT_STR_h101_unpack_t unpack;
    EXT_STR_h101_CALIFA_t califa;
} EXT_STR_h101;

void califa_calibParFinder()
{
    TStopwatch timer;
    timer.Start();

    const Int_t nev = -1; /* number of events to read, -1 - until CTRL+C */

    /* Create source using ucesb for input ------------------ */
    // TString filename = "--stream=lxg0898:6002";
    //TString filename = "/d/land5/202105_s494/califa/lmd_leyla/main0671_*.lmd";// path to 22Na lmd files for s494 June 2021
    //TString filename = "/d/land4/202105_s494/lmd/main0676_*.lmd";// path to 22Na lmd files for s494 June 2021
    TString filename = "/lustre/land/202002_s467/stitched/main0383_0001.lmd";// path to 22Na lmd files for s467 
    TString outputFileName = "calib_latar_s467_22Na_root.root";

    /* Create source using ucesb for input ------------------ */
    // UCESB paths
    TString ntuple_options = "RAW";
    TString ucesb_dir = getenv("UCESB_DIR");
    //TString upexps_dir = "/u/land/fake_cvmfs/9.13/upexps/20204_s467/";
    TString upexps_dir = "/u/land/latar/";

    TString ucesb_path;
    ucesb_path = upexps_dir + "upexps/202002_s467_jentob/202002_s467 --allow-errors --input-buffer=100Mi";
    ucesb_path.ReplaceAll("//", "/");

    // Parameters for CALIFA
    TString dir = gSystem->Getenv("VMCWORKDIR");
    //TString califamapdir = dir + "/sofia/macros/s444/parameters/";
    TString califamapfilename = "/u/land/latar/s467/R3BRoot/macros/r3b/unpack/s467/califa/parameters/CALIFA_mapping.par";//wrong
    //TString califamapfilename = "CALIFA_mapping_s467_Nov2021.par";//new mapping for s467 November 2021
    califamapfilename.ReplaceAll("//", "/");

    // CALIFA output file with parameters for calibrating in keV
    TString outputCalFile = "calibpar_latar_s467_22Na_root.root";

    /* Definition of reader --------------------------------- */
    EXT_STR_h101 ucesb_struct;

    R3BUcesbSource* source = new R3BUcesbSource(filename, ntuple_options, ucesb_path, &ucesb_struct, sizeof(ucesb_struct));
    source->SetMaxEvents(nev);

    source->AddReader(new R3BUnpackReader((EXT_STR_h101_unpack*)&ucesb_struct, offsetof(EXT_STR_h101, unpack)));
    source->AddReader(new R3BCalifaFebexReader((EXT_STR_h101_CALIFA*)&ucesb_struct.califa, offsetof(EXT_STR_h101, califa)));

    /* Create online run ------------------------------------ */
    FairRunOnline* run = new FairRunOnline(source);
    run->SetRunId(1);
    run->SetSink(new FairRootFileSink(outputFileName));

    /* Runtime data base ------------------------------------ */
     FairRuntimeDb* rtdb = run->GetRuntimeDb();
	    // CALIFA mapping
      FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo(); // Ascii
	             parIo1->open(califamapfilename, "in");
	                 rtdb->setFirstInput(parIo1);
	                     rtdb->print();
	    
 
 
    /* Add analysis task ------------------------------------ */

    // R3BCalifaMapped2CrystalCalPar ----
    TArrayF* EnergythePeaks = new TArrayF();
    Float_t e1 = 1274.5;//511.0; 22Na
    Float_t e2 = 511.0;//1274.5; 22Na
    //Float_t e1 = 1173.2;// 60Co
    //Float_t e2 = 1332.5;// 60Co
    EnergythePeaks->Set(2);
    EnergythePeaks->AddAt(e1, 0);
    EnergythePeaks->AddAt(e2, 1);

    R3BCalifaMapped2CrystalCalPar* CalPar = new R3BCalifaMapped2CrystalCalPar();
    CalPar->SetMinStadistics(1000);
    CalPar->SetNumParameterFit(2); // OPTIONAL by default 2
    // Gamma range
    CalPar->SetCalRange_left(380);
    CalPar->SetCalRange_right(1600);
    CalPar->SetCalRange_bins(122);
    /*CalPar->SetCalRange_left(750);
    CalPar->SetCalRange_right(1500);
    CalPar->SetCalRange_bins(150);*/
    // particle range
    CalPar->SetCalRangeP_left(38);
    CalPar->SetCalRangeP_right(160);
    CalPar->SetCalRangeP_bins(122);
    CalPar->SetSigma(3.0);
    CalPar->SetThreshold(0.0001);
    CalPar->SetEnergyPeaks(EnergythePeaks);
    CalPar->SetDebugMode(1);
    run->AddTask(CalPar);

    /* Runtime data base ------------------------------------ */
    //FairRuntimeDb* rtdb = run->GetRuntimeDb();
    // CALIFA mapping
    //FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo(); // Ascii
    //parIo1->open(califamapfilename, "in");
    //rtdb->setFirstInput(parIo1);
    //rtdb->print();

    /* Initialize ------------------------------------------- */
    run->Init();
    //    FairLogger::GetLogger()->SetLogScreenLevel("WARNING");
    //    FairLogger::GetLogger()->SetLogScreenLevel("DEBUG");
    FairLogger::GetLogger()->SetLogScreenLevel("INFO");

    // Choose Root or Ascii file
    // 1-Root file with the Calibration Parameters
    Bool_t kParameterMerged = kTRUE;
    FairParRootFileIo* parOut = new FairParRootFileIo(kParameterMerged);
    parOut->open(outputCalFile);
    rtdb->setOutput(parOut);
    // 2-Ascii file with the Calibration Parameters
    /*FairParAsciiFileIo* parout = new FairParAsciiFileIo();
    parout->open("calibpar_jluis_s494_60Co_ascii.par","out");
    rtdb->setOutput(parout);*/
      

    /* Run -------------------------------------------------- */
    run->Run((nev < 0) ? nev : 0, (nev < 0) ? 0 : nev);

    /* Save parameters (if needed) -------------------------- */
    rtdb->saveOutput();

    /* Finish ----------------------------------------------- */
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    cout << endl << endl;
    cout << "Macro finished succesfully." << endl;
    cout << "Output file is " << outputFileName << endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl << endl;
    gApplication->Terminate();
}
