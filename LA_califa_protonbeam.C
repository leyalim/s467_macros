// root -l -b -q 'LA_califa_protonbeam.C(23)'; root -l run23_histos.root
// /u/land/s444_2020/r3brootdir/r3broot/macros/r3b/unpack/s467/califa/fastChecking.C
// /u/land/s444_2020/r3brootdir/leylar3b/macros/r3b/unpack/s444/Sofia_Califa_p2p_LA.C
// /d/land4/202006_s444/califa/rootfiles/califa_main0023.root
// root files are generated using Generate_rootfile_califa.C 
// lmd files of the proton beam time /d/land4/202006_s444 stiched data :
//  run28 Ebeam=250 MeV --> 75 < opangle < 85
//  run24,27 Ebeam=300 MeV --> 76 < opangle < 86 
//  run23 Ebeam=400 MeV --> 75 < opangle < 85
//  Ebeam=500 MeV --> 74 < opangle < 84
//  Ebeam=600 MeV --> 75 < opangle < 85
//  Ebeam=800 MeV --> 75 < opangle < 85
//  Ebeam=1 GeV --> 75 < opangle < 85
//
//


// define the Classes

TChain* chain;
TClonesArray* hitCA; 
TBranch* bhitCA;
TClonesArray* mapCA; 
TBranch* bmapCA;
TClonesArray* calCA; 
TBranch* bcalCA;
TClonesArray* evntHeader; // event Header
TBranch* bevntHeader;

// CALIFA parameters
const double THETA_SPREAD = 1.5; // in degrees, estimate for now
const double PHI_SPREAD = 2.8125; //2.8125; // in degrees, 360 deg/2 halves/16 preAmp/4 crystals  for Barrel
const double BETA = 0.75; // for 46-48K(p,2p)45-47Ar for 488 MeV/u before target, beta before target: 0.755 and after target 0.728 --> avarage 0.7415 for  run 291

const double EPSUMMIN = 0;
const double EPSUMMAX = 1e7;
const int NPCALIFA = 2; //0: all, 1: 1p, 2: 2p (p,2p)
const double EPMIN = 20e3; // 50 MeV for proton threshold
const int NMAXCRYSPERHIT = 20;
double opangmin = 75, opangmax=85;
//const double Ebeam = 600; 

char fname[200];
TH2D* hewixvsmes0;
TH2D* hewixvsmes1;
TH2D* hewixvsmes2;
TH2D* hmaxEwixvsmes;
TH2D* hmaxEwixvsmescut;
TH2D* hmaxEvsTwix;
TH2D* hmaxEvsTmes;
TH2D* hmaxEvsTwixopangcut;
TH2D* hmaxEvsTmesopangcut;
TH2D* hmaxEvsTmesTwixcut;
TH2D* hmaxEvsTwixTmescut;
TH2D* hmaxEvsTwixTmescut1;
TH2D* hmaxEvsTmesTwixcut1;
TH1D* hmaxEwixTmescut1;
TH1D* hmaxEmesTwixcut1;
TH2D* hmaxEvsTwixTmescut2;
TH2D* hmaxEvsTmesTwixcut2;
TH1D* hmaxEwixTmescut2;
TH1D* hmaxEmesTwixcut2;
TH1D* hmaxTwixTmescut1;
TH1D* hmaxTmesTwixcut1;
TH2D* hmaxEvsClusterwix;
TH2D* hmaxEvsClustermes;
TH2D* hwixevst;
TH2D* hmesevst0;
TH2D* hmesevsopang;
TH2D* hwixevsopang;
TH2D* hmesevst1;
TH1D* hwixthetalab2;
TH1D* hmesthetalab2;
TH1D* hwixthetacm2;
TH1D* hmesthetacm2;
TH1D* hwixelab2;
TH1D* hmeselab2;
TH2D* hmesevst2;
TH2D* hmesevstcut1;
TH2D* hmesevstcut2;
TH2D* hesumvsopangle;
TH1D* hepsum; 
TH1D* hhitmult;
TH1D* hmapmult;
TH1D* hcalmult;
TH2D* hthetawixvsmes;
TH1D* hopanglewixmes;
TH1D* hopeningangle;
TH2D* hewixvsmescut;
TH1D* hsumthetawixmes;
TH1D* hphiwixmes;
TH1D* henergy_undop;
TH1D* henergy_dop;
TH2D* hphiwixvsmes;
TH2D* hthetaegamma;
TH2D* hmultegamma;
TH2D* hopeningegamma;
TH2D* hmapEvscrystalID;
TH2D* hcalEvscrystalID;
  
TCanvas* can = new TCanvas("can", "can", 1600, 1200);

double GetLabToRestFrameDopplerCorrectionFactor(double beta, double theta){
  return (1.-beta*cos(theta))/sqrt(1.-beta*beta);
}

void AnalyzeRun(int runNum){

  chain->Reset();
  //sprintf(fname, "run%d_sorted.root", runNum);
  sprintf(fname, "/d/land4/202006_s444/califa/rootfiles/califa_main%04d.root", runNum);
  chain->Add(fname);

  if(runNum==23 || runNum==38){ opangmin = 75.; opangmax = 85.;}// 400 MeV
  if(runNum>23 && runNum<28){ opangmin = 76.; opangmax = 86.;} // 300 MeV
  if(runNum==28){ opangmin = 76.; opangmax = 86.;} // 250 MeV
  if(runNum>28 && runNum<32){ opangmin = 74.; opangmax = 84.;} // 500 MeV
  if(runNum>31 && runNum<35){ opangmin = 73.; opangmax = 83.;} // 600 MeV
  if(runNum==36){ opangmin = 72.; opangmax = 82.;}// 800 MeV
  if(runNum==35){ opangmin = 70.; opangmax = 81.;}// 1 GeV

  cout<<" run no: "<<runNum<<" opangmin: "<<opangmin<<" opangmax: "<<opangmax<<endl;

  hitCA = new TClonesArray("R3BCalifaHitData", 5);
  bhitCA = chain->GetBranch("CalifaHitData");
  if (bhitCA!=NULL)
    bhitCA->SetAddress(&hitCA);
 
  mapCA = new TClonesArray("R3BCalifaMappedData", 5);
  bmapCA = chain->GetBranch("CalifaMappedData");
  if (bmapCA!=NULL)
    bmapCA->SetAddress(&mapCA);

  calCA = new TClonesArray("R3BCalifaCrystalCalData", 5);
  bcalCA = chain->GetBranch("CalifaCrystalCalData");
  if (bcalCA!=NULL)
    bcalCA->SetAddress(&calCA);

  evntHeader = new TClonesArray("R3BEventHeader", 5);
  //char evntHeader[100];
  bevntHeader = chain->GetBranch("R3BEventHeader");
  if (bevntHeader!=NULL)
	  bevntHeader->SetAddress(&evntHeader);


  Long64_t nevents = chain->GetEntries();
  int nHitsCA = 0; // number of hits per event
  int nMapCA = 0; // number of map hits per event
  int nCalCA = 0; // number of cal hits per event
  int nevntHeader = 0;

  vector<double> ve_lab; // hit energy, in lab frame
  vector<double> ve_cor; // hit energy, Doppler corrected
  vector<double> vtheta; // polar angle of a hit
  vector<double> vmult; // crystal multiplicity in a hit
  
  for(Long64_t i = 0; i <50; i++){
    hitCA->Clear();
    mapCA->Clear();
    calCA->Clear();
    evntHeader->Clear();

    if(i%100000==0)
      cout<<"Event "<<i<<"/"<<nevents<<endl;
    chain->GetEvent(i);
    nHitsCA = hitCA->GetEntries();
    hhitmult->Fill(nHitsCA);

    nMapCA = mapCA->GetEntries();
    hmapmult->Fill(nMapCA);  
    
    nCalCA = calCA->GetEntries();
    hcalmult->Fill(nCalCA); 

    double califa_e[nHitsCA];
    double califa_theta[nHitsCA];
    double califa_phi[nHitsCA];

    double califa_emap[nMapCA];
    double califa_mapcrystalID[nMapCA];

    double califa_ecal[nCalCA];
    double califa_calcrystalID[nCalCA];

   //if((Tpat&2)!=2) continue;

   if (nHitsCA > 0 ) 
    {
      
      ve_lab.clear();
      ve_cor.clear();
      vtheta.clear();
      vmult.clear();
      double emes = 0.;
      double ewix = 0.;
      double thetames, thetawix, phimes, phiwix, theta, phi;
      TVector3 vecmes, vecwix;


       R3BEventHeader *evntdat = (R3BEventHeader*)evntHeader->At(0);
       cout<<" Tpat: ";
       //cout << ' ' << (int)evntHeader[i];
       //for (int i = 0; i < 100; ++i) cout << ' ' << (int)evntHeader[i];
       cout<<evntdat.GetTpat();
       cout << '\n';


      for(int j = 0; j < nHitsCA; j++){

  	R3BCalifaHitData *cdat = (R3BCalifaHitData*)hitCA->At(j);

	califa_e[j] = cdat->GetEnergy();
	califa_theta[j] = cdat->GetTheta()*TMath::RadToDeg();// + gRandom->Uniform(-THETA_SPREAD, THETA_SPREAD);;
	califa_phi[j] = cdat->GetPhi()*TMath::RadToDeg();

	if (cdat->GetEnergy() > EPMIN){ // proton hit candidate

		phi = cdat->GetPhi()*TMath::RadToDeg();

		if (abs(phi) > 90. && cdat->GetEnergy() > ewix){ // wix
		  ewix = cdat->GetEnergy();
                  //thetawix = cdat->GetTheta()*TMath::RadToDeg()+gRandom->Uniform(-THETA_SPREAD, THETA_SPREAD);
		  thetawix = cdat->GetTheta()*TMath::RadToDeg();
	          //cout<<" theta 1: "<<cdat->GetTheta()*TMath::RadToDeg()<<" random: "<<gRandom->Uniform(-THETA_SPREAD, THETA_SPREAD)<<" theta 2: "<<thetawix<<endl;
	          //phiwix = cdat->GetPhi()*TMath::RadToDeg()+gRandom->Uniform(-PHI_SPREAD, PHI_SPREAD);
	          phiwix = phi;
	          vecwix.SetMagThetaPhi(1., thetawix*TMath::DegToRad(), phiwix*TMath::DegToRad());
	  }

	  else if (abs(phi) < 90. && cdat->GetEnergy() > emes){ // mes
	    emes = cdat->GetEnergy();
	    //thetames = cdat->GetTheta()*TMath::RadToDeg()+gRandom->Uniform(-THETA_SPREAD, THETA_SPREAD);
	    thetames = cdat->GetTheta()*TMath::RadToDeg();
	    //phimes = cdat->GetPhi()*TMath::RadToDeg()+gRandom->Uniform(-PHI_SPREAD, PHI_SPREAD);;
	    phimes = phi;
	    vecmes.SetMagThetaPhi(1., thetames*TMath::DegToRad(), phimes*TMath::DegToRad());
	  }  
	}

	// ************************  gammas ******************************
	else if (cdat->GetNbOfCrystalHits()<=NMAXCRYSPERHIT){ // save gamma-ray candidate information
	  double elab = cdat->GetEnergy();
	  double theta  = cdat->GetTheta()*TMath::RadToDeg()+gRandom->Uniform(-THETA_SPREAD, THETA_SPREAD);
	  double doppler_factor = GetLabToRestFrameDopplerCorrectionFactor(BETA, theta);
	  ve_lab.push_back(elab);
	  ve_cor.push_back(elab*doppler_factor);
	  vtheta.push_back(theta);
	  vmult.push_back(cdat->GetNbOfCrystalHits());
	}

	// ***************************************************************
      }

      if (NPCALIFA==0){ // take all gamma rays, based on the fragment PID gate
	for(int j = 0; j < ve_lab.size(); j++){
	  henergy_undop->Fill(ve_lab[j]);	    
	  henergy_dop->Fill(ve_cor[j]);
	  hthetaegamma->Fill(ve_cor[j], vtheta[j]*TMath::RadToDeg());
	  hmultegamma->Fill(ve_cor[j], vmult[j]);
	}
      }

      else if (NPCALIFA==1){ // ask for at least 1 high-energy hit in either side of CALIFA
	if (ewix > EPMIN || emes > EPMIN){
	  for(int j = 0; j < ve_lab.size(); j++){
	    henergy_undop->Fill(ve_lab[j]);	    
	    henergy_dop->Fill(ve_cor[j]);
	    hthetaegamma->Fill(ve_cor[j], vtheta[j]*TMath::RadToDeg());
	    hmultegamma->Fill(ve_cor[j], vmult[j]);
	  }
	}	
      }

      //******************************** finding max energy ***********************************
    
        double maxEWix=-1., maxEMes=-1.; 
	TVector3 master[2];

        for(int j = 0; j < nHitsCA; j++){

        if(califa_e[j]>maxEWix && TMath::Abs(califa_phi[j])>=90.){

		maxEWix = califa_e[j];
		master[0].SetMagThetaPhi(1.,califa_theta[j]*TMath::DegToRad(),califa_phi[j]*TMath::DegToRad());
	}

	if(califa_e[j]>maxEMes && TMath::Abs(califa_phi[j])<90.){

		maxEMes = califa_e[j];
		master[1].SetMagThetaPhi(1.,califa_theta[j]*TMath::DegToRad(),califa_phi[j]*TMath::DegToRad());
	}
	}
       //cout<<"event no: "<<i<<" hit: "<<nHitsCA<<" maxEwix: "<<maxEWix<<" theta wix: "<<master[0].Theta()*TMath::RadToDeg()<<" phi wix: "<<master[0].Phi()*TMath::RadToDeg()<<" maxEMes: "<<maxEMes<<" theta mes: "<<master[1].Theta()*TMath::RadToDeg()<<" phi mes: "<<master[1].Phi()*TMath::RadToDeg()<<"\n"<<endl; 	


        if(maxEWix>EPMIN && maxEMes>EPMIN){

		double openingangle = master[0].Angle(master[1])*TMath::RadToDeg();
                hopeningangle->Fill(openingangle);
		hmaxEwixvsmes->Fill(maxEWix,maxEMes);
		hmaxEvsTwix->Fill(master[0].Theta()*TMath::RadToDeg(),maxEWix);
		hmaxEvsTmes->Fill(master[1].Theta()*TMath::RadToDeg(),maxEMes);
		if(openingangle > opangmin && openingangle < opangmax) hmaxEwixvsmescut->Fill(maxEWix,maxEMes);
                if(master[0].Theta()*TMath::RadToDeg()>60.0) hmaxEvsTmesTwixcut->Fill(master[1].Theta()*TMath::RadToDeg(),maxEMes);
		if(master[1].Theta()*TMath::RadToDeg()>60.0) hmaxEvsTwixTmescut->Fill(master[0].Theta()*TMath::RadToDeg(),maxEWix);
                if(master[0].Theta()*TMath::RadToDeg()>60.0 && (maxEWix>75000. && maxEWix<125000.)){
			hmaxEvsTmesTwixcut1->Fill(master[1].Theta()*TMath::RadToDeg(),maxEMes);
			hmaxEmesTwixcut1->Fill(maxEMes);
			hmaxTmesTwixcut1->Fill(master[1].Theta()*TMath::RadToDeg());
		}
                if(master[1].Theta()*TMath::RadToDeg()>60.0 && (maxEMes>75000. && maxEMes<125000.)){
			hmaxEvsTwixTmescut1->Fill(master[0].Theta()*TMath::RadToDeg(),maxEWix);
			hmaxEwixTmescut1->Fill(maxEWix);
			hmaxTwixTmescut1->Fill(master[0].Theta()*TMath::RadToDeg());
		}
                if(master[0].Theta()*TMath::RadToDeg()>60.0 && maxEWix>125000.) {hmaxEvsTmesTwixcut2->Fill(master[1].Theta()*TMath::RadToDeg(),maxEMes);hmaxEmesTwixcut2->Fill(maxEMes);}
                if(master[1].Theta()*TMath::RadToDeg()>60.0 && maxEMes>125000.) {hmaxEvsTwixTmescut2->Fill(master[0].Theta()*TMath::RadToDeg(),maxEWix);hmaxEwixTmescut2->Fill(maxEWix);}
		hmaxEvsClusterwix->Fill(nHitsCA,maxEWix);
		hmaxEvsClustermes->Fill(nHitsCA,maxEMes);
		if(openingangle>100.){
			hmaxEvsTwixopangcut->Fill(master[0].Theta()*TMath::RadToDeg(),maxEWix);
		        hmaxEvsTmesopangcut->Fill(master[1].Theta()*TMath::RadToDeg(),maxEMes);
		}

	}

      //*********************************** (p,2p) events *******************************************************

      if (NPCALIFA==2 && nHitsCA==2){ // require both sides of CALIFA to be hit with protons
	if (ewix > EPMIN && emes > EPMIN){ // && (thetawix+thetames)>70. && (thetawix+thetames)<90.){
	  hewixvsmes2->Fill(ewix, emes);
	  hwixevst->Fill(thetawix, ewix);
	  hmesevst2->Fill(thetames, emes);
	  hwixthetalab2->Fill(thetawix);
	  hmesthetalab2->Fill(thetames);
	  hwixelab2->Fill(ewix);
	  hmeselab2->Fill(emes);
	  hepsum->Fill(ewix+emes);
	  hthetawixvsmes->Fill(thetawix, thetames);
	  hphiwixvsmes->Fill(phiwix, phimes);	
	  hsumthetawixmes->Fill(thetawix+thetames);
	  hphiwixmes->Fill(abs(phiwix-phimes));
	  hopanglewixmes->Fill(vecwix.Angle(vecmes)*TMath::RadToDeg());
	  double opangle = vecwix.Angle(vecmes)*TMath::RadToDeg();
	  hmesevsopang->Fill(opangle, emes);
	  hwixevsopang->Fill(opangle, ewix);
	  if(opangle > opangmin && opangle < opangmax) hewixvsmescut->Fill(ewix, emes);
	  if(opangle > opangmin && opangle < opangmax) hmesevstcut2->Fill(thetames,emes);
	  hesumvsopangle->Fill(vecwix.Angle(vecmes)*TMath::RadToDeg(), ewix+emes);
	  if (ewix+emes > EPSUMMIN && ewix+emes < EPSUMMAX){
	    for(int j = 0; j < ve_lab.size(); j++){
	      henergy_undop->Fill(ve_lab[j]);	    
	      henergy_dop->Fill(ve_cor[j]);
	      hthetaegamma->Fill(ve_cor[j], vtheta[j]*TMath::RadToDeg());
	      hmultegamma->Fill(ve_cor[j], vmult[j]);
	      hopeningegamma->Fill(ve_cor[j], thetawix+thetames);
	    }
	  }
	}
      } //**************************** (p,2p) event loop ******************************************************

	//****************************** map level events *****************************************************
	
	if(nMapCA>0){
	
	for(int imap=0; imap<nMapCA; imap ++){
        
		R3BCalifaMappedData *cmapdat = (R3BCalifaMappedData*)mapCA->At(imap);
		hmapEvscrystalID->Fill(cmapdat->GetCrystalId(),cmapdat->GetEnergy());


	}
	}
	
	//****************************** map level events *****************************************************

	//****************************** cal level events *****************************************************
	
	if(nCalCA>0){
	
	for(int ical=0; ical<nCalCA; ical ++){
        
		R3BCalifaCrystalCalData *ccaldat = (R3BCalifaCrystalCalData*)calCA->At(ical);
		hcalEvscrystalID->Fill(ccaldat->GetCrystalId(),ccaldat->GetEnergy());


	}
	}
	//****************************** cal level events *****************************************************

     if(nHitsCA>0){
	 hmesevst0->Fill(thetames,emes); 
	 hewixvsmes0->Fill(ewix, emes);
     }

     if(nHitsCA==1){

       hmesevstcut1->Fill(thetames,emes); 
       hmesevst1->Fill(thetames,emes); 
       hewixvsmes1->Fill(emes, ewix+emes);

     }
    }
  }
  sprintf(fname, "run%d_histos.root", runNum);
  TFile* outfile = new TFile(fname, "RECREATE");
  hewixvsmes0->Write();
  hewixvsmes1->Write();
  hewixvsmes2->Write();
  hmaxEwixvsmes->Write();
  hmaxEwixvsmescut->Write();
  hewixvsmescut->Write();
  hmaxEvsTwix->Write();
  hmaxEvsTmes->Write();
  hmaxEvsTwixopangcut->Write();
  hmaxEvsTmesopangcut->Write();
  hmaxEvsTmesTwixcut->Write();
  hmaxEvsTwixTmescut->Write();
  hmaxEvsTmesTwixcut1->Write();
  hmaxEvsTwixTmescut1->Write();
  hmaxEmesTwixcut1->Write();
  hmaxEwixTmescut1->Write();
  hmaxEvsTmesTwixcut2->Write();
  hmaxEvsTwixTmescut2->Write();
  hmaxEmesTwixcut2->Write();
  hmaxEwixTmescut2->Write();
  hmaxEvsClusterwix->Write();
  hmaxTmesTwixcut1->Write();
  hmaxTwixTmescut1->Write();
  hmaxEvsClustermes->Write();
  hesumvsopangle->Write();
  hwixevsopang->Write();
  hmesevsopang->Write();
  hwixthetalab2->Write();
  hmesthetalab2->Write();
  hwixelab2->Write();
  hmeselab2->Write();
  hepsum->Write();
  hhitmult->Write();
  hmapmult->Write();
  hcalmult->Write();
  hwixevst->Write();
  hmesevst0->Write();
  hmesevst1->Write();
  hmesevst2->Write();
  hmesevstcut1->Write();
  hmesevstcut2->Write();
  hthetawixvsmes->Write();
  hphiwixvsmes->Write();
  hopanglewixmes->Write();
  hopeningangle->Write();
  hsumthetawixmes->Write();
  hphiwixmes->Write();
  hphiwixvsmes->Write();
  hmapEvscrystalID->Write();
  hcalEvscrystalID->Write();

  //henergy_undop->Write();
  //henergy_dop->Write();
  //hthetaegamma->Write();
  //hmultegamma->Write();
  //hopeningegamma->Write();
  outfile->Close();

} // end of event loop

void DrawResults(){    
  can->Divide(2, 2);
  can->cd(1);
  hewixvsmes2->Draw("colz");
  can->cd(2);
  hthetawixvsmes->Draw("colz");
  can->cd(3);
  hopanglewixmes->Draw();
  hsumthetawixmes->Draw("same");
  hphiwixmes->Draw("same");
  hopanglewixmes->GetYaxis()->SetRangeUser(0, hopanglewixmes->GetMaximum()*2.0);
  TLegend* lh = new TLegend(0.65, 0.6, 0.9, 0.9);
  lh->AddEntry(hopanglewixmes, "Opening angle", "l");
  lh->AddEntry(hsumthetawixmes, "#theta_{wix} + #theta_{mes}", "l");
  lh->AddEntry(hphiwixmes, "|#phi_{wix} #minus #phi_{mes}|", "l");
  lh->Draw();
  can->cd(4);
  hphiwixvsmes->Draw("colz");

}
void DefineHistograms(){
 
  hewixvsmes0 = new TH2D("ewixvsmes0", "E_{wix} vs E_{mes} for pmul>0;E_{wix} (keV); E_{mes} (keV)", 200, 0, 5e5, 200, 0, 5e5);
  hewixvsmes1 = new TH2D("ewixvsmes1", "E_{wix}+E_{mes} vs E_{mes} for pmul=1;E_{mes} (keV); E_{wix+mes} (keV)", 200, 0, 5e5, 200, 0, 5e5);
  hewixvsmes2 = new TH2D("ewixvsmes2", "E_{wix} vs E_{mes} for pmul=2;E_{wix} (keV); E_{mes} (keV)", 200, 0, 5e5, 200, 0, 5e5);
  hmaxEwixvsmes = new TH2D("maxEwixvsmes", "max E_{wix} vs E_{mes};E_{wix} (keV); E_{mes} (keV)", 200, 0, 5e5, 200, 0, 5e5);
  hmaxEvsTwix = new TH2D("maxEvsTwix", "max E_{wix} vs #theta_{wix};#theta_{wix} (deg); E_{wix} (deg)", 90, 0, 90, 200, 0, 5e5);
  hmaxEvsTmes = new TH2D("maxEvsTmes", "max E_{mes} vs #theta_{mes};#theta_{mes} (deg); E_{mes} (keV)", 90, 0, 90, 200, 0, 5e5);
  hmaxEvsTwixopangcut = new TH2D("maxEvsTwixopangcut", "max E_{wix} vs #theta_{wix} for opening angle > 100;#theta_{wix} (deg); E_{wix} (deg)", 90, 0, 90, 200, 0, 5e5);
  hmaxEvsTmesopangcut = new TH2D("maxEvsTmesopangcut", "max E_{mes} vs #theta_{mes}for opening angle > 100;#theta_{mes} (deg); E_{mes} (keV)", 90, 0, 90, 200, 0, 5e5);
  hmaxEvsTmesTwixcut = new TH2D("maxEvsTmesTwixcut", "max E_{mes} vs #theta_{mes} for #theta_{wix}>60 ;#theta_{mes} (deg); E_{mes} (keV)", 90, 0, 90, 200, 0, 5e5);
  hmaxEvsTwixTmescut = new TH2D("maxEvsTwixTmescut", "max E_{wix} vs #theta_{wix} for #theta_{mes}>60 ;#theta_{wix} (deg); E_{wix} (keV)", 90, 0, 90, 200, 0, 5e5);
  hmaxEvsTmesTwixcut1 = new TH2D("maxEvsTmesTwixcut1", "max E_{mes} vs #theta_{mes} for #theta_{wix}>60 & 75<E_{wix}<125;#theta_{mes} (deg); E_{mes} (keV)", 90, 0, 90, 200, 0, 5e5);
  hmaxEvsTwixTmescut1 = new TH2D("maxEvsTwixTmescut1", "max E_{wix} vs #theta_{wix} for #theta_{mes}>60 & 75<E_{mes}<125;#theta_{wix} (deg); E_{wix} (keV)", 90, 0, 90, 200, 0, 5e5);
  hmaxEmesTwixcut1 = new TH1D("maxEmesTwixcut1", "max E_{mes} for #theta_{wix}>60 & 75<E_{wix}<125;E_{mes} (keV)", 200, 0, 5e5);
  hmaxEwixTmescut1 = new TH1D("maxEwixTmescut1", "max E_{wix} for #theta_{mes}>60 & 75<E_{mes}<125;E_{wix} (keV)", 200, 0, 5e5);
  hmaxEvsTmesTwixcut2 = new TH2D("maxEvsTmesTwixcut2", "max E_{mes} vs #theta_{mes} for #theta_{wix}>60 & E_{wix}>125;#theta_{mes} (deg); E_{mes} (keV)", 90, 0, 90, 200, 0, 5e5);
  hmaxEvsTwixTmescut2 = new TH2D("maxEvsTwixTmescut2", "max E_{wix} vs #theta_{wix} for #theta_{mes}>60 & E_{mes}>125;#theta_{wix} (deg); E_{wix} (keV)", 90, 0, 90, 200, 0, 5e5);
  hmaxEmesTwixcut2 = new TH1D("maxEmesTwixcut2", "max E_{mes} for #theta_{wix}>60 & E_{wix}>125;E_{mes} (keV)", 200, 0, 4e5);
  hmaxEwixTmescut2 = new TH1D("maxEwixTmescut2", "max E_{wix} for #theta_{mes}>60 & E_{mes}>125;E_{wix} (keV)", 200, 0, 4e5);
  hmaxTmesTwixcut1 = new TH1D("maxTmesTwixcut1", "max E_{mes} for #theta_{wix}>60 & 75<E_{wix}<125;#theta_{mes} (deg)", 100, 0, 100);
  hmaxTwixTmescut1 = new TH1D("maxTwixTmescut1", "max E_{wix} for #theta_{mes}>60 & 75<E_{mes}<125;#theta_{wix} (deg)", 100, 0, 100);
  hmaxEvsClusterwix = new TH2D("maxEvsClusterwix", "max E_{wix} vs Clusters multiplicity;Cluster;E_{wix} (keV)", 50, 0, 50,200, 0, 4e5);
  hmaxEvsClustermes = new TH2D("maxEvsClustermes", "max E_{mes} vs Clusters multiplicity;Cluster;E_{mes} (keV)", 50, 0, 50,200, 0, 4e5);
  hmaxEwixvsmescut = new TH2D("maxEwixvsmescut", "max E_{wix} vs E_{mes} for cut 75 < opening angle < 85;E_{wix} (keV); E_{mes} (keV)", 200, 0, 5e5, 200, 0, 5e5);
  hwixevsopang = new TH2D("wixevsopang", "Wix (Right) E vs opening angle for pmul=2; #theta_{opang} (deg);E_{mes} (keV)", 200, 0, 200, 1000, 0, 5e5);
  hmesevsopang = new TH2D("mesevsopang", "Messel (Left) E vs opening angle for pmul=2; #theta_{opang} (deg);E_{mes} (keV)", 200, 0, 200, 1000, 0, 5e5);
  hewixvsmescut = new TH2D("ewixvsmescut", "E_{wix} vs E_{mes} for cut 75 < opening angle < 85 ;E_{wix} (keV); E_{mes} (keV)", 200, 0, 5e5, 200, 0, 5e5);
  hwixevst = new TH2D("wixevst", "Wix (Right) E vs #theta; #theta_{wix} (deg);E_{wix} (keV)", 90, 0, 90, 1000, 0, 5e5);
  hmesevst0 = new TH2D("mesevst0", "Messel (Left) E vs #theta for pmul>0; #theta_{mes} (deg);E_{mes} (keV)", 90, 0, 90, 1000, 0, 5e5);
  hmesevst1 = new TH2D("mesevst1", "Messel (Left) E vs #theta for pmul=1; #theta_{mes} (deg);E_{mes} (keV)", 90, 0, 90, 1000, 0, 5e5);
  hwixthetalab2 = new TH1D("wixthetalab2", "Wix (Right) #theta for pmul=2; #theta_{mes} (deg)", 90, 0, 90);
  hmesthetalab2 = new TH1D("mesthetalab2", "Messel (Left) #theta for pmul=2; #theta_{mes} (deg)", 90, 0, 90);
  hwixelab2 = new TH1D("wixelab2", "Wix (Right) E (keV) for pmul=2; E_{wix} (keV)", 1000, 0, 5e5);
  hmeselab2 = new TH1D("meselab2", "Messel (Left) E (keV) for pmul=2; E_{mes} (keV)", 1000, 0, 5e5);
  hmesevst2 = new TH2D("mesevst2", "Messel (Left) E vs #theta for pmul=2; #theta_{mes} (deg);E_{mes} (keV)", 90, 0, 90, 1000, 0, 5e5);
  hmesevstcut1 = new TH2D("mesevstcut1", "Messel (Left) E vs #theta for pmul=1; #theta_{mes} (deg);E_{mes} (keV)", 90, 0, 90, 1000, 0, 5e5);
  hmesevstcut2 = new TH2D("mesevstcut2", "Messel (Left) E vs #theta for pmul=2 & cut 75 < opening angle < 85; #theta_{mes} (deg);E_{mes} (keV)", 90, 0, 90, 1000, 0, 5e5);
  hepsum = new TH1D("hepsum", "sum energy ;E_{wix} + E_{mes} (keV);", 1000, 0, 1e6);
  hthetawixvsmes = new TH2D("thetawixvsmes", "#theta_{wix} vs #theta_{mes};#theta_{wix} (deg); #theta_{mes} (deg)", 90, 0, 90, 90, 0, 90);
  hopanglewixmes = new TH1D("opanglewixmes", "opening angle (deg)", 200, 0, 200);
  hopeningangle = new TH1D("openingangle", "opening angle (deg)", 200, 0, 200);
  hsumthetawixmes = new TH1D("sumthetawixmes", "sum theta;#theta_{wix} + #theta_{mes} (deg)", 180, 0, 180);
  hesumvsopangle = new TH2D("esumvsopangle", "sum energy vs opening angle; #theta_{opangle}; E_{sum} (keV)", 180, 0, 180, 1000,0,1e6);
  hphiwixmes = new TH1D("hphiwixmes", ";|#phi_{wix} - #phi_{mes}| (deg)", 270, 0, 270);
  henergy_undop = new TH1D("henergy_undop", "henergy_undop", 200, 0, 10000);
  henergy_dop = new TH1D("henergy_dop", "henergy_dop", 200, 0, 10000);
  hopanglewixmes->SetLineColor(1);
  hsumthetawixmes->SetLineColor(2);
  hphiwixmes->SetLineColor(4);
  hhitmult = new TH1D("hitmult", "cluster multiplicity", 100, 0, 100);
  hmapmult = new TH1D("mapmult", "hit  multiplicity in map level;", 500, 0, 500);
  hcalmult = new TH1D("calmult", "hit  multiplicity in cal level;", 500, 0, 500);
  hphiwixvsmes = new TH2D("phiwixvsmes", ";#phi_{wix} (deg) vs #phi_{mes} (deg)", 100, -200, 200, 100, -200, 200);
  
  hthetaegamma = new TH2D("thetaeggamma",";E_{#gamma} (Dopp. corrected) (keV);#theta (deg)", 400, 0, 20000, 90, 0, 90);
  hmultegamma = new TH2D("multegamma",";E_{#gamma} (Dopp. corrected) (keV);Crystal hits in cluster", 400, 0, 20000, 10, 0, 10);
  hopeningegamma = new TH2D("openingeggamma",";E_{#gamma} (Dopp. corrected) (keV);#theta_mes + #theta_wix protons (deg)", 200, 0, 10000, 180, 0, 180);
  
  hmapEvscrystalID = new TH2D("mapEvscrystalID", "Energy in map level vs Crystal ID; crystal ID;energy (arb. units))", 3500, 1500, 5000, 1500,0,30000);
  hcalEvscrystalID = new TH2D("calEvscrystalID", "Energy in cal level vs Crystal ID; crystal ID;energy (keV))", 3500, 1500, 5000, 2000,0,500000);

}
void LA_califa_protonbeam(int runNum){
  chain = new TChain("evt");
  DefineHistograms();
  AnalyzeRun(runNum);
  DrawResults();

}
