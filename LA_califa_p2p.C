// root -l -b -q 'LA_califa_protonbeam.C(23)'; root -l run23_histos.root
// /u/land/s444_2020/r3brootdir/r3broot/macros/r3b/unpack/s467/califa/fastChecking.C
// /u/land/s444_2020/r3brootdir/leylar3b/macros/r3b/unpack/s444/Sofia_Califa_p2p_LA.C
// /d/land4/202006_s444/califa/rootfiles/califa_main0023.root
// root files are generated using Generate_rootfile_califa.C 
// lmd files of the proton beam time /d/land4/202006_s444 stiched data :
// s467: root -l -b -q filltree_new_upexps.C'(305)'
// s467: 
#include "preamp.h"

// define the Classes
TChain* chain;
TClonesArray* hitCA; 
TBranch* bhitCA;
TClonesArray* mapCA; 
TBranch* bmapCA;
TClonesArray* calCA; 
TBranch* bcalCA;
TClonesArray* wrMasterCA; 
TBranch* bwrMasterCA;

// Califa parameters
const double EPMIN = 20e3; //MeV

// Histograms 
TH1D* htpat;
TH1D* htpatspillon;
TH1D* htpatspilloff;
TH1D* htpatshift;
TH1D* htpatshiftspillon;
TH1D* htpatshiftspilloff;
TH1D* htpatshiftcalmap;
TH1D* htpatshiftcalmap0;
TH1D* htpatshiftcalmapspillon;
TH1D* htpatshiftcalmap0spillon;
TH1D* htpatshiftcalmapspilloff;
TH1D* htpatshiftcalmap0spilloff;
TH1D* htrig;
TH1D* hwrMasterTS;
TH1D* hwrMasterTSspilloff;
TH1D* hmapmult;
TH1D* hcalmult;
TH1D* hhitmult;
TH2D* hmapEvscrystalID;
TH2D* hmapEvscrystalIDspillon;
TH2D* hmapEvscrystalIDspilloff;
TH2D* hmapEvscrystalIDspillofftpat0;
TH2D* hmapNfvscrystalID;
TH2D* hmapNsvscrystalID;
TH2D* hmapNfvscrystalIDspilloff;
TH2D* hmapNsvscrystalIDspilloff;
TH2D* hcalEvscrystalID;
TH2D* hcalEvscrystalIDspillofftpat0;
TH1D* hcalEspillofftpat0;
TH2D* hcalNfvscrystalID;
TH2D* hcalNsvscrystalID;
TH2D* hcalNfvscrystalIDspilloff;
TH2D* hcalNsvscrystalIDspilloff;
TH2D* hmapNfNs;
TH2D* hmapNfNsspillon;
TH2D* hmapNfNsspilloff;
TH2D* hmapNfNsGamBar;
TH2D* hmapNfNsGamIph;
TH2D* hmapNfNsProBar;
TH2D* hmapNfNsProIph;
TH2D* hcalNfNs;
TH2D* hcalNfNsgoodtpat;
TH2D* hcalNfNstpat0;
TH2D* hcalNfNsGamBar;
TH2D* hcalNfNsGamIph;
TH2D* hcalNfNsProBar;
TH2D* hcalNfNsProIph;
TH1D* hmapcrystalID;

TCanvas* can = new TCanvas("can", "can", 1600, 1200);

char fname[200];

void AnalyzeRun(int runNum){
  bool onspill{}, tpat0{};
  int64_t duration[2]{}; // duration[0]: offspill, duration[1]: onspill

  int64_t wrts_last{}, wrts_cur{};
  chain->Reset();
  //sprintf(fname, "./rootfiles/rootfiletmp/fragment_Sep2021/s467_filltree_Setting13_%04d_29Oct.root", runNum);
  //sprintf(fname, "./rootfiles/rootfiletmp/fragment_Sep2021/s467_filltree_Setting13_%04d_19Nov.root", runNum);
  sprintf(fname, "./rootfiles/rootfiletmp/fragment_Sep2021/s467_filltree_Setting13_%04d_1Dec.root", runNum);
  chain->Add(fname);

  cout<<" run no: "<<runNum<<endl;

  wrMasterCA = new TClonesArray("R3BWRData", 5);
  bwrMasterCA = chain->GetBranch("WRMasterData");
  if (bwrMasterCA!=NULL)
    bwrMasterCA->SetAddress(&wrMasterCA);
  
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

  //Header Data
  R3BEventHeader* evntHeader = new R3BEventHeader();
  TBranch* bevntHeader = chain->GetBranch("EventHeader.");
  bevntHeader->SetAddress(&evntHeader);
 
  Long64_t nevents = chain->GetEntries();
  int nHitsCA = 0; // number of hits per event
  int nMapCA = 0; // number of map hits per event
  int nCalCA = 0; // number of cal hits per event
  int nWRMasterCA = 0; //number of WR Master hits per event
  int nevntHeader = 0; //number of event header per event

  for(Long64_t i = 0; i <100000; i++){
	  
    hitCA->Clear();
    mapCA->Clear();
    calCA->Clear();
    wrMasterCA->Clear();

    if(i%100000==0)
      cout<<"Event "<<i<<"/"<<nevents<<endl;
    chain->GetEvent(i);
    
    nHitsCA = hitCA->GetEntries();
    hhitmult->Fill(nHitsCA);

    nMapCA = mapCA->GetEntries();
    hmapmult->Fill(nMapCA);  
    
    nCalCA = calCA->GetEntries();
    hcalmult->Fill(nCalCA); 

    nWRMasterCA = wrMasterCA->GetEntries(); 

    //nevntHeader = evntHeader->GetEntries()cout<<" evnt header: "<<nevntHeader<<endl; ;
    
    double califa_e[nHitsCA];
    double califa_theta[nHitsCA];
    double califa_phi[nHitsCA];

    double califa_emap[nMapCA];
    double califa_mapcrystalID[nMapCA];

    double califa_ecal[nCalCA];
    double califa_calcrystalID[nCalCA];

    double mapcrystalID[nMapCA], mapcrystalIDp[nMapCA], mapcrystalIDg[nMapCA];
    double calcrystalID[nCalCA], calcrystalIDp[nCalCA], calcrystalIDg[nCalCA];
     
    //if((Tpat&2)!=2) continue;

    //if (nHitsCA > 0 ) 
    //{

    //if(evntHeader->GetTpat() > 0) cout<<"Tpat: "<<evntHeader->GetTpat()<<" Trigger: "<<evntHeader->GetTrigger()<<endl;
    htpat->Fill(evntHeader->GetTpat());
    htrig->Fill(evntHeader->GetTrigger());
    if(evntHeader->GetTpat()==0) tpat0=true;

    int tpatbin, fTrigger;
    
    if ((evntHeader->GetTpat())&1)// all detectors & Trigger=1 the start detector are fired
      {
        onspill=true;
        htpatspillon->Fill(evntHeader->GetTpat());
      }
    if  ((evntHeader->GetTpat())& 0b1111111111110000)// offspill events;all trigger bits except the first 4 bits (0000) representing onspill events
      {
        onspill=false;
        htpatspilloff->Fill(evntHeader->GetTpat());
      }
        
    if(evntHeader->GetTpat()>0) 
      for(int itpat=0; itpat<16; itpat++)
	{
	  tpatbin = (evntHeader->GetTpat()) & (1<<itpat);
	  if (tpatbin != 0)
            {
              fTrigger = itpat+1;
              //cout<<"tpat: "<<evntHeader->GetTpat()<<" 1<<itpat: "<<(1<<itpat)<<" tpat bin: "<<tpatbin<<" fTrigger: "<<fTrigger<<endl;
              htpatshift->Fill(fTrigger);
              if (nMapCA!=0) htpatshiftcalmap->Fill(fTrigger);  // a bit strange; the main trigger pattern did not see anything
              if (nMapCA==0) htpatshiftcalmap0->Fill(fTrigger);  // a bit strange; the main trigger pattern did not see anything
	      if (onspill==true) htpatshiftspillon->Fill(fTrigger);
	      if (onspill==false) htpatshiftspilloff->Fill(fTrigger);
              if (nMapCA!=0 && onspill==true) htpatshiftcalmapspillon->Fill(fTrigger);  
              if (nMapCA==0 && onspill==true) htpatshiftcalmap0spillon->Fill(fTrigger);  
              if (nMapCA!=0 && onspill==false) htpatshiftcalmapspilloff->Fill(fTrigger);  
              if (nMapCA==0 && onspill==false) htpatshiftcalmap0spilloff->Fill(fTrigger);  
	    }
        }
    else
      {
          htpatshift->Fill(0);  // a bit strange; the main trigger pattern did not see anything
          if (nMapCA!=0) htpatshiftcalmap->Fill(0);  // a bit strange; the main trigger pattern did not see anything
          if (nMapCA==0) htpatshiftcalmap0->Fill(0);  // a bit strange; the main trigger pattern did not see anything
	  if (onspill) htpatshiftspillon->Fill(0);
	  if (!onspill) htpatshiftspilloff->Fill(0);
          if (nMapCA!=0 && onspill==true) htpatshiftcalmapspillon->Fill(0);  
          if (nMapCA==0 && onspill==true) htpatshiftcalmap0spillon->Fill(0);  
          if (nMapCA!=0 && onspill==false) htpatshiftcalmapspilloff->Fill(0);  
          if (nMapCA==0 && onspill==false) htpatshiftcalmap0spilloff->Fill(0);  
      }

    wrts_last=wrts_cur;
    for(int iwrm=0; iwrm<nWRMasterCA; iwrm ++){
      R3BWRData* wrmdat = (R3BWRData*)wrMasterCA->At(iwrm);
      hwrMasterTS->Fill(wrmdat->GetTimeStamp());
      wrts_cur=wrmdat->GetTimeStamp();
      if(onspill==false) hwrMasterTSspilloff->Fill(wrmdat->GetTimeStamp());
    }

    if (wrts_last)
      {
        int64_t diff_to_last=wrts_cur-wrts_last;
        duration[onspill]+=diff_to_last;
      }

    //****************************** map level events *****************************************************
    
    for(int imap=0; imap<nMapCA; imap ++){
        
      R3BCalifaMappedData *cmapdat = (R3BCalifaMappedData*)mapCA->At(imap);
      hmapEvscrystalID->Fill(cmapdat->GetCrystalId(),cmapdat->GetEnergy());
      hmapcrystalID->Fill(cmapdat->GetCrystalId());
      hmapNfNs->Fill(cmapdat->GetNf(),cmapdat->GetNs());
      auto id=cmapdat->GetCrystalId();
      hmapNfvscrystalID->Fill(id,cmapdat->GetNf());
      hmapNsvscrystalID->Fill(id,cmapdat->GetNs());
      if (onspill==true)//spill-on
        {
          hmapEvscrystalIDspillon->Fill(cmapdat->GetCrystalId(),cmapdat->GetEnergy());
          hmapNfNsspillon->Fill(cmapdat->GetNf(),cmapdat->GetNs());
        }
      else//spill-off
        {
        hmapEvscrystalIDspilloff->Fill(cmapdat->GetCrystalId(),cmapdat->GetEnergy());
        if(tpat0==true) hmapEvscrystalIDspillofftpat0->Fill(cmapdat->GetCrystalId(),cmapdat->GetEnergy());
        hmapNfNsspilloff->Fill(cmapdat->GetNf(),cmapdat->GetNs());
        hmapNfvscrystalIDspilloff->Fill(id,cmapdat->GetNf());
        hmapNsvscrystalIDspilloff->Fill(id,cmapdat->GetNs());
	}
      //if (preamp_ranges.at(id)==1) // gamma
        {
          //if(id<1960) hmapNfNsGamBar->Fill(cmapdat->GetNf(),cmapdat->GetNs());
          //if(id>1960 && id<2432) hmapNfNsGamIph->Fill(cmapdat->GetNf(),cmapdat->GetNs());
          if(id<2432) hmapNfNsGamBar->Fill(cmapdat->GetNf(),cmapdat->GetNs());
          if(id>1960 && id<2432) hmapNfNsGamIph->Fill(cmapdat->GetNf(),cmapdat->GetNs());
        }
      //if (preamp_ranges.at(id)==2) //proton
        {
           if(id>2432 && id<4380) hmapNfNsProBar->Fill(cmapdat->GetNf(),cmapdat->GetNs());
           if(id>4380) hmapNfNsProIph->Fill(cmapdat->GetNf(),cmapdat->GetNs());
           
        }
    }

    //for(int imapcr=945; imapcr<1905; imapcr+64){
    //	    if(imapcr<imapcr+32)
    //  }
    //std::cout << preamp_ranges[945];
    //****************************** map level events *****************************************************

    //****************************** cal level events *****************************************************
    for(int ical=0; ical<nCalCA; ical ++){
        
      R3BCalifaCrystalCalData *ccaldat = (R3BCalifaCrystalCalData*)calCA->At(ical);
      hcalEvscrystalID->Fill(ccaldat->GetCrystalId(),ccaldat->GetEnergy());
      hcalNfNs->Fill(ccaldat->GetNf(),ccaldat->GetNs());
      if(evntHeader->GetTpat()>0)
	      hcalNfNsgoodtpat->Fill(ccaldat->GetNf(),ccaldat->GetNs());
	      else
	      hcalNfNstpat0->Fill(ccaldat->GetNf(),ccaldat->GetNs());
      
      auto id=ccaldat->GetCrystalId();
      hcalNfvscrystalID->Fill(id,ccaldat->GetNf());
      hcalNsvscrystalID->Fill(id,ccaldat->GetNs());
      if(onspill==false && tpat0==true)
        {
          hcalEvscrystalIDspillofftpat0->Fill(ccaldat->GetCrystalId(),ccaldat->GetEnergy());
          hcalEspillofftpat0->Fill(ccaldat->GetEnergy());
          hcalNfvscrystalIDspilloff->Fill(id,ccaldat->GetNf());  
          hcalNsvscrystalIDspilloff->Fill(id,ccaldat->GetNs());  
	}
      //if (preamp_ranges.at(id)==1) // gamma
        {
          //if(id<1960) hcalNfNsGamBar->Fill(ccaldat->GetNf(),ccaldat->GetNs());
          //if(id>1960 && id<2432) hcalNfNsGamIph->Fill(ccaldat->GetNf(),ccaldat->GetNs());
          if(id<2432) hcalNfNsGamBar->Fill(ccaldat->GetNf(),ccaldat->GetNs());
          if(id>1960 && id<2432) hcalNfNsGamIph->Fill(ccaldat->GetNf(),ccaldat->GetNs());
        }
      //if (preamp_ranges.at(id)==2) //proton
        {
           if(id>2432 && id<4380) hcalNfNsProBar->Fill(ccaldat->GetNf(),ccaldat->GetNs());
           if(id>4380) hcalNfNsProIph->Fill(ccaldat->GetNf(),ccaldat->GetNs());
           
        }
    }
    //****************************** cal level events *****************************************************
   
    if (nHitsCA > 0 ) 
      {
 
        double emes = 0.;
        double ewix = 0.;
        double thetames, thetawix, phimes, phiwix, theta, phi;
        TVector3 vecmes, vecwix;

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
        //******************************** finding max energy ***********************************

      }
  }
  sprintf(fname, "run%d_histos.root", runNum);
  TFile* outfile = new TFile(fname, "RECREATE");
 
  htpat->Write();
  htpatspillon->Write();
  htpatspilloff->Write();
  htpatshift->Write();
  htpatshiftspillon->Write();
  htpatshiftspilloff->Write(); 
  htpatshiftcalmap->Write();
  htpatshiftcalmap0->Write();
  htpatshiftcalmapspillon->Write();
  htpatshiftcalmap0spillon->Write();
  htpatshiftcalmapspilloff->Write();
  htpatshiftcalmap0spilloff->Write();
  htrig->Write();
  hwrMasterTS->Write();
  hwrMasterTSspilloff->Write();

  hhitmult->Write();
  hmapmult->Write();
  hcalmult->Write();
  
  hmapEvscrystalID->Write();
  hmapEvscrystalIDspillon->Write();
  hmapEvscrystalIDspilloff->Write();
  hmapEvscrystalIDspillofftpat0->Write();
  hmapNfvscrystalID->Write();
  hmapNsvscrystalID->Write();
  hmapNfvscrystalIDspilloff->Write();
  hmapNsvscrystalIDspilloff->Write();
  hmapcrystalID->Write();
  hcalEvscrystalID->Write();
  hcalEvscrystalIDspillofftpat0->Write();
  hcalEspillofftpat0->Write();
  hcalNfvscrystalID->Write();
  hcalNsvscrystalID->Write();
  hcalNfvscrystalIDspilloff->Write();
  hcalNsvscrystalIDspilloff->Write();
  
  hmapNfNs->Write();
  hmapNfNsspillon->Write();
  hmapNfNsspilloff->Write();
  hmapNfNsGamBar->Write();
  hmapNfNsGamIph->Write();
  hmapNfNsProBar->Write();
  hmapNfNsProIph->Write();
  hcalNfNs->Write();
  hcalNfNsgoodtpat->Write();
  hcalNfNstpat0->Write();
  hcalNfNsGamBar->Write();
  hcalNfNsGamIph->Write();
  hcalNfNsProBar->Write();
  hcalNfNsProIph->Write();

  outfile->Close();
  std::cout << "onspill duration: " << duration[1] <<"\noffspill duration: " << duration[0]<<"\n";
} // end of event loop

void DrawResults(){

  gStyle->SetOptStat(11111111);
  
  can->Divide(4, 2);
  can->SetTicks();
  can->cd(1);
  gPad->SetLogy();
  htpatshiftspillon->Draw();
  htpatshiftspillon->SetLineWidth(2);
  can->cd(2);
  gPad->SetLogy();
  htpatshiftspilloff->Draw();
  htpatshiftspilloff->SetLineWidth(2);
  can->cd(3);
  gPad->SetLogy();
  htpatshiftcalmap->Draw();
  htpatshiftcalmap->SetLineWidth(2);
  can->cd(4);
  gPad->SetLogy();
  htpatshiftcalmap0->Draw();
  htpatshiftcalmap0->SetLineWidth(2);
  can->cd(5);
  gPad->SetLogy();
  htpatshiftcalmapspillon->Draw();
  htpatshiftcalmapspillon->SetLineWidth(2);
  can->cd(6);
  gPad->SetLogy();
  htpatshiftcalmap0spillon->Draw();
  htpatshiftcalmap0spillon->SetLineWidth(2);
  can->cd(7);
  gPad->SetLogy();
  htpatshiftcalmapspilloff->Draw();
  htpatshiftcalmapspilloff->SetLineWidth(2);
  can->cd(8);
  gPad->SetLogy();
  htpatshiftcalmap0spilloff->Draw();
  htpatshiftcalmap0spilloff->SetLineWidth(2);

  //hopanglewixmes->Draw();
  //hopanglewixmes->GetYaxis()->SetRangeUser(0, hopanglewixmes->GetMaximum()*2.0);
  //TLegend* lh = new TLegend(0.65, 0.6, 0.9, 0.9);
  //lh->AddEntry(hopanglewixmes, "Opening angle", "l");
  //lh->AddEntry(hsumthetawixmes, "#theta_{wix} + #theta_{mes}", "l");
  //lh->AddEntry(hphiwixmes, "|#phi_{wix} #minus #phi_{mes}|", "l");
  //lh->Draw();

}
void DefineHistograms(){
 
  htpat = new TH1D("Tpat", "Tpat", 1001, -0.5, 1000.5);
  htpatspillon = new TH1D("Tpatspillon", "Tpat for onspill events", 1001, -0.5, 1000.5);
  htpatspilloff = new TH1D("Tpatspilloff", "Tpat for offspill events", 1001, -0.5, 1000.5);
  htpatshift = new TH1D("Tpatshift", "Tpat shifted", 1001, -0.5, 1000.5);
  htpatshiftspillon  = new TH1D("Tpatshiftspillon", "Tpat shifted for onspill events", 13, -0.5, 12.5);
  htpatshiftspilloff  = new TH1D("Tpatshiftspilloff", "Tpat shifted for offspill events", 13, -0.5, 12.5);
  htpatshiftcalmap  = new TH1D("Tpatshiftcalmap", "Tpat shifted for califa map !=0 events", 13, -0.5, 12.5);
  htpatshiftcalmap0  = new TH1D("Tpatshiftcalmap0", "Tpat shifted for califa map ==0 events", 13, -0.5, 12.5);
  htpatshiftcalmapspillon  = new TH1D("Tpatshiftcalmapspillon", "Tpat shifted for califa map !=0 onspill events", 13, -0.5, 12.5);
  htpatshiftcalmap0spillon  = new TH1D("Tpatshiftcalmap0spillon", "Tpat shifted for califa map ==0 onspill events", 13, -0.5, 12.5);
  htpatshiftcalmapspilloff  = new TH1D("Tpatshiftcalmapspilloff", "Tpat shifted for califa map !=0 offspill events", 13, -0.5, 12.5);
  htpatshiftcalmap0spilloff  = new TH1D("Tpatshiftcalmap0spilloff", "Tpat shifted for califa map ==0 offspill events", 13, -0.5, 12.5);
  htrig = new TH1D("Trigger", "Trigger", 1001, -0.5, 1000.5);
  hwrMasterTS = new TH1D("wrMasterTS", "WR master TS", 10000, 1582e15, 1583e15);
  hwrMasterTSspilloff = new TH1D("wrMasterTSspilloff", "WR master TS while spill-off", 10000, 1582e15, 1583e15);
  
  hhitmult = new TH1D("hitmult", "cluster multiplicity", 100, 0, 100);
  hmapmult = new TH1D("mapmult", "hit  multiplicity in map level;", 500, 0, 500);
  hcalmult = new TH1D("calmult", "hit  multiplicity in cal level;", 500, 0, 500);
  
  hmapEvscrystalID = new TH2D("mapEvscrystalID", "Energy in map level vs Crystal ID; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);
  hmapEvscrystalIDspillon = new TH2D("mapEvscrystalIDspillon", "Energy in map level vs Crystal ID for onspill; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);
  hmapEvscrystalIDspilloff = new TH2D("mapEvscrystalIDspilloff", "Energy in map level vs Crystal ID for offspill; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);
  hmapEvscrystalIDspillofftpat0 = new TH2D("mapEvscrystalIDspillofftpat0", "Energy in map level vs Crystal ID for offspill & tpat=0; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);
  hmapNfvscrystalID = new TH2D("mapNfvscrystalID", "Nf in map level vs Crystal ID; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);
  hmapNsvscrystalID = new TH2D("mapNsvscrystalID", "Ns in map level vs Crystal ID; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);
  hmapNfvscrystalIDspilloff = new TH2D("mapNfvscrystalIDspilloff", "Nf in map level vs Crystal ID while spill-off; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);
  hmapNsvscrystalIDspilloff = new TH2D("mapNsvscrystalIDspilloff", "Ns in map level vs Crystal ID while spill-off; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);
  hmapcrystalID = new TH1D("mapcrystalID", "Crystal ID; crystal ID", 5000, -0.5, 4999.5);

  hcalEvscrystalID = new TH2D("calEvscrystalID", "Energy in cal level vs Crystal ID; crystal ID;energy (keV))", 5000, -0.5, 4999.5, 3000,0,300000);
  hcalEvscrystalIDspillofftpat0 = new TH2D("calEvscrystalIDspillofftpat0", "Energy in cal level vs Crystal ID while spill-off & tpat=0; crystal ID;energy (keV))", 5000, -0.5, 4999.5, 3000,0,300000);
  hcalEspillofftpat0 = new TH1D("calEspillofftpat0", "Energy in cal level while spill-off & tpat=0;energy (keV))", 10000,0,10000);
  hcalNfvscrystalID = new TH2D("calNfvscrystalID", "Nf in cal level vs Crystal ID; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);
  hcalNsvscrystalID = new TH2D("calNsvscrystalID", "Ns in cal level vs Crystal ID; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);
  hcalNfvscrystalIDspilloff = new TH2D("calNfvscrystalIDspilloff", "Nf in cal level vs Crystal ID while spill-off; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);
  hcalNsvscrystalIDspilloff = new TH2D("calNsvscrystalIDspilloff", "Ns in cal level vs Crystal ID while spill-off; crystal ID;energy (arb. units))", 5000, -0.5, 4999.5, 1000,-50000,50000);

  hmapNfNs = new TH2D("mapNfNs", "QPID in map level: Nf vs Ns; Nf (arb. u.);Ns (arb. units))", 1000, -50000, 50000, 1000,-50000,50000);
  hmapNfNsspillon = new TH2D("mapNfNsspillon", "QPID in map level while spillon: Nf vs Ns; Nf (arb. u.);Ns (arb. units))", 1000, -50000, 50000, 1000,-50000,50000);
  hmapNfNsspilloff = new TH2D("mapNfNsspilloff", "QPID in map level while spilloff: Nf vs Ns; Nf (arb. u.);Ns (arb. units))", 1000, -50000, 50000, 1000,-50000,50000);
  hmapNfNsGamBar = new TH2D("mapNfNsGamBar", "QPID in map level for Barrel Gammas: Nf vs Ns; Nf (arb. u.);Ns (arb. units))", 1000, -40000, 40000, 1000,-40000,40000);
  hmapNfNsGamIph = new TH2D("mapNfNsGamIph", "QPID in map level for Iphos Gammas: Nf vs Ns; Nf (arb. u.);Ns (arb. units))", 1000, -40000, 40000, 1000,-40000,40000);
  hmapNfNsProBar = new TH2D("mapNfNsProBar", "QPID in map level for Barrel Protons: Nf vs Ns; Nf (arb. u.);Ns (arb. units))", 1000, -40000, 40000, 1000,-40000,40000);
  hmapNfNsProIph = new TH2D("mapNfNsProIph", "QPID in map level for Iphos Protons: Nf vs Ns; Nf (arb. u.);Ns (arb. units))", 1000, -40000, 40000, 1000,-40000,40000);
  hcalNfNs = new TH2D("calNfNs", "QPID in cal level: Nf vs Ns; Nf (keV);Ns (keV))", 1000, -100000, 100000, 1000,-100000,100000);
  hcalNfNsgoodtpat = new TH2D("calNfNsgoodtpat", "QPID in cal level for Tpat>0: Nf vs Ns; Nf (keV);Ns (keV))", 1000, -100000, 100000, 1000,-100000,100000);
  hcalNfNstpat0 = new TH2D("calNfNstpat0", "QPID in cal level for Tpat=0: Nf vs Ns; Nf (keV);Ns (keV))", 1000, -100000, 100000, 1000,-100000,100000);
  hcalNfNsGamBar = new TH2D("calNfNsGamBar", "QPID in cal level for Barrel Gammas: Nf vs Ns; Nf (keV);Ns (keV))", 1000, -50000, 50000, 1000,-40000,50000);
  hcalNfNsGamIph = new TH2D("calNfNsGamIph", "QPID in cal level for Iphos Gammas: Nf vs Ns; Nf (keV);Ns (keV))", 1000, -50000, 50000, 1000,-50000,50000);
  hcalNfNsProBar = new TH2D("calNfNsProBar", "QPID in cal level for Barrel Protons: Nf vs Ns; Nf (keV);Ns (keV)", 1000, -50000, 100000, 1000,-50000,50000);
  hcalNfNsProIph = new TH2D("calNfNsProIph", "QPID in cal level for Iphos Protons: Nf vs Ns; Nf (keV);Ns (keV)", 1000, -50000, 100000, 1000,-50000,50000);

}
void LA_califa_p2p(int runNum){
  //std::cout << preamp_ranges[1040] << "\n";
  gStyle->SetOptStat(111111);
  chain = new TChain("evt");
  DefineHistograms();
  AnalyzeRun(runNum);
  DrawResults();
  can->SaveAs("tpat.pdf");
}
