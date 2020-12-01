#include "TCanvas.h"
#include "TProof.h"
#include "TString.h"
#include "TChain.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TLine.h"
#include <fstream>

using namespace std;

//TProof *p =  TProof::Open("");
TProof *p =  TProof::Open("workers=10");

TCanvas *c;
TString Zgate = "abs(MusicZ-TwimZ)<0.5";
TString fragbeta[NUMPADDLE] = {""};
TString fragaoq[NUMPADDLE] = {""};
Double_t mc_e = 3.1071; // m_u * c_0 / e
Double_t min_brho = 8.7, max_brho=9.3;
TString brho = Form("Beta_S2_Cave / sqrt(1-Beta_S2_Cave*Beta_S2_Cave) * AoQ_S2_Cave * %f", mc_e);
TString brhocave = "";
TChain *ch;
Double_t par[NUMPADDLE][2] = {{0.}};

TH2D *h_beta_tof[NUMPADDLE], *h_beta_tof_cut[NUMPADDLE], *h_music_twim[2], *h_beta_beta[NUMPADDLE], *h_beta_beta_cut[NUMPADDLE];
TH2D *h_mw_brho[NUMPADDLE+1], *h_fragaoq[NUMPADDLE+1], *h_fragpid[NUMPADDLE+1], *h_frspid;
TProfile *prof_tof[NUMPADDLE], *prof_mw_brho;
TH1D *h_tofw_paddle, *h_median_tof[NUMPADDLE];
TF1 *f_tof[NUMPADDLE], *f_mw3brho;

//For MWPC ana
#define NUMCOND 20
int cond = 0;
TString beamcondition = "";//"abs(MusicZ-20.)<0.4 && abs(Beta_S2_Cave - 0.765)<0.005 &&";
TString conditionwithbetacut = "";
TString mwcondition = "";
TString beta_tofw="";
TString beta_tofw_mod="";
TString fragbrhostring="";
TString fragaoqstring="";
//
Double_t max_mw1[4][2]={{0.}}, range_mw = 60.;
TLine *line_mw[4][2][2];
TH2D *h_frspid_mw[NUMCOND], *h_music_twim_mw[NUMCOND], *h_beta_beta_mw[NUMCOND][NUMPADDLE+1];
TString axis_mw_h[4] = {"Mw2_X", "Mw2_X", "Mw2_Y", "(Mw2_X-Mw1_X)"};
TString axis_mw_v[4] = {"Mw2_Y", "(Mw2_X-Mw1_X)", "(Mw2_Y-Mw1_Y)", "(Mw2_Y-Mw1_Y)"};
TString axis_mw3[4] = {"Mw2_X", "Mw1_Y", "(Mw2_X-Mw1_X)", "(Mw2_Y-Mw1_Y)"};
//TString axis_mw_v[4] = {"Mw1_Y", "(Mw2_X-Mw1_X)", "(Mw3_Y-Mw1_Y)", "(Mw3_Y-Mw1_Y)"};
//TString axis_mw3[4] = {"Mw1_X", "Mw1_Y", "(Mw2_X-Mw1_X)", "(Mw3_Y-Mw1_Y)"};
TString Mw3_X_mod = "Mw3_X";
TH2D *h_mw12[NUMCOND][4], *h_mw3[NUMCOND][4], *h_mwbeta[NUMCOND][4];
TProfile *prof_mw3[4], *prof_mwbeta[NUMCOND][4];
TH1D *h_median_mwbeta[NUMCOND][4];
TF1 *f_mw3[NUMCOND][4], *f_mwbeta[NUMCOND][4];
Double_t range_fit_mw3_low[4]={-25., -10., -8., -50.};
Double_t range_fit_mw3_high[4]={15., 20., 0., -35.};
Double_t range_cut_mw3_low[4]={-50., -60., -30., -60.};
Double_t range_cut_mw3_high[4]={50., 60., 30., -10.};
//Double_t range_cut_mw3_low[4]={-30., -60., -30., -60.};
//Double_t range_cut_mw3_high[4]={25., 60., 30., -10.};
//TCut cut_mw[4];
TString cut_mw = "";

//For Brho reconstruction
TH2D *h_beta_mw3[NUMCOND][4], *h_brho_mw3[NUMCOND][4], *h_brhobrho, *h_aoqaoq, *h_pid;
TProfile *prof_beta_mw3[NUMCOND][4], *prof_brho_mw3[NUMCOND][4];
TF1 *f_beta_mw3[NUMCOND][4], *f_brho_mw3[NUMCOND][4];
