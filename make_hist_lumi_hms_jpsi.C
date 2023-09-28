#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void make_hist_lumi_hms_jpsi(TString basename="",Int_t nrun=2043,Double_t mean_current=2.0){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_lumi_hms_jpsi_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
TTree *tscal = (TTree*) fsimc->Get("TSH");
//
 Double_t  Scal_evNumber;
   tscal->SetBranchAddress("evNumber",&Scal_evNumber);
 Double_t  Scal_BCM4A_charge;
   tscal->SetBranchAddress("H.BCM4A.scalerCharge",&Scal_BCM4A_charge);
 Double_t  Scal_BCM4A_current;
   tscal->SetBranchAddress("H.BCM4A.scalerCurrent",&Scal_BCM4A_current);
 Double_t  Scal_time;
   tscal->SetBranchAddress("H.1MHz.scalerTime",&Scal_time);
 Double_t  Scal_TRIG2;
   tscal->SetBranchAddress("H.pTRIG2.scaler",&Scal_TRIG2);
 Double_t  Scal_TRIG3;
   tscal->SetBranchAddress("H.pTRIG3.scaler",&Scal_TRIG3);
 Double_t  Scal_TRIG1;
   tscal->SetBranchAddress("H.pTRIG1.scaler",&Scal_TRIG1);
 Double_t  Scal_TRIG4;
   tscal->SetBranchAddress("H.pTRIG4.scaler",&Scal_TRIG4);
 Double_t  Scal_TRIG5;
   tscal->SetBranchAddress("H.pTRIG5.scaler",&Scal_TRIG5);
 Double_t  Scal_Splane[4];
   tscal->SetBranchAddress("H.S1X.scaler",&Scal_Splane[0]);
   tscal->SetBranchAddress("H.S1Y.scaler",&Scal_Splane[1]);
   tscal->SetBranchAddress("H.S2X.scaler",&Scal_Splane[2]);
   tscal->SetBranchAddress("H.S2Y.scaler",&Scal_Splane[3]);
//loop through scalers
    Double_t tot_scal_Splane[4]={0,0,0,0};
     Double_t prev_Splane[4]={0,0,0,0};
    Double_t tot_scal_cut_Splane[4]={0,0,0,0};
     Int_t nscal_reads=0;
     Int_t nscal_reads_cut=0;
     Double_t prev_read=-1;
     Double_t ave_current=0;
     Double_t ave_current_cut=0;
     Double_t charge_sum=0;
     Double_t charge_sum_cut=0;
     Double_t prev_charge=0;
     Double_t event_flag[10000];
     Double_t scal_event_number[10000];
     Double_t tot_scal_TRIG2=0;
     Double_t tot_scal_TRIG3=0;
     Double_t prev_TRIG2=0;
     Double_t prev_TRIG3=0;
     Double_t tot_scal_cut_TRIG2=0;
     Double_t tot_scal_cut_TRIG3=0;
     Double_t tot_scal_TRIG1=0;
     Double_t tot_scal_TRIG4=0;
     Double_t prev_TRIG1=0;
     Double_t prev_TRIG4=0;
     Double_t tot_scal_cut_TRIG1=0;
     Double_t tot_scal_cut_TRIG4=0;
     Double_t tot_scal_cut_time=0;
     Double_t tot_scal_TRIG5=0;
     Double_t prev_TRIG5=0;
    Double_t tot_scal_cut_TRIG5=0;
      Double_t tot_scal_time=0;
     Double_t prev_time=0;
     //
   Double_t threshold_cut=3.;
Long64_t scal_entries = tscal->GetEntries();
Long64_t data_entries = tsimc->GetEntries();
	for (int i = 0; i < scal_entries; i++) {
      		tscal->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
          event_flag[nscal_reads] = 0;
             scal_event_number[nscal_reads] = Scal_evNumber;
          ave_current+=Scal_BCM4A_current;
	  if (TMath::Abs(Scal_BCM4A_current-mean_current) < threshold_cut) {
             event_flag[nscal_reads] = 1;
             ave_current_cut+=Scal_BCM4A_current;
 	     tot_scal_cut_time+=(Scal_time-prev_time);
 	     tot_scal_cut_TRIG2+=(Scal_TRIG2-prev_TRIG2);
 	     tot_scal_cut_TRIG3+=(Scal_TRIG3-prev_TRIG3);
 	     tot_scal_cut_TRIG1+=(Scal_TRIG1-prev_TRIG1);
 	     tot_scal_cut_TRIG4+=(Scal_TRIG4-prev_TRIG4);
 	     tot_scal_cut_TRIG5+=(Scal_TRIG5-prev_TRIG5);
 	     for (Int_t s=0;s<4;s++) tot_scal_cut_Splane[s]+=(Scal_Splane[s]-prev_Splane[s]);
            charge_sum_cut+=(Scal_BCM4A_charge-prev_charge);
             nscal_reads_cut++;
	  }
	  prev_charge = Scal_BCM4A_charge;
	  prev_time = Scal_time;
	  prev_TRIG2 = Scal_TRIG2;
	  prev_TRIG3 = Scal_TRIG3;
	  prev_TRIG1 = Scal_TRIG1;
	  prev_TRIG4 = Scal_TRIG4;
	  prev_TRIG5 = Scal_TRIG5;
	     for (Int_t s=0;s<4;s++) prev_Splane[s]=Scal_Splane[s];
	  // cout <<  nscal_reads <<  " " << Scal_BCM4A_current << " " << event_flag[nscal_reads] << " " << Scal_TRIG1 << endl;
          nscal_reads++;
          charge_sum=Scal_BCM4A_charge;
	  tot_scal_TRIG2=Scal_TRIG2;
	  tot_scal_TRIG3=Scal_TRIG3;
	  tot_scal_TRIG1=Scal_TRIG1;
	  tot_scal_TRIG4=Scal_TRIG4;
	  tot_scal_TRIG5=Scal_TRIG5;
          tot_scal_time=Scal_time;
	  for (Int_t s=0;s<4;s++) tot_scal_Splane[s]=Scal_Splane[s];
	}
   //
	Double_t gevtyp;
	tsimc->SetBranchAddress("g.evtyp",&gevtyp);
	Double_t gevnum;
	tsimc->SetBranchAddress("g.evnum",&gevnum);
	Double_t ntrack;
	tsimc->SetBranchAddress("H.dc.ntrack",&ntrack);
	Double_t goodscinhit;
	tsimc->SetBranchAddress("H.hod.goodscinhit",&goodscinhit);
	Double_t starttime;
	tsimc->SetBranchAddress("H.hod.starttime",&starttime);
 	Double_t betanotrack;
	tsimc->SetBranchAddress("H.hod.betanotrack",&betanotrack);
   Double_t hgcer_npeSum;
   tsimc->SetBranchAddress("H.cer.npeSum",&hgcer_npeSum);
   Double_t delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
   Double_t etotnorm;
   tsimc->SetBranchAddress("H.cal.etotnorm",&etotnorm);
   Double_t etottracknorm;
   tsimc->SetBranchAddress("H.cal.etottracknorm",&etottracknorm);
	//
	TH1F* h_starttime = new TH1F("h_starttime"," ; Startime (elec good scinhit)",60,0.,120.);
	TH1F* h_ch1_nhit = new TH1F("h_ch1_nhit"," ; CH1 nhits",30,0.,30.);
	TH1F* h_ch2_nhit = new TH1F("h_ch2_nhit"," ; CH2 nhits",30,0.,30.);
	TH1F* h_starttime_ntrack = new TH1F("h_starttime_ntrack"," ; Startime (elec good scinhit ntrack)",60,0.,120.);
	TH1F* h_etotnorm_all = new TH1F("h_etotnorm_all"," All events; etot  norm",100,0.,3.);
	TH1F* h_etottracknorm_trig2 = new TH1F("h_etottracknorm_trig1"," Trig 2 ; etot track norm",100,0.,3.);
	TH1F* h_etottracknorm_trig4 = new TH1F("h_etottracknorm_trig4"," Trig 4 ; etot track norm",100,0.,3.);
	TH1F* h_etottracknorm_all = new TH1F("h_etottracknorm_all"," All events; etot track norm",100,0.,3.);
	TH1F* h_hgcernpeSum_all = new TH1F("h_ngcernpeSum_all"," All events; HG cer npe sum",120,0.,30.);
	TH1F* h_hgcernpeSum_all_etotcut = new TH1F("h_ngcernpeSum_all_etotcut"," All events; HG cer npe sum",120,0.,30.);
	TH1F* h_etottracknorm_curcut = new TH1F("h_etottracknorm_curcut"," Current cut; etot track norm",100,0.,3.);
	//
   Int_t nscal_reads_2=0;
   prev_read=-1;
   Int_t ps3;
   Double_t trackeff;
   if (nrun==7570) ps3=17;
   if (nrun==7571) ps3=3;
   if (nrun==7572) ps3=5;
   if (nrun==7574) ps3=17;
   if (nrun==7499) ps3=2;
   if (nrun==7500) ps3=2;
   if (nrun==7501) ps3=2;
   if (nrun==7419) ps3=2;
   if (nrun==7420) ps3=2;
   //Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < data_entries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		h_hgcernpeSum_all->Fill(hgcer_npeSum);
		if (etottracknorm > 0.97) h_hgcernpeSum_all_etotcut->Fill(hgcer_npeSum);
		if (hgcer_npeSum>1&& delta>-8&& delta<8) h_etottracknorm_all->Fill(etottracknorm);
 		if (hgcer_npeSum>1&& delta>-8&& delta<8) h_etotnorm_all->Fill(etottracknorm);
		if (event_flag[nscal_reads_2]==1&&etotnorm>0.97&&hgcer_npeSum>1&&goodscinhit==1) {
		  h_starttime->Fill(starttime);
		  if (ntrack >0) h_starttime_ntrack->Fill(starttime);
		}
                if (event_flag[nscal_reads_2]==1&&etottracknorm>0.97&&hgcer_npeSum>1.&& delta>-8&& delta<8) {
		  h_etottracknorm_curcut->Fill(etottracknorm);
		  if (gevtyp==2) {
		    h_etottracknorm_trig2->Fill(etottracknorm);
		      }
		  if (gevtyp==4) {
		    h_etottracknorm_trig4->Fill(etottracknorm);
		      }
 		}
	       	//  cout << nscal_reads_2 << " " << gevtyp << " " << gevnum << " " << scal_event_number[nscal_reads_2] << endl;
		
		if (gevnum>scal_event_number[nscal_reads_2]) {
                   nscal_reads_2++;
		}
	}
	//
	Double_t calc_treff= float(h_starttime_ntrack->Integral())/float(h_starttime->Integral()) ;
	Int_t good_ev2 = h_etottracknorm_trig2->Integral();
	Double_t good_ev2_err = TMath::Sqrt(good_ev2);
	Int_t good_ev4 = h_etottracknorm_trig4->Integral();
	Double_t good_ev4_err = TMath::Sqrt(good_ev4);
	Double_t good_ev = (good_ev2*ps3+good_ev4)/calc_treff;
	Double_t good_ev_err = TMath::Sqrt( (good_ev2_err*ps3)*(good_ev2_err*ps3) + (good_ev4_err)*(good_ev4_err))/calc_treff;
	//
        Double_t err_trig2 = 1./TMath::Sqrt(tot_scal_cut_TRIG2);
        Double_t err_trig3 = 1./TMath::Sqrt(tot_scal_cut_TRIG3);
        Double_t err_trig1 = 1./TMath::Sqrt(tot_scal_cut_TRIG1);
        Double_t err_trig4 = 1./TMath::Sqrt(tot_scal_cut_TRIG4);
        Double_t err_trig5 = 1./TMath::Sqrt(tot_scal_cut_TRIG5);
	cout << nrun << " " << charge_sum_cut << " " << threshold_cut << " "<< charge_sum_cut/tot_scal_cut_time << " " << tot_scal_cut_time/tot_scal_time << " " << " " << tot_scal_cut_TRIG3/tot_scal_cut_time << " " << tot_scal_cut_TRIG4/tot_scal_cut_time << " " << good_ev/tot_scal_cut_time  << " " << tot_scal_cut_TRIG3/charge_sum_cut << " " << err_trig3*tot_scal_cut_TRIG3/charge_sum_cut << " " << tot_scal_cut_TRIG4/charge_sum_cut << " " << err_trig4*tot_scal_cut_TRIG4/charge_sum_cut<< " " << good_ev/charge_sum_cut<< " " << good_ev_err/charge_sum_cut << " " << endl;
cout << nrun << " " << mean_current << " "<< charge_sum_cut/tot_scal_cut_time << " "<< tot_scal_cut_Splane[0]/tot_scal_cut_time<< " "<< tot_scal_cut_Splane[1]/tot_scal_cut_time<< " "<< tot_scal_cut_Splane[2]/tot_scal_cut_time<< " "<< tot_scal_cut_Splane[3]/tot_scal_cut_time<< " " << calc_treff<< endl;
// loop through data
//
}
