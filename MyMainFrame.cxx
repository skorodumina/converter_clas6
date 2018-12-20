#include "TROOT.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include <math.h>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TText.h"
#include "TStyle.h"
#include "TGObject.h"
#include "TObject.h"
#include "TSystem.h"
#include "TMinuit.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGButtonGroup.h>
#include <RQ_OBJECT.h>
#include <TGNumberEntry.h>
#include <TGProgressBar.h>
#include <TGLabel.h>
#include <stdio.h>
#include <dlfcn.h>
#include "MyMainFrame.h"
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <TGFileDialog.h>
#include <GuiTypes.h>
#include <TGDoubleSlider.h>
#include <TGComboBox.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TError.h> 
#include <auto_ptr.h>
#ifndef __CINT__
#include <cstdlib>
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include <RooLandau.h> 
#include <RooNumConvPdf.h>
#include <RooDataHist.h>
#include "RooBinning.h"
#include <sys/types.h>
#include <wait.h>
#include <unistd.h>
#include <cstring>
#include <getopt.h>
#include <cstdlib>



 using namespace std; 
 
#define _USE_MATH_DEFINES
 
   




    void MyMainFrame::MainFrame(UChar_t flag, Float_t E_beam, Short_t nfiles, string inp_files[], string outfile_in) { 


	inpfile_inp = inp_files[0];
	outfile_inp = outfile_in;
	
	E0 = E_beam;

        n_files = nfiles;

        file = new string[nfiles];

        file = inp_files;
	data_sim = flag;
	

DoDraw();


   
}

void MyMainFrame::DoDraw() {

 
 
       t20tot21(); 
   
}








    void MyMainFrame::t20tot21() { 
    

    ostringstream adc_num;
    ostringstream tdc_num;
    ostringstream ref_tdc;
    
    
	 
    Long64_t j;
    
 Double_t * adc_offset;
 adc_offset = new Double_t [12];
 Double_t * adc_cut;
 adc_cut = new Double_t [12]; 
 Short_t m,l;
 Long64_t i,nstart,nstop,n_incl,n_elast;

 
TH1I *hist_adc_off[12]; 

 TFile *finp;
 Int_t block_total = 0;
 Int_t block_last = 0;
 Float_t Qfull = 0.;
 
 
 
 
 
 
 
 
 
 
  	t21 = new TTree("t21","Tree 21"); 
  t21->SetDirectory(0);
  
  
	Int_t         npart,segment,sector,indtype,pdhit,n_PIp,n_PIm,n_P;
	Short_t pmt_hit;
	Int_t PdHit_EL,PdHit_PIp,PdHit_PIm,PdHit_P;
	Float_t * p; 
        p = new Float_t [20];
        Float_t  P_EL,P_EL_new,ph_EL,th_EL,th_EL_new,W,Q2,m_proton,NpheCC_EL,ECtot_EL,hadr_mom;
	Float_t  x_EL,y_EL,z_EL,pion_low,pion_high,proton_low,proton_high;
	Float_t  z_PIp,z_PIm,z_P;
	Float_t  dc_x_EL,dc_y_EL,dc_z_EL;
	Float_t  dc_x_P,dc_y_P,dc_z_P;
	Float_t  dc_x_PIp,dc_y_PIp,dc_z_PIp;
	Float_t  dc_x_PIm,dc_y_PIm,dc_z_PIm;
	Float_t  sigma,px_fermi,py_fermi,pz_fermi;
	Float_t  ECin_EL,ECout_EL,P_PIp,P_PIm,P_P,beta,th_PIp,th_PIm,th_P;
	Int_t block;
        Long64_t gpart,k,last_i,last_k;
	Float_t  q_l,Qdiff,Qcurr,Qprev,Qtotal,deltaQ,ph_PIp,ph_P,ph_PIm;
	Int_t sc_part_local;
	Float_t sc_pd_local,sc_sect_local;
	Float_t sc_x,sc_y,sc_z;
	//Float_t sc_x_pip,sc_y_pip,sc_z_pip;
	//Float_t sc_x_pim,sc_y_pim,sc_z_pim;
	//Float_t sc_x_p,sc_y_p,sc_z_p;
	Float_t nx,ny,nz,t,theta_cc;
	Float_t sx,sy,sz,px,py,pz,ph_cc,LiveTime;

	Float_t beta_PIm,beta_PIp,beta_P,NpheCC_PIp;
	Float_t beta_PIp_time,beta_PIm_time,beta_P_time;
	Float_t EL_dist,PIp_dist,PIm_dist,P_dist;
	Float_t PIp_time,PIm_time,P_time, c;
	Float_t EL_tof,PIp_tof,P_tof,PIm_tof;
	Float_t pf_x,pf_y,pf_z;
	
	sx=-0.000784;
	sy=0.;
	sz=-0.00168;
	c =30.;	
	
//	E0 = 2.039;
	m_proton = 0.93821;
  //ob'yavlyayutsya vetki v novom dereve t21 
t21->Branch("indtype",&indtype);  
t21->Branch("npart",&npart); 
t21->Branch("P_EL",&P_EL);
t21->Branch("block",&block_total);
t21->Branch("deltaQ",&deltaQ);
t21->Branch("LiveTime",&LiveTime);
t21->Branch("n_incl",&n_incl);
t21->Branch("n_elast",&n_elast);
t21->Branch("ph_EL",&ph_EL);
t21->Branch("th_EL",&th_EL);
t21->Branch("W",&W);
t21->Branch("Q2",&Q2);
t21->Branch("NpheCC_EL",&NpheCC_EL);
t21->Branch("ECtot_EL",&ECtot_EL);
t21->Branch("ECin_EL",&ECin_EL);
t21->Branch("ECout_EL",&ECout_EL);
t21->Branch("x_EL",&x_EL);
t21->Branch("y_EL",&y_EL);
t21->Branch("z_EL",&z_EL);
t21->Branch("z_P",&z_P);
t21->Branch("z_PIp",&z_PIp);
t21->Branch("z_PIm",&z_PIm);

t21->Branch("dc_x_EL",&dc_x_EL);
t21->Branch("dc_y_EL",&dc_y_EL);
t21->Branch("dc_z_EL",&dc_z_EL);

t21->Branch("dc_x_P",&dc_x_P);
t21->Branch("dc_y_P",&dc_y_P);
t21->Branch("dc_z_P",&dc_z_P);

t21->Branch("dc_x_PIp",&dc_x_PIp);
t21->Branch("dc_y_PIp",&dc_y_PIp);
t21->Branch("dc_z_PIp",&dc_z_PIp);

t21->Branch("dc_x_PIm",&dc_x_PIm);
t21->Branch("dc_y_PIm",&dc_y_PIm);
t21->Branch("dc_z_PIm",&dc_z_PIm);

t21->Branch("PdHit_EL",&pdhit);
t21->Branch("PdHit_PIp",&PdHit_PIp);
t21->Branch("PdHit_PIm",&PdHit_PIm);
t21->Branch("PdHit_P",&PdHit_P);
t21->Branch("sc_x",&sc_x);
t21->Branch("sc_y",&sc_y);
t21->Branch("sc_z",&sc_z);
//t21->Branch("sc_x_pip",&sc_x_pip);
//t21->Branch("sc_y_pip",&sc_y_pip);
//t21->Branch("sc_z_pip",&sc_z_pip);
//t21->Branch("sc_x_p",&sc_x_p);
//t21->Branch("sc_y_p",&sc_y_p);
//t21->Branch("sc_z_p",&sc_z_p);
//t21->Branch("sc_x_pim",&sc_x_pim);
//t21->Branch("sc_y_pim",&sc_y_pim);
//t21->Branch("sc_z_pim",&sc_z_pim);


t21->Branch("pmt_hit",&pmt_hit);
t21->Branch("segment",&segment);
t21->Branch("theta_cc",&theta_cc);
t21->Branch("ph_cc",&ph_cc);
t21->Branch("sector",&sector); 
t21->Branch("n_PIp",&n_PIp);
t21->Branch("n_PIm",&n_PIm);
t21->Branch("n_P",&n_P); 
t21->Branch("P_PIp",&P_PIp,"P_PIp/F");
t21->Branch("P_PIm",&P_PIm,"P_PIm/F");
t21->Branch("P_P",&P_P,"P_P/F");  
t21->Branch("th_PIp",&th_PIp,"th_PIp/F");
t21->Branch("th_PIm",&th_PIm,"th_PIm/F");
t21->Branch("th_P",&th_P,"th_P/F"); 
t21->Branch("ph_PIp",&ph_PIp,"ph_PIp/F");
t21->Branch("ph_PIm",&ph_PIm,"ph_PIm/F");
t21->Branch("ph_P",&ph_P,"ph_P/F");  
t21->Branch("beta_PIm",&beta_PIm,"beta_PIm/F");
t21->Branch("beta_PIp",&beta_PIp,"beta_PIp/F");  
t21->Branch("beta_P",&beta_P,"beta_P/F"); 
t21->Branch("NpheCC_PIp",&NpheCC_PIp);

t21->Branch("PIp_time",&PIp_time,"PIp_time/F"); 
t21->Branch("PIp_dist",&PIp_dist,"PIp_dist/F"); 
t21->Branch("beta_PIp_time",&beta_PIp_time,"beta_PIp_time/F"); 

t21->Branch("P_time",&P_time,"P_time/F"); 
t21->Branch("P_dist",&P_dist,"P_dist/F"); 
t21->Branch("beta_P_time",&beta_P_time,"beta_P_time/F"); 

t21->Branch("PIm_time",&PIm_time,"PIm_time/F"); 
t21->Branch("PIm_dist",&PIm_dist,"PIm_dist/F");
// t21->Branch("pf_x",&pf_x,"pf_x/F");
//  t21->Branch("pf_y",&pf_y,"pf_y/F");
// t21->Branch("pf_z",&pf_z,"pf_z/F"); 

t21->Branch("beta_PIm_time",&beta_PIm_time,"beta_PIm_time/F");
//t21->Branch("tr_time",&tr_time,"tr_time/F");
 
// t21->Branch("EL_time_fr_vert",&EL_time_fr_vert,"EL_time_fr_vert/F");
// t21->Branch("EL_tof",&EL_tof,"EL_tof/F");
 
 
 t21->Branch("sigma",&sigma,"sigma/F");
 t21->Branch("px_fermi",&px_fermi,"px_fermi/F");
 t21->Branch("py_fermi",&py_fermi,"py_fermi/F");
 t21->Branch("pz_fermi",&pz_fermi,"pz_fermi/F");

 
 if (data_sim == 1) {
 indtype = 1;
 };
 
//tsikl po failam
  for (m=1; m<=n_files; m++) {
  
  finp = new TFile(file[m-1].c_str()); 

  
  TTree *t20 = (TTree*)finp->Get("h10");
  TBranch *br_npart = t20->GetBranch("npart");
  TBranch *br_gpart = t20->GetBranch("gpart");
   TBranch *br_stat = t20->GetBranch("stat");
  TBranch *br_b = t20->GetBranch("b");
  TBranch *br_p = t20->GetBranch("p");
  TBranch *br_q = t20->GetBranch("q");
   TBranch *br_dc = t20->GetBranch("dc");
   TBranch *br_dc_stat = t20->GetBranch("dc_stat"); 
   TBranch *br_sc = t20->GetBranch("sc");
   TBranch *br_ec = t20->GetBranch("ec");
   TBranch *br_cc = t20->GetBranch("cc");
   TBranch *br_cx = t20->GetBranch("cx");
   TBranch *br_cy = t20->GetBranch("cy");
   TBranch *br_cz = t20->GetBranch("cz"); 
   TBranch *br_nphe = t20->GetBranch("nphe"); 
   TBranch *br_etot = t20->GetBranch("etot");
   TBranch *br_ec_ei = t20->GetBranch("ec_ei");
   TBranch *br_ec_eo = t20->GetBranch("ec_eo");
   TBranch *br_vx = t20->GetBranch("vx");
   TBranch *br_vy = t20->GetBranch("vy");
   TBranch *br_vz = t20->GetBranch("vz"); 
   TBranch *br_sc_sect = t20->GetBranch("sc_sect");
   TBranch *br_sc_pd = t20->GetBranch("sc_pd"); 
   TBranch *br_sc_part = t20->GetBranch("sc_part"); 
   TBranch *br_dc_xsc = t20->GetBranch("dc_xsc"); 
   TBranch *br_dc_ysc = t20->GetBranch("dc_ysc");
   TBranch *br_dc_zsc = t20->GetBranch("dc_zsc");  
   TBranch *br_cc_segm = t20->GetBranch("cc_segm");  
   TBranch *br_dc_cxsc = t20->GetBranch("dc_cxsc"); 
   TBranch *br_dc_cysc = t20->GetBranch("dc_cysc");
   TBranch *br_dc_czsc = t20->GetBranch("dc_czsc");     
          
   
   TBranch *br_q_l = t20->GetBranch("q_l");
   TBranch *br_t_l = t20->GetBranch("t_l"); 
   TBranch *br_sc_r = t20->GetBranch("sc_r");
   TBranch *br_sc_t = t20->GetBranch("sc_t");
   
   
TBranch *br_dc_vx;
   TBranch *br_dc_vy;
   TBranch *br_dc_vz;

TBranch *br_mcnentr;
TBranch *br_mcnpart;
TBranch *br_mcst;
TBranch *br_mcid;
TBranch *br_mcpid;
TBranch *br_mctheta;
TBranch *br_mcphi;
TBranch *br_mcp;
TBranch *br_mcm;
TBranch *br_mcvx_x_el;
TBranch *br_mcvx_y_el;
TBranch *br_mcvx_z_el;
TBranch *br_sigma_total;
TBranch *br_px_Fermi;
TBranch *br_py_Fermi;
TBranch *br_pz_Fermi;

//TBranch *br_pf_x;
//TBranch *br_pf_y;
//TBranch *br_pf_z;

if (data_sim == 1) {

  br_dc_vx = t20->GetBranch("dc_vx");
  br_dc_vy = t20->GetBranch("dc_vy");
  br_dc_vz = t20->GetBranch("dc_vz");


};



 if (data_sim == 2) { //to esti tol'ko dlya simulyatsii
br_mcnentr = t20->GetBranch("mcnentr");
br_mcnpart = t20->GetBranch("mcnpart");
br_mcst = t20->GetBranch("mcst");
br_mcid = t20->GetBranch("mcid");
br_mcpid = t20->GetBranch("mcpid");
br_mctheta = t20->GetBranch("mctheta");
br_mcphi = t20->GetBranch("mcphi");
br_mcp = t20->GetBranch("mcp");
br_mcm= t20->GetBranch("mcm");
br_mcvx_x_el= t20->GetBranch("mcvx_x_el");
br_mcvx_y_el= t20->GetBranch("mcvx_y_el");
br_mcvx_z_el= t20->GetBranch("mcvx_z_el");
br_sigma_total = t20->GetBranch("sigma_total");
br_px_Fermi = t20->GetBranch("px_Fermi");
br_py_Fermi = t20->GetBranch("py_Fermi");
br_pz_Fermi = t20->GetBranch("pz_Fermi");


//br_pf_x= t20->GetBranch("pf_x");
//br_pf_y= t20->GetBranch("pf_y");
//br_pf_z= t20->GetBranch("pf_z");
 };
                         
//  cout << "E0 = " << E0 << "\n";
  
  cout << "Processing file " << m << "  N enteries = " << br_npart->GetEntries() << "\n";
  adc_num << "adc" << j+1;
    
  // Creating and filling new tree
  



// Loop over old tree 

  Bool_t adc_cut_switch,tdc_cut_switch;
  
  Qdiff = 0.;
  Qcurr = 0.;
  Qprev = 0.;
  Qtotal = 0.;
  k = 0;
  block = 0;
  last_k = 0;
  nstart = 0;
  nstop = 0;
  n_incl = 0;
  n_elast = 0;
  
  
//  cout << M_PI << "\n";
// cout << "npart = " << br_npart->GetEntries() << "\n";

//tsikl po sobitiyam/ mogno ispol'zovat' lyubuyu peremennuyu, tk v nih odinakovoe kol-vo sobitij 
  for (i=0; i<br_npart->GetEntries(); i++) { 
  
  Qprev = Qcurr;



  br_npart-> GetEntry(i);
  br_gpart-> GetEntry(i);
  br_stat->GetEntry(i); 
  br_p-> GetEntry(i);
  br_b-> GetEntry(i);
  br_q-> GetEntry(i);
  br_dc-> GetEntry(i);
  br_sc-> GetEntry(i);
  br_cc-> GetEntry(i);
  br_ec-> GetEntry(i);
  br_cx-> GetEntry(i);
  br_cy-> GetEntry(i);
  br_cz-> GetEntry(i);  
  br_dc_stat-> GetEntry(i);  
  br_q_l-> GetEntry(i); 
  //br_tr_time-> GetEntry(i);
  br_t_l-> GetEntry(i); 
  br_nphe-> GetEntry(i); 
  br_etot-> GetEntry(i); 
  br_ec_ei-> GetEntry(i);
  br_ec_eo-> GetEntry(i);  
  br_vx-> GetEntry(i);
  br_vy-> GetEntry(i);
  br_vz-> GetEntry(i);
  br_sc_sect-> GetEntry(i);
  br_sc_part-> GetEntry(i);  
  br_sc_pd-> GetEntry(i);    
  br_dc_xsc-> GetEntry(i);
  br_dc_ysc-> GetEntry(i);   
  br_dc_zsc-> GetEntry(i);  
  br_cc_segm-> GetEntry(i);
  br_dc_cxsc-> GetEntry(i);
  br_dc_cysc-> GetEntry(i);   
  br_dc_czsc-> GetEntry(i);   
  br_sc_t-> GetEntry(i);
  br_sc_r-> GetEntry(i);

 if (data_sim==2){
  
 br_mcnentr-> GetEntry(i); 
br_mcnpart-> GetEntry(i);
br_mcst-> GetEntry(i);
br_mcid-> GetEntry(i);
br_mcpid-> GetEntry(i);
br_mctheta-> GetEntry(i); 
br_mcphi-> GetEntry(i);
br_mcp -> GetEntry(i);
br_mcm-> GetEntry(i);
br_mcvx_x_el-> GetEntry(i);
br_mcvx_y_el-> GetEntry(i);
br_mcvx_z_el-> GetEntry(i);
br_sigma_total->GetEntry(i);
br_px_Fermi->GetEntry(i);
br_py_Fermi->GetEntry(i);
br_pz_Fermi->GetEntry(i);


//br_pf_x-> GetEntry(i);
//br_pf_y-> GetEntry(i);
//br_pf_z-> GetEntry(i);
Qdiff = 1.;
indtype = 1.;
}; 
//------------------------
if (data_sim==1){
 //nazivaem pervuyu chastitsu electronom (predvavritel'no) - chtobi podschitat' n_incl & n_elast
 if (br_gpart->GetLeaf("gpart")->GetValue() > 0){
    if (br_stat->GetLeaf("stat")->GetValue(0) > 0){
     if (br_dc->GetLeaf("dc")->GetValue(0) > 0) {
    if (br_sc->GetLeaf("sc")->GetValue(0) > 0) {
    if (br_ec->GetLeaf("ec")->GetValue(0) > 0) {
    if (br_cc->GetLeaf("cc")->GetValue(0) > 0) {  
    if (br_q->GetLeaf("q")->GetValue(0) == -1) {    
    
   n_incl = n_incl + 1;   

P_EL = br_p->GetLeaf("p")->GetValue(0);  
th_EL = (180./M_PI)*acos(br_cz->GetLeaf("cz")->GetValue(0));
		
Q2 = 4.*E0*P_EL*(sin(th_EL*M_PI/2./180.))*(sin(th_EL*M_PI/2./180.));
W = m_proton*m_proton+2.*m_proton*(E0-P_EL)-Q2;	
if (W > 0.) { 
W = sqrt(W);
if ( (W > 0.85) && (W < 1.15) ) {

n_elast = n_elast + 1;

};
};    
    
    };
    };
    };
    };
    };
  };
  };//end of n_incl & n_elast selection
  
  Qcurr = br_q_l->GetLeaf("q_l")->GetValue();
  Qdiff = Qcurr - Qprev;

  
  if (block == 0) {
  n_incl = 0;
  n_elast = 0;
  };
  };
 //---------------------
 
  
  nstart = nstop;
  if (Qdiff > 0.) {
  nstop = i;
  block = block + 1;
  
 if (data_sim==2){
 block = 2;
 nstart = 1;
 nstop = 2;
 };
     if (block > 1) {
   deltaQ = Qdiff;
//    cout << "deltaQ = " << deltaQ << "\n";
//nachinaem vsegda so 2go blocka
    for (j=nstart; j<nstop; j++) {
// mi snachala opredeliaem gde block konchaetsia (t.e. kak zariad pomenialsia tak novi' block nachalsia), a potom kogda nam razmer blocka uze izvesten mi delaem loop po sobitiam vnutri blocka uze.    


//Mi kazdomu sobitiu prisvaevaem nomer blocka, chtobi potom delat' cuti na live time, chitat' zariad i t.d. A poka mi ne doshli do mesta gde zariad pomenialsia nam nomer bloka neizvesten, poetomu kak on nam stal izvesten mi vozvraschaemsia v ego nachalo i pishem vse sobitia uze s nomerom bloka i deltaQ dlia nih
//pervi' raz mi schitaem vihod inclusivnih i uprugih sobiti', a dlia etogo nam elektron nuzen. Potom mi etot vihod toze kazdomu bloku svoi' prisvaevaem...
//vo vneshnem mi nahodim gde nachinaetsia i konchaetsia blok, zatem vozvraschaemsia k nachalu bloka i chitaem tam vse chastici, vrode vse ok


// cikl po j sluchaetsia tol'ko esli nstart ne ravno nstop, t.e. tol'ko kogda blok pomenialsia
//dlia simuliacii blocki voobsche ne nuzni, tam takogo poniatia voobsche net


//nstart i nstop eto nomer pervogo i poslednego sobitya v bloke

//dlia simuliacii vetki berutsia tol'ko v cikle po i
//dal'se vsegda dlia simuliacii cikl po j vipolniaetsia odin raz t.k. tam vsegda nstart = 1, a nstop =2
//t.e. fakticheski etogo cikla dlia simuliacii net, a prosto dlia kazdogo sobitia vse schitaetsia
//v to vremia kak dlia dannih cikl po j srabativaet ne kazdi' raz, a tol'ko kogda nstart ne ravno nstop i tam mi uze vozvraschaemsia na skol'ko-to sobiti' nazad

  if (data_sim==1){
   //dlya dannih nuzno eshe raz vzyat' cobitiya t.k. mi vernulis' nazad - k nachalu blocka. esly ne vzyat' zanovo, to u peremennoj budet to znachenie, chto bilo u nee v kontse blocka
  br_npart-> GetEntry(j);
  br_gpart-> GetEntry(j);
  br_stat->GetEntry(j); 
  br_p-> GetEntry(j);
  br_b-> GetEntry(j);
  br_q-> GetEntry(j);
  br_dc-> GetEntry(j);
  br_sc-> GetEntry(j);
  br_ec-> GetEntry(j); 
  br_cc-> GetEntry(j);
  br_cx-> GetEntry(j);
  br_cy-> GetEntry(j); 
  br_cz-> GetEntry(j);  
  br_dc_stat-> GetEntry(j);  
  br_q_l-> GetEntry(j);
 // br_tr_time-> GetEntry(j);
  br_t_l-> GetEntry(j); 
  br_nphe-> GetEntry(j); 
  br_etot-> GetEntry(j);  
  br_ec_ei-> GetEntry(j); 
  br_ec_eo-> GetEntry(j);   
  br_vx-> GetEntry(j);
  br_vy-> GetEntry(j);
  br_vz-> GetEntry(j);      
  br_sc_sect-> GetEntry(j);
  br_sc_part-> GetEntry(j);  
  br_sc_pd-> GetEntry(j); 
  br_dc_xsc-> GetEntry(j);
  br_dc_ysc-> GetEntry(j);   
  br_dc_zsc-> GetEntry(j); 
  br_cc_segm-> GetEntry(j);     
  br_dc_cxsc-> GetEntry(j);
  br_dc_cysc-> GetEntry(j);   
  br_dc_czsc-> GetEntry(j);  
  br_sc_t-> GetEntry(j);
  br_sc_r-> GetEntry(j);
  
  br_dc_vx -> GetEntry(j);
  br_dc_vy -> GetEntry(j);
  br_dc_vz -> GetEntry(j);
 }; 
 
 //--------------------------
 if (data_sim==2){
 indtype = 2;
 
/* th_EL=br_mctheta->GetLeaf("mctheta")->GetValue(0);
 th_PIm=br_mctheta->GetLeaf("mctheta")->GetValue(1);
 th_PIp=br_mctheta->GetLeaf("mctheta")->GetValue(2);
 th_P=br_mctheta->GetLeaf("mctheta")->GetValue(3);
 
 ph_EL=br_mcphi->GetLeaf("mcphi")->GetValue(0);
 ph_PIm=br_mcphi->GetLeaf("mcphi")->GetValue(1);
 ph_PIp=br_mcphi->GetLeaf("mcphi")->GetValue(2);
 ph_P=br_mcphi->GetLeaf("mcphi")->GetValue(3);
 
 P_EL=br_mcp->GetLeaf("mcp")->GetValue(0);
 P_PIm=br_mcp->GetLeaf("mcp")->GetValue(1);
 P_PIp=br_mcp->GetLeaf("mcp")->GetValue(2);
 P_P=br_mcp->GetLeaf("mcp")->GetValue(3);*/
 for (l=0; l<br_mcnpart->GetLeaf("mcnpart")->GetValue(); l++){
 
 if (br_mcid->GetLeaf("mcid")->GetValue(l) == 11){
  th_EL=br_mctheta->GetLeaf("mctheta")->GetValue(l);
  ph_EL=br_mcphi->GetLeaf("mcphi")->GetValue(l);
  P_EL=br_mcp->GetLeaf("mcp")->GetValue(l);
 };
 
 if (br_mcid->GetLeaf("mcid")->GetValue(l) == -211){ 
 th_PIm=br_mctheta->GetLeaf("mctheta")->GetValue(l); 
 ph_PIm=br_mcphi->GetLeaf("mcphi")->GetValue(l); 
 P_PIm=br_mcp->GetLeaf("mcp")->GetValue(l); 
 };

 if (br_mcid->GetLeaf("mcid")->GetValue(l) == 211){
  th_PIp=br_mctheta->GetLeaf("mctheta")->GetValue(l);
  ph_PIp=br_mcphi->GetLeaf("mcphi")->GetValue(l);
  P_PIp=br_mcp->GetLeaf("mcp")->GetValue(l);
 }; 

 if (br_mcid->GetLeaf("mcid")->GetValue(l) == 2212){
 th_P=br_mctheta->GetLeaf("mctheta")->GetValue(l);
 ph_P=br_mcphi->GetLeaf("mcphi")->GetValue(l);
 P_P=br_mcp->GetLeaf("mcp")->GetValue(l);
};

}; 
 x_EL = br_mcvx_x_el ->GetLeaf("mcvx_x_el")->GetValue();
 y_EL = br_mcvx_y_el ->GetLeaf("mcvx_y_el")->GetValue();
 z_EL = br_mcvx_z_el ->GetLeaf("mcvx_z_el")->GetValue();
 
 sigma = br_sigma_total ->GetLeaf("sigma_total")->GetValue();
 px_fermi = br_px_Fermi ->GetLeaf("px_Fermi")->GetValue();
 py_fermi = br_py_Fermi ->GetLeaf("py_Fermi")->GetValue();
 pz_fermi = br_pz_Fermi ->GetLeaf("pz_Fermi")->GetValue();
 
// pf_x = br_pf_x ->GetLeaf("pf_x")->GetValue();
//  pf_y = br_pf_y ->GetLeaf("pf_y")->GetValue();
//  pf_z = br_pf_z ->GetLeaf("pf_z")->GetValue();
 
 if ((ph_EL < 30.) || (ph_EL > 330.)) {
sector = 1;
};
if ((ph_EL > 30.) && (ph_EL < 90.)) {
sector = 2;
};
if ((ph_EL > 90.) && (ph_EL < 150.)) {
sector = 3;
};
if ((ph_EL > 150.) && (ph_EL < 210.)) {
sector = 4;
};
if ((ph_EL > 210.) && (ph_EL < 270.)) {
sector = 5;
};
if ((ph_EL > 270.) && (ph_EL < 330.)) {
sector = 6;
};
 
 npart = 4;
 block = -1000;
 deltaQ = -1000.;
 LiveTime = -1000.;
 n_incl = -1000;
 n_elast = -1000;
 PdHit_PIp = -1000;
 PdHit_PIm = -1000;
 PdHit_P = -1000; 
 PdHit_EL = -1000.; 
 pdhit = -1000.;
 beta_PIm = -1000.;
 beta_PIp = -1000.;
 beta_P = -1000.;  
 NpheCC_EL = -1000;
 NpheCC_PIp = -1000;
 ECtot_EL =-1000;
 ECin_EL =-1000;
 ECout_EL = -1000;
 
  z_P = -1000.;
 z_PIm = -1000.;
 z_PIp = -1000.;
 sc_x= -1000.; 
 sc_y= -1000.;
 sc_z= -1000.;
 pmt_hit = -1000;
 segment = -1000;  
 theta_cc = -1000;
 ph_cc = -1000;
// sector = -1000;
 n_PIp = -1000;
 n_PIm = -1000;
 n_P = -1000;
 PIp_dist = -1000;
 PIp_time = -1000;
 PIm_dist = -1000;
 PIm_time = -1000;
 P_dist = -1000;
 P_time = -1000;
 dc_x_EL=-1000.;
 dc_y_EL=-1000.;
 dc_z_EL=-1000.;
 
  dc_x_P=-1000.;
 dc_y_P=-1000.;
 dc_z_P=-1000.;
 
 dc_x_PIp=-1000.;
 dc_y_PIp=-1000.;
 dc_z_PIp=-1000.;
 
  dc_x_PIm=-1000.;
 dc_y_PIm=-1000.;
 dc_z_PIm=-1000.;
   
     Q2 = 4.*E0*P_EL*(sin(th_EL*M_PI/2./180.))*(sin(th_EL*M_PI/2./180.));
     W = m_proton*m_proton+2.*m_proton*(E0-P_EL)-Q2;
     W = sqrt(W);
        
      
 t21->Fill();//zapolnyaem derevo dlya indtype==2 (dlya sgenerirovannih). 
 indtype = 1;
 
  dc_x_EL=-1000.;
 dc_y_EL=-1000.;
 dc_z_EL=-1000.;
 
  
  dc_x_P=-1000.;
 dc_y_P=-1000.;
 dc_z_P=-1000.;
 
 dc_x_PIp=-1000.;
 dc_y_PIp=-1000.;
 dc_z_PIp=-1000.;
 
  dc_x_PIm=-1000.;
 dc_y_PIm=-1000.;
 dc_z_PIm=-1000.;
 
 
  };
 //--------------------------  
           
   npart = br_npart->GetLeaf("npart")->GetValue();
   gpart = br_gpart->GetLeaf("gpart")->GetValue();
  if (npart> gpart) cout << npart << "  "<<gpart <<"\n";
//   cout << br_dc->GetLeaf("dc")->GetValue(0) << "\n";

   P_EL = br_p->GetLeaf("p")->GetValue(0);
 
     n_PIp = 0;
     n_PIm = 0;
     n_P = 0;
     
    // for (k=0; k<10; k++) {
     P_PIp = -1000.;
     P_PIm = -1000.;
     P_P = -1000.;  
     th_PIp = -1000.;
     th_PIm = -1000.;
     th_P = -1000.; 
     ph_PIp = -1000.;
     ph_PIm = -1000.;
     ph_P = -1000.; 
     PdHit_PIp = 0;
     PdHit_PIm = 0;
     PdHit_P = 0; 
     PdHit_EL = 0;
     pdhit = 0;
     beta_PIm = -1000.;
     beta_PIp = -1000.;
     beta_PIp_time = -1000.;
     beta_PIm_time = -1000.;
     beta_P_time = -1000.;
     beta_P = -1000.;   
     NpheCC_PIp = -1000; 
     z_P = -1000.;
     z_PIm = -1000.;
     z_PIp = -1000.;    
     PIp_dist = -1000;
     PIp_time = -1000;
     PIm_dist = -1000;
     PIm_time = -1000;
     P_dist = -1000;
     P_time = -1000; 
     
  //   };
    
 // cout << "grart = " << gpart << "\n";  
    for (k=1; k < gpart; k++) {
    
    if (br_gpart->GetLeaf("gpart")->GetValue() > 0){
    if (br_stat->GetLeaf("stat")->GetValue(k) > 0){
    if (br_dc->GetLeaf("dc")->GetValue(k) > 0) {
     if (br_sc->GetLeaf("sc")->GetValue(k) > 0) {
     
     if (br_q->GetLeaf("q")->GetValue(k) == 1) {
     beta = br_b->GetLeaf("b")->GetValue(k);
     hadr_mom = br_p->GetLeaf("p")->GetValue(k); 
     
      if (hadr_mom <= E0) {
     
    // pion_low = (1+5*1.4*(hadr_mom-0.07))/(1+5*(hadr_mom-0.07));
    // pion_low = pion_low*(hadr_mom-0.07)/sqrt((hadr_mom-0.07)*(hadr_mom-0.07)+0.138*0.138);
     //pion_low = pion_low - 0.4;
     pion_low =hadr_mom/sqrt(hadr_mom*hadr_mom+0.938*0.938) + 0.03;
     pion_low = pion_low*(1.2+0.92*hadr_mom)/(1+hadr_mom);
     pion_high = 2.;
     
     
   //---------------******* PI+ ID*****-------------------------------------
   
   
     if ((beta > pion_low) && (beta < pion_high)) {
     
     n_PIp = n_PIp + 1;
     if (n_PIp == 1 ) {
     P_PIp = hadr_mom;
     z_PIp = br_vz->GetLeaf("vz")->GetValue(k);
     beta_PIp = beta;
     
     PIp_dist=br_sc_r->GetLeaf("sc_r")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     EL_dist = br_sc_r->GetLeaf("sc_r")->GetValue(br_sc->GetLeaf("sc")->GetValue(0)-1);
     EL_tof = br_sc_t->GetLeaf("sc_t")->GetValue(br_sc->GetLeaf("sc")->GetValue(0)-1);
     //EL_dist = br_sc_r->GetLeaf("sc_r")->GetValue(0);//el_dist/c - vremya vileta iz misheni
     //sEL_tof = br_sc_t->GetLeaf("sc_t")->GetValue(0);
     PIp_tof = br_sc_t->GetLeaf("sc_t")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     PIp_time = EL_dist/c -EL_tof + PIp_tof;
     //PIp_time = PIp_tof-tr_time;
   
     beta_PIp_time = PIp_dist/PIp_time/c;
     
     sc_pd_local = br_sc_pd->GetLeaf("sc_pd")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     sc_sect_local = br_sc_sect->GetLeaf("sc_sect")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     sc_part_local = br_sc_part->GetLeaf("sc_part")->GetValue();
     PdHit_PIp = 10000.*sc_sect_local + 100.*sc_pd_local;
     PdHit_PIp = (PdHit_PIp-10000.*int(PdHit_PIp/10000.))/100.;
     
     NpheCC_PIp = br_nphe->GetLeaf("nphe")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     
     th_PIp = (180./M_PI)*acos(br_cz->GetLeaf("cz")->GetValue(k));
     //sc_x_pip = br_dc_xsc->GetLeaf("dc_xsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(k)-1);

     //sc_y_pip = br_dc_ysc->GetLeaf("dc_ysc")->GetValue(br_dc->GetLeaf("dc")->GetValue(k)-1);

    // sc_z_pip = br_dc_zsc->GetLeaf("dc_zsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(k)-1);
     
     				  if (br_cx->GetLeaf("cx")->GetValue(k) != 0.) {
			         ph_PIp = (180./M_PI)*atan(br_cy->GetLeaf("cy")->GetValue(k)/br_cx->GetLeaf("cx")->GetValue(k));
			         } else {
		           	  if(br_cy->GetLeaf("cy")->GetValue(k) > 0.) ph_PIp = 90.;
			           if(br_cy->GetLeaf("cy")->GetValue(k) < 0.) ph_PIp = 270.;
		        	   };
				   

if ((br_cx->GetLeaf("cx")->GetValue(k) < 0.) && (br_cy->GetLeaf("cy")->GetValue(k) > 0)) ph_PIp = ph_PIp+180.;
if ((br_cx->GetLeaf("cx")->GetValue(k) < 0.) && (br_cy->GetLeaf("cy")->GetValue(k) < 0)) ph_PIp = ph_PIp+180.;
if ((br_cx->GetLeaf("cx")->GetValue(k) > 0.) && (br_cy->GetLeaf("cy")->GetValue(k) < 0)) ph_PIp = ph_PIp+360.;	

if (data_sim ==1){
dc_x_PIp = br_dc_vx->GetLeaf("dc_vx")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1); 
dc_y_PIp = br_dc_vy->GetLeaf("dc_vy")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
dc_z_PIp = br_dc_vz->GetLeaf("dc_vz")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);  

};

     }; //end if n_PIp=1            
     
     };//end if PIp low_high cut
     
     
 //---------------******* PROTON ID*****-------------------------------------
     
     //proton_low = hadr_mom/sqrt(hadr_mom*hadr_mom+0.938*0.938) - 0.075;
     //proton_low = proton_low*(0.9+1.06*hadr_mom)/(1+hadr_mom);
     proton_low = 0.;
    // proton_high =hadr_mom/sqrt(hadr_mom*hadr_mom+0.938*0.938) + 0.05;//+0.03
    // proton_high = proton_high*(1.2+0.92*hadr_mom)/(1+hadr_mom);
     proton_high  = (1+5*1.4*(hadr_mom-0.07))/(1+5*(hadr_mom-0.07));
     proton_high = proton_high *(hadr_mom-0.07)/sqrt((hadr_mom-0.07)*(hadr_mom-0.07)+0.138*0.138);
     proton_high = proton_high - 0.4;
     
       if ((beta > proton_low) && (beta < proton_high)) {
     n_P = n_P + 1;
     if (n_P == 1) {
    P_P = hadr_mom;
    z_P = br_vz->GetLeaf("vz")->GetValue(k); 
    beta_P = beta;
    
     //EL_dist = br_sc_r->GetLeaf("sc_r")->GetValue(0);//el_dist/c - vremya vileta iz misheni
    // EL_tof = br_sc_t->GetLeaf("sc_t")->GetValue(0);
     EL_dist = br_sc_r->GetLeaf("sc_r")->GetValue(br_sc->GetLeaf("sc")->GetValue(0)-1);
     EL_tof = br_sc_t->GetLeaf("sc_t")->GetValue(br_sc->GetLeaf("sc")->GetValue(0)-1);
     P_dist=br_sc_r->GetLeaf("sc_r")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     P_tof = br_sc_t->GetLeaf("sc_t")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
    P_time = EL_dist/c -EL_tof + P_tof;
    //P_time = P_tof-tr_time;
    
    beta_P_time = P_dist/P_time/c;
     
    sc_pd_local = br_sc_pd->GetLeaf("sc_pd")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     sc_sect_local = br_sc_sect->GetLeaf("sc_sect")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     sc_part_local = br_sc_part->GetLeaf("sc_part")->GetValue();
     PdHit_P = 10000.*sc_sect_local + 100.*sc_pd_local;
PdHit_P = (PdHit_P-10000.*int(PdHit_P/10000.))/100.;
    
    th_P = (180./M_PI)*acos(br_cz->GetLeaf("cz")->GetValue(k));
//   if (data_sim == 1) 
    P_P = P_P + corrfunc.correct_energy_theta_pf(P_P, th_P);
   // sc_x_p = br_dc_xsc->GetLeaf("dc_xsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(k)-1);

    // sc_y_p = br_dc_ysc->GetLeaf("dc_ysc")->GetValue(br_dc->GetLeaf("dc")->GetValue(k)-1);

    // sc_z_p = br_dc_zsc->GetLeaf("dc_zsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(k)-1);
    if (br_cx->GetLeaf("cx")->GetValue(k) != 0.) {
			         ph_P = (180./M_PI)*atan(br_cy->GetLeaf("cy")->GetValue(k)/br_cx->GetLeaf("cx")->GetValue(k));
			         } else {
		           	  if(br_cy->GetLeaf("cy")->GetValue(k) > 0.) ph_P = 90.;
			           if(br_cy->GetLeaf("cy")->GetValue(k) < 0.) ph_P = 270.;
		        	   };
				   

if ((br_cx->GetLeaf("cx")->GetValue(k) < 0.) && (br_cy->GetLeaf("cy")->GetValue(k) > 0)) ph_P = ph_P+180.;
if ((br_cx->GetLeaf("cx")->GetValue(k) < 0.) && (br_cy->GetLeaf("cy")->GetValue(k) < 0)) ph_P = ph_P+180.;
if ((br_cx->GetLeaf("cx")->GetValue(k) > 0.) && (br_cy->GetLeaf("cy")->GetValue(k) < 0)) ph_P = ph_P+360.;


if (data_sim ==1){
dc_x_P = br_dc_vx->GetLeaf("dc_vx")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1); 
dc_y_P = br_dc_vy->GetLeaf("dc_vy")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
dc_z_P = br_dc_vz->GetLeaf("dc_vz")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);  

};


     }; //end if n_P=1   
     
     };//end if P low_high cut
     
     };//end if mom < E0
     
     };//end if q=1



//cout << indtype << "\n";


 //---------------******* PI-  ID*****-------------------------------------
 
     if (br_q->GetLeaf("q")->GetValue(k) == -1) { 
     beta = br_b->GetLeaf("b")->GetValue(k);
     hadr_mom = br_p->GetLeaf("p")->GetValue(k); 

      if (hadr_mom <= E0) {
     
     //pion_low = (1+5*1.4*(hadr_mom-0.07))/(1+5*(hadr_mom-0.07));
     //pion_low = pion_low*(hadr_mom-0.07)/sqrt((hadr_mom-0.07)*(hadr_mom-0.07)+0.138*0.138);
     //pion_low = pion_low - 0.4;
     //pion_high = 1.4;
     pion_low = hadr_mom/sqrt(hadr_mom*hadr_mom+0.938*0.938) + 0.03;
     pion_low = pion_low*(1.2+0.92*hadr_mom)/(1+hadr_mom);
     pion_high = 2.;
     
     if ((beta > pion_low) && (beta < pion_high)) {
     
     n_PIm = n_PIm + 1;
     if (n_PIm == 1) {
     P_PIm = hadr_mom;
     z_PIm = br_vz->GetLeaf("vz")->GetValue(k);
     beta_PIm = beta; 
     
     //EL_dist = br_sc_r->GetLeaf("sc_r")->GetValue(0);//el_dist/c - vremya vileta iz misheni
     //EL_tof = br_sc_t->GetLeaf("sc_t")->GetValue(0);
     EL_dist = br_sc_r->GetLeaf("sc_r")->GetValue(br_sc->GetLeaf("sc")->GetValue(0)-1);
     EL_tof = br_sc_t->GetLeaf("sc_t")->GetValue(br_sc->GetLeaf("sc")->GetValue(0)-1);
     PIm_dist = br_sc_r->GetLeaf("sc_r")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     PIm_tof = br_sc_t->GetLeaf("sc_t")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     PIm_time = EL_dist/c -EL_tof + PIm_tof;
    //PIm_time = PIm_tof-tr_time;
     
     beta_PIm_time = PIm_dist/PIm_time/c;
     
     sc_pd_local = br_sc_pd->GetLeaf("sc_pd")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     sc_sect_local = br_sc_sect->GetLeaf("sc_sect")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
     sc_part_local = br_sc_part->GetLeaf("sc_part")->GetValue();
     PdHit_PIm = 10000.*sc_sect_local + 100.*sc_pd_local;
PdHit_PIm = (PdHit_PIm-10000.*int(PdHit_PIm/10000.))/100.;
     
     th_PIm = (180./M_PI)*acos(br_cz->GetLeaf("cz")->GetValue(k));
    // sc_x_pim = br_dc_xsc->GetLeaf("dc_xsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(k)-1);

    // sc_y_pim = br_dc_ysc->GetLeaf("dc_ysc")->GetValue(br_dc->GetLeaf("dc")->GetValue(k)-1);

    // sc_z_pim = br_dc_zsc->GetLeaf("dc_zsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(k)-1);
     
     				  if (br_cx->GetLeaf("cx")->GetValue(k) != 0.) {
			         ph_PIm = (180./M_PI)*atan(br_cy->GetLeaf("cy")->GetValue(k)/br_cx->GetLeaf("cx")->GetValue(k));
			         } else {
		           	  if(br_cy->GetLeaf("cy")->GetValue(k) > 0.) ph_PIm = 90.;
			           if(br_cy->GetLeaf("cy")->GetValue(k) < 0.) ph_PIm = 270.;
		        	   };
				   

if ((br_cx->GetLeaf("cx")->GetValue(k) < 0.) && (br_cy->GetLeaf("cy")->GetValue(k) > 0)) ph_PIm = ph_PIm+180.;
if ((br_cx->GetLeaf("cx")->GetValue(k) < 0.) && (br_cy->GetLeaf("cy")->GetValue(k) < 0)) ph_PIm = ph_PIm+180.;
if ((br_cx->GetLeaf("cx")->GetValue(k) > 0.) && (br_cy->GetLeaf("cy")->GetValue(k) < 0)) ph_PIm = ph_PIm+360.;	



if (data_sim ==1){
dc_x_PIm = br_dc_vx->GetLeaf("dc_vx")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1); 
dc_y_PIm = br_dc_vy->GetLeaf("dc_vy")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);
dc_z_PIm = br_dc_vz->GetLeaf("dc_vz")->GetValue(br_sc->GetLeaf("sc")->GetValue(k)-1);  

}; 


     };//end if n_PIm=1  
     };// end if PIm low high cut       
     
     };//end if mom < E0
     
     };//end if q=-1    






     
     };//end if sc signal
     };//end if dc signal
     };//end if stat
     };//end if gpart
     };//end cycle over k?? over particles
     
     
     
     
 //---------------******* ELECTRON ID*****-------------------------------------     
 
 
    if (br_gpart->GetLeaf("gpart")->GetValue() > 0){
    if (br_stat->GetLeaf("stat")->GetValue(0) > 0){
    if (br_dc->GetLeaf("dc")->GetValue(0) > 0) {
    if (br_sc->GetLeaf("sc")->GetValue(0) > 0) {
    if (br_ec->GetLeaf("ec")->GetValue(0) > 0) {
    if (br_cc->GetLeaf("cc")->GetValue(0) > 0) {
    if (br_q->GetLeaf("q")->GetValue(0) == -1) {    
    
    EL_dist = br_sc_r->GetLeaf("sc_r")->GetValue(br_sc->GetLeaf("sc")->GetValue(0)-1);
   // EL_time_fr_vert=EL_dist/c;
     EL_tof = br_sc_t->GetLeaf("sc_t")->GetValue(br_sc->GetLeaf("sc")->GetValue(0)-1);
    pmt_hit = int((br_cc_segm->GetLeaf("cc_segm")->GetValue(br_cc->GetLeaf("cc")->GetValue(0)-1))/1000)-1;
    segment = ((int(br_cc_segm->GetLeaf("cc_segm")->GetValue(br_cc->GetLeaf("cc")->GetValue(0)-1)))%1000/10);
    if ((segment<0)||(segment>18)) cout << "segment="<< segment <<"\n";
     
 //   cout << "cc_segm = " <<  int((br_cc_segm->GetLeaf("cc_segm")->GetValue(0))/1000) << "\n";
    
    
   
    
    NpheCC_EL = br_nphe->GetLeaf("nphe")->GetValue(br_cc->GetLeaf("cc")->GetValue(0)-1);
    ECtot_EL = br_etot->GetLeaf("etot")->GetValue(br_ec->GetLeaf("ec")->GetValue(0)-1);
    
    ECin_EL = br_ec_ei->GetLeaf("ec_ei")->GetValue(br_ec->GetLeaf("ec")->GetValue(0)-1); 
    
    ECout_EL = br_ec_eo->GetLeaf("ec_eo")->GetValue(br_ec->GetLeaf("ec")->GetValue(0)-1);        
    
				  if (br_cx->GetLeaf("cx")->GetValue(0) != 0.) {
			         ph_EL = (180./M_PI)*atan(br_cy->GetLeaf("cy")->GetValue(0)/br_cx->GetLeaf("cx")->GetValue(0));
			         } else {
		           	  if(br_cy->GetLeaf("cy")->GetValue(0) > 0.) ph_EL = 90.;
			           if(br_cy->GetLeaf("cy")->GetValue(0) < 0.) ph_EL = 270.;
		        	   };
				   

if ((br_cx->GetLeaf("cx")->GetValue(0) < 0.) && (br_cy->GetLeaf("cy")->GetValue(0) > 0)) ph_EL = ph_EL+180.;
if ((br_cx->GetLeaf("cx")->GetValue(0) < 0.) && (br_cy->GetLeaf("cy")->GetValue(0) < 0)) ph_EL = ph_EL+180.;
if ((br_cx->GetLeaf("cx")->GetValue(0) > 0.) && (br_cy->GetLeaf("cy")->GetValue(0) < 0)) ph_EL = ph_EL+360.;				   

th_EL = (180./M_PI)*acos(br_cz->GetLeaf("cz")->GetValue(0));

if (data_sim == 1){
P_EL_new = corrfunc.correct_pel_e1_2039_2250_feb09(P_EL,th_EL,ph_EL);
th_EL_new = corrfunc.correct_thel_e1_2039_2250_feb09(P_EL,th_EL,ph_EL);
P_EL = P_EL_new;
th_EL = th_EL_new;
};	
			   
Q2 = 4.*E0*P_EL*(sin(th_EL*M_PI/2./180.))*(sin(th_EL*M_PI/2./180.));
W = m_proton*m_proton+2.*m_proton*(E0-P_EL)-Q2;


sc_pd_local = br_sc_pd->GetLeaf("sc_pd")->GetValue(br_sc->GetLeaf("sc")->GetValue(0)-1);

sc_sect_local = br_sc_sect->GetLeaf("sc_sect")->GetValue(br_sc->GetLeaf("sc")->GetValue(0)-1);

sc_part_local = br_sc_part->GetLeaf("sc_part")->GetValue();



PdHit_EL = 10000.*sc_sect_local + 100.*sc_pd_local;

//cout << "sc_sect_local = " << sc_sect_local << " sc_pd_local = " << sc_pd_local << "\n";


sc_x = br_dc_xsc->GetLeaf("dc_xsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1);

sc_y = br_dc_ysc->GetLeaf("dc_ysc")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1);

sc_z = br_dc_zsc->GetLeaf("dc_zsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1);


if (abs(br_dc_cxsc->GetLeaf("dc_cxsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1)) > 0) {
nx=(br_dc_cxsc->GetLeaf("dc_cxsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1));
} else {
nx = 0;
};

if (abs(br_dc_cysc->GetLeaf("dc_cysc")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1)) > 0) {
ny=(br_dc_cysc->GetLeaf("dc_cysc")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1));
} else {
ny = 0;
};

if (abs(br_dc_czsc->GetLeaf("dc_czsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1)) > 0) {
nz=(br_dc_czsc->GetLeaf("dc_czsc")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1));
} else {
nz = 0;
};

t=abs((sx*sc_x+sy*sc_y+sz*sc_z+1.)/(sx*nx+sy*ny+sz*nz));



px=sc_x+t*nx;
py=sc_y+t*ny;
pz=sc_z+t*nz;

//cout << pz<<" "<< -sx/sz*px -1./sz<<" \n";

if ((ph_EL < 30.) || (ph_EL > 330.)) {
sector = 1;
};
if ((ph_EL > 30.) && (ph_EL < 90.)) {
sector = 2;
};
if ((ph_EL > 90.) && (ph_EL < 150.)) {
sector = 3;
};
if ((ph_EL > 150.) && (ph_EL < 210.)) {
sector = 4;
};
if ((ph_EL > 210.) && (ph_EL < 270.)) {
sector = 5;
};
if ((ph_EL > 270.) && (ph_EL < 330.)) {
sector = 6;
};

if (abs(pz/(sqrt(px*px+py*py+pz*pz))) <= 1.) {
	       theta_cc=acos(pz/(sqrt(px*px+py*py+pz*pz)))*180./M_PI;
} else {
	       theta_cc = 0.;
};


ph_cc = atan2(py,px)*180./M_PI;


PdHit_EL = (PdHit_EL-10000.*int(PdHit_EL/10000.))/100.;
pdhit = PdHit_EL;

				   
x_EL = br_vx->GetLeaf("vx")->GetValue(0); 
y_EL = br_vy->GetLeaf("vy")->GetValue(0);
z_EL = br_vz->GetLeaf("vz")->GetValue(0);  
  
if (data_sim ==1){
dc_x_EL = br_dc_vx->GetLeaf("dc_vx")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1); 
dc_y_EL = br_dc_vy->GetLeaf("dc_vy")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1);
dc_z_EL = br_dc_vz->GetLeaf("dc_vz")->GetValue(br_dc->GetLeaf("dc")->GetValue(0)-1);  

};    
//    if (br_dc_stat->GetLeaf("dc_stat")->GetValue(0) > 0 ) {
    if (br_q->GetLeaf("q")->GetValue(0) == -1) { 
   if ((P_EL > 0.) && (P_EL < E0) && (W > 0.)) { 
    k=k+1;
    W = sqrt(W);
    
    if (data_sim==1){
    block = block - 1;
    block_total = block + block_last;
//    cout << "block numer = " << block_total << "\n";
    };
    segment = segment - 1;
    if (data_sim==1){
     br_t_l-> GetEntry(nstop); 
//    cout << "l_time = " <<  br_t_l->GetLeaf("t_l")->GetValue(0)  << " block = " <<  block << "\n";
    LiveTime = br_t_l->GetLeaf("t_l")->GetValue(0);
    };
    t21->Fill();//zapolnyaem derevo dlya vosstanovlennih. toliko esli est' electron
    segment = segment + 1;
    
    if (data_sim==1){
    block = block + 1;
    };
    
    };
    };
    
    
    };//end if q=-1
    };//end if cc signal
    };//end if ec signal
    
    };//end if sc signal
    };//end if dc signal
    };//end if stat
    };//end if gpart
    
    };//konets tsikla po j (po sobitiyam vnutri blocka)
    
    if (data_sim==1){
    if (i == last_i) {
    last_k = k;
    };    
    };

    }; //end if block >1
    
  if (data_sim==1){
  if (block > 1) {
  Qtotal = Qtotal + Qdiff;
  last_i = i;
  deltaQ = Qdiff;
  n_incl = 0;
  n_elast = 0;
  };
  };
   }; //end if Qdiff >0 
 
  
   };//end tsikl po i - po sobitiyam
 
   
  Qfull = Qfull + Qtotal-Qdiff;
  t20->Delete();
  finp->Close();
  block_last = block_last + block - 1;
 

 
 };//konets tsikla po failam (m)
 
 
 
     cout << "Q full = " << Qfull << "\n";
     TFile *outFile;
     outFile = new TFile(outfile_inp.c_str(),"recreate");
     outFile->cd();
     t21->Write("", TObject::kOverwrite);
     outFile->Write();
     outFile->Close();
     t21->Delete();
     
        
    
    };//end void

    
     
