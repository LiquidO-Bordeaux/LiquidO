#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include "TROOT.h"
#include "TStyle.h"

#include <stdio.h>

TFile *current_file = NULL;
TTree *event_tree;

float waveform[64][1024];
bool ch_enable[64];


void wavecatcher_analysis()
{
  // gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  gStyle->SetPalette(55);
}

void open_file (const char *filename)
{
  if (current_file != NULL) current_file->Close();
  current_file = new TFile (filename, "READ");
  event_tree = (TTree*)(current_file->Get("event_tree"));

  for (int ch=0; ch<64; ++ch) {
    if (event_tree->FindBranch(Form("ch%02d_waveform",ch)) != NULL) {
      event_tree->SetBranchAddress(Form("ch%02d_waveform",ch), waveform[ch]);
      printf("ch %02d detected\n", ch);
      ch_enable[ch] = true;}
    else ch_enable[ch] = false;}
}


void analyse_baseline()
{
  TH1F *baseline_histo = new TH1F ("baseline_histo", "", 18, 0, 18);

  for (int ch=0; ch<18; ++ch) {
    
    if (!ch_enable[ch]) continue;

    TH1F tmp_histo ("tmp_histo", "", 1000, -100, 100);
    event_tree->Project("tmp_histo", Form("ch%02d_waveform",ch));

    TF1 tmp_fit ("tmp_fit", "gaus", -100, 100);
    tmp_fit.SetParameters(tmp_histo.GetMaximum(), tmp_histo.GetMean(), tmp_histo.GetRMS());
    tmp_histo.Fit("tmp_fit", "0QR");

    baseline_histo->SetBinContent(1+ch, tmp_fit.GetParameter(1));
    baseline_histo->SetBinError(1+ch, tmp_fit.GetParameter(2));}

  baseline_histo->Draw("e");

}

TCanvas *Cevent = NULL;
TH1F *waveform_histo[64];

void draw_event (int event_id, int rebin=0)
{
  if (current_file == NULL) return;
  
  if (Cevent == NULL) {
    Cevent = new TCanvas ("Cevent", "", 1250, 500);
    Cevent->SetLeftMargin(0.05);}
  
  for (int ch=0; ch<64; ++ch) {

    if (!ch_enable[ch]) continue;

    waveform_histo[ch] = new TH1F (Form("ch%02d_waveform_histo",ch), "", 1024, 0, 0.3125*1024);
    waveform_histo[ch]->SetLineColor(1+ch);}
  
  event_tree->GetEntry(event_id);
  
  Cevent->cd();
  
  for (int ch=0; ch<64; ++ch) {

    if (!ch_enable[ch]) continue;
    
    for (int s=0; s<1024; ++s) {
      waveform_histo[ch]->SetBinContent(1+s, waveform[ch][s]);}
    waveform_histo[ch]->Draw(ch==0 ? "" : "same");}


  if (rebin == 0) return;
  
  for (int ch=0; ch<64; ++ch) {

    if (!ch_enable[ch]) continue;

    waveform_histo[ch]->Rebin(rebin);
    
    for (int s=0; s<1024/rebin; ++s)
      waveform_histo[ch]->SetBinContent(1+s, waveform_histo[ch]->GetBinContent(1+s));
  }
  
    
}


TCanvas *Cpersistency = NULL;
TH2F *persistency_waveform = NULL;

void channel_persistency (int ch, int biny, float ymin, float ymax)
{
  if (current_file == NULL) return;
  
  if (persistency_waveform != NULL) delete persistency_waveform;
  persistency_waveform = new TH2F ("persistency_waveform", "", 1024, 0, 1024*0.3125, biny, ymin, ymax);
  
  for (int entry=0; entry<event_tree->GetEntries(); ++entry) {
    event_tree->GetEntry(entry); for (int s=0; s<1024; ++s)
      persistency_waveform->Fill(0.3125*(s+0.5), waveform[ch][s]);}

  if (Cpersistency == NULL) {
    Cpersistency = new TCanvas ("Cpersistency", "", 1250, 500);
    Cpersistency->SetLeftMargin(0.05);
    Cpersistency->SetLogz();}

  Cpersistency->cd();
  persistency_waveform->Draw("colz");
}


const float baseline_range = 50; // ns
// int baseline_sample = 80; // 25 ns

const float charge_start =   75; // ns
const float charge_stop  =  300; // ns

void reconstruction(const char *filename)
{
  if (current_file != NULL) current_file->Close();
  current_file = new TFile (filename, "UPDATE");
  event_tree = (TTree*)(current_file->Get("event_tree"));

  bool ch_enable[64];
  float waveforms[64][1024];

  float baselines[64];
  float amplitudes_min[64];
  float amplitudes_max[64];
  float charges[64];
  float tstarts[64];

  TTree *reco_tree = new TTree ("reco_tree", "");

  ////////////////////////////////

  for (int ch=0; ch<64; ++ch) {

    if (event_tree->FindBranch(Form("ch%02d_waveform",ch)) == NULL) {
      ch_enable[ch] = false; continue;}

    printf("ch %02d detected\n", ch);
    ch_enable[ch] = true;

    event_tree->SetBranchAddress(Form("ch%02d_waveform",ch), waveforms[ch]);

    reco_tree->Branch(Form("ch%02d_baseline",ch), &baselines[ch]);
    reco_tree->Branch(Form("ch%02d_amplitude_min",ch), &amplitudes_min[ch]);
    reco_tree->Branch(Form("ch%02d_amplitude_max",ch), &amplitudes_max[ch]);
    reco_tree->Branch(Form("ch%02d_charge",ch), &charges[ch]);
    reco_tree->Branch(Form("ch%02d_tstart",ch), &tstarts[ch]);
    reco_tree->Branch(Form("ch%02d_tmat",ch), &tstarts[ch]);
  }

  const int baseline_sample = baseline_range/0.3125;
  printf("+++ baseline : sample 0 -> %d\n", baseline_sample);
  
  const int charge_start_sample = charge_start/0.3125;
  const int charge_stop_sample = charge_stop/0.3125;
  printf("+++ charge : sample %d -> %d\n", charge_start_sample, charge_stop_sample);
      
      ////////////////////////////////
      
  const int event_entries = event_tree->GetEntries();

  for (int entry=0; entry<event_entries; ++entry) {

    event_tree->GetEntry(entry);

    for (int ch=0; ch<64; ++ch) {

      if (!ch_enable[ch]) continue;

      baselines[ch] = 0;
      amplitudes_min[ch] = 1250;
      amplitudes_max[ch] = -1250;
      charges[ch] = 0;
      tstarts[ch] = -1;

      for (int s=0; s<1024; ++s) {
	if (waveforms[ch][s] < amplitudes_min[ch]) amplitudes_min[ch] = waveforms[ch][s];
	if (waveforms[ch][s] > amplitudes_max[ch]) amplitudes_max[ch] = waveforms[ch][s];}
      
      for (int s=0; s<baseline_sample; ++s)
	baselines[ch] += waveforms[ch][s];
      baselines[ch] /= baseline_sample;

      for (int s=charge_start_sample; s<charge_stop_sample; ++s)
	charges[ch] += waveforms[ch][s];

      charges[ch] -= (charge_stop_sample-charge_start_sample) * baselines[ch];
	
    } // for ch

    reco_tree->Fill();
	  
  } // for entry

  reco_tree->Write("", TObject::kOverwrite);
  current_file->Close();

}


/////////
// FIT //
/////////


const int NPE_MAX = 30;

double fit_poisson(double *x, double *par)
{
  const double X = x[0];
  
  const double NORM    = par[0];
  const double MEAN    = par[1];
  const double SIGMA_A = par[2]; // a
  const double SIGMA_B = par[3]; // b/sqrt(E)
  // const double SIGMA_C = par[4]; // c x E
  const double MU      = par[4];

  double VALUE = 0;
  
  for (int PE=1; PE<NPE_MAX; ++PE)
    {
      double POISSON = TMath::Poisson(PE, MU);

      double SIGMA = SIGMA_A + SIGMA_B*TMath::Sqrt(PE*1.0);
      double GAUSS = TMath::Gaus(X, MEAN*PE, SIGMA, true);

      VALUE += POISSON*GAUSS;
    }
  
  return NORM*VALUE;
}

double fit_poisson_free(double *x, double *par)
{
  const double X = x[0];
  
  const double MEAN    = par[0];
  const double SIGMA_A = par[1]; // a
  const double SIGMA_B = par[2]; // b/sqrt(E)
  // const double SIGMA_C = par[4]; // c x E

  const int NPE = par[3];
  
  double VALUE = 0;
  
  for (int PE=0; PE<NPE; ++PE)
    {
      double SIGMA = SIGMA_A + SIGMA_B*TMath::Sqrt(PE*1.0);
      double GAUSS = TMath::Gaus(X, MEAN*PE, SIGMA);

      VALUE += par[4+PE]*GAUSS;
    }
  
  return VALUE;
}


void fit ()
{
  const int npe = 28;
  
  TF1 *fitf = new TF1("fitf", fit_poisson_free, 0, 20000, 4+npe);

  fitf->SetParNames("MEAN", "SIGMA_a", "SIGMA_b", "NPE");
  for (int pe=0; pe<npe; ++pe) fitf->SetParName(4+pe, Form("N%02d",pe));

   fitf->SetParameters(500, 130, 20);
   fitf->FixParameter(3, npe);

   for (int pe=0; pe<npe; ++pe) fitf->SetParameter(4+pe, 50);

   fitf->SetLineColor(kRed);
   fitf->SetLineWidth(2);
   fitf->SetNpx(1000);
   
}
