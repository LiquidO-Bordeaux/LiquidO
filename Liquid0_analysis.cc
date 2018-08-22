
//////////////////////////////////////////////////////////////////////////:///////////////////////////////////////////////
//																														//
//		LIQUIDO_ANALYSIS																								//
//		by MSP & AP																										//
//----------------------------------------------------------------------------------------------------------------------//
//		2018-08-03																										//
//		Call Function:	spectre_SiPM()																					//
//																														//
//	*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*	//
//	*	- channel #01 -> SiPM #1 (µ-LiquidO position #0 (low)         							  					*	//
//	*	- channel #02 -> SiPM #2 (µ-LiquidO position #1 (middle)													*	//
//	*	- channel #03 -> SiPM #3 (µ-LiquidO position #2 (high)														*	//
//	*	- channel #04 -> PM 																						*	//
//	*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*	//
//																														//
//----------------------------------------------------------------------------------------------------------------------//
//		2018-08-04		Added:		µ error calculation																	//
//						Corrected: 	µ calculation -> int to float definition for "somme"							 	//
//						Modified: 	Call Function includes arguments -> spectre_SiPM(filename,E_cut,npe01,npe02,npe03)	//
//----------------------------------------------------------------------------------------------------------------------//
//		2018-08-05		Modified: 	Call Function includes arguments -> spectre_SiPM(filename,E_cut,npe01,npe02,npe03)	//
//----------------------------------------------------------------------------------------------------------------------//
//		2018-08-06		Added:		Write µ results on last Canvas														//
//						Modified:	New function (as a replacement of several lines) used in SiPM fits					//
//----------------------------------------------------------------------------------------------------------------------//
//		2018-08-07		Added:		Alternate calculation of mean µ														//
//----------------------------------------------------------------------------------------------------------------------//
//		2018-08-10		Added:		NoWash & Fibers Draw routines														//
//----------------------------------------------------------------------------------------------------------------------//
//		2018-08-13		Added:		Charge vs Time Draw routine 														//
//////////////////////////////////////////////////////////////////////////:///////////////////////////////////////////////

#include <math.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include <fstream>
#include <string>
#include <cstring>
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TApplication.h"
#include "TMultiGraph.h"
using namespace std;

//
//	Fit functions (from Manu)
//
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
//
//
//
TF1* fit (const int npe)
{
//  const int npe = 28;

  TF1 *fitf = new TF1("fitf", fit_poisson_free, 0, 20000, 4+npe);

  fitf->SetParNames("MEAN", "SIGMA_a", "SIGMA_b", "NPE");
  for (int pe=0; pe<npe; ++pe) fitf->SetParName(4+pe, Form("N%02d",pe));

   fitf->SetParameters(500, 130, 20);
   fitf->FixParameter(3, npe);

   for (int pe=0; pe<npe; ++pe) fitf->SetParameter(4+pe, 50);

   fitf->SetLineColor(kRed);
   fitf->SetLineWidth(2);
   fitf->SetNpx(1000);
   return fitf;
}

void delete_canvas(TCanvas* c1, TCanvas* c2, TCanvas* c3, TCanvas* c4, TCanvas* c99)
{
	c1->Close();
	c2->Close();
	c3->Close();
	c4->Close();
	c99->Close();
}

//
//	MEAN µ CALCULATION
//		- Weighted Arithmetic Mean, where:
//			somme			= Sum[i=5,npe+i-1)](amplitude(i))
//			amplitude 		= (Sum[i=5,npe+i-1)](amplitude(i)*µ(i)))/somme
//		- Error is computed as:
//			numerateur 		= Sum[i=5,npe+i-1)](error(i)*µ(i))
//			denominateur	= Sum[i=5,npe+i-1)](error(i))
//			incertitude 	= sqrt((numerateur/somme)^2+((amplitude*denominateur)/somme)^2)
//			(Note: this is a crude estimation)
//
//	New vector used in SiPM fits (2018-08-06 by Axel)
//
//

vector<double> *calcul_amplitude(const int npe, TF1* myfit)
{
  double amplitude 	  	= 0.;
  double somme 		 	= 0.;
  double numerateur 	= 0.;
  double denominateur	= 0.;
  double incertitude 	= 0.;
  for (int i=5; i<npe+4; i++)
  {
  	amplitude    += myfit->GetParameter(i)*(i-4);
  	somme 		 += myfit->GetParameter(i);
  	numerateur   += myfit->GetParError(i)*(i-4);
  	denominateur += myfit->GetParError(i);
  }
  amplitude /= somme;
  incertitude = sqrt((numerateur/somme)*(numerateur/somme)+((amplitude*denominateur)/somme)*((amplitude*denominateur)/somme));
  vector<double>* amp = new vector<double>;
  amp->clear();
  amp->push_back(amplitude);
  amp->push_back(incertitude);
  return amp;
}

//																														//
//	*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*	//
//	*																												*	//
//	*	         	####   #######	 ####   #####   #######                     									*	//
//	*			   #    #  #  #  #	#	 #  #    #  #  #  #															*	//
//	*			 	# 		  #	    #	 #  #    #     #															*	//
//	*				  #       #     ######  #####      #															*	//
//  *              #	#     #     #    #  #  #       #															*	//
//  *				####      #     #    #  #    #     #															*	//
//	*																												*	//
//	*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*	//
//																														//

void spectre_SiPM(const char *filename, const int E_cut, const int npe01, const int npe02, const int npe03)
{
	TFile *_file0 		= new TFile(filename,"READ");
	TTree* reco_tree 	= (TTree*)_file0->Get("reco_tree");

	float charge01, charge02, charge03, charge04;

	reco_tree->SetBranchAddress("ch01_charge", &charge01);
	reco_tree->SetBranchAddress("ch02_charge", &charge02);
	reco_tree->SetBranchAddress("ch03_charge", &charge03);
	reco_tree->SetBranchAddress("ch04_charge", &charge04);

	TH1F *charge01_histo  = new TH1F("charge01_histo", "SiPM 0", 400, 0, 20000);
	TH1F *charge02_histo  = new TH1F("charge02_histo", "SiPM 1", 400, 0, 20000);
	TH1F *charge03_histo  = new TH1F("charge03_histo", "SiPM 2", 400, 0, 20000);
	TH1F *charge04_histo  = new TH1F("charge04_histo", "PM"    , 1000, 0, 30000);

	for (int i=0; i<reco_tree->GetEntries(); i++)
	{
		reco_tree->GetEntry(i);
		if (-charge04>E_cut)
		{
			charge01_histo->Fill(charge01);
			charge02_histo->Fill(charge02);
			charge03_histo->Fill(charge03);
		}
		charge04_histo->Fill(-charge04);
	}

//
//	Fit SiPM position 0
//
TF1* myfit01 = fit(npe01);
TCanvas *c1 = new TCanvas;
c1->SetLogy();
charge01_histo->Draw();
charge01_histo->Fit("fitf","q");
myfit01->Draw("same");
vector<double>* amplitude01 = calcul_amplitude(npe01, myfit01);
cout << "µ position 0 = " << amplitude01->at(0) << " +- " << amplitude01->at(1) << " (chi2 = " << myfit01->GetChisquare()/myfit01->GetNDF() << ")" << endl;
//
//	Fit SiPM position 1
//
TF1* myfit02 	= fit(npe02);
TCanvas *c2 	= new TCanvas;
c2->SetLogy();
charge02_histo->Draw();
charge02_histo->Fit("fitf","q");
myfit02->Draw("same");
vector<double>* amplitude02 = calcul_amplitude(npe02, myfit02);
cout << "µ position 1 = " << amplitude02->at(0) << " +- " << amplitude02->at(1) << " (chi2 = " << myfit02->GetChisquare()/myfit02->GetNDF() << ")" << endl;
//
//	Fit SiPM position 3
//
TF1* myfit03 	= fit(npe03);
TCanvas *c3 	= new TCanvas;
c3->SetLogy();
charge03_histo->Draw();
charge03_histo->Fit("fitf","q");
myfit03->Draw("same");
vector<double>* amplitude03 = calcul_amplitude(npe03, myfit03);
cout << "µ position 2 = " << amplitude03->at(0) << " +- " << amplitude03->at(1) << " (chi2 = " << myfit03->GetChisquare()/myfit03->GetNDF() << ")" << endl;
//
//	PM
//
TCanvas *c4 = new TCanvas;
c4->SetLogy();
charge04_histo->Draw();

TF1* fit_pm = new TF1("fit_pm","gaus",2000,20000);
charge04_histo->Fit("fit_pm","RQ");
for (int i=0; i<3; i++)
{
	double mean_gaussian  = fit_pm->GetParameter(1);
	double sigma_gaussian = fit_pm->GetParameter(2);	
	fit_pm->SetRange(mean_gaussian-0.8*sigma_gaussian, mean_gaussian+2.5*sigma_gaussian );
	charge04_histo->Fit("fit_pm","RQ");
}
fit_pm->Draw("same");
cout << "Mean PM = " << fit_pm->GetParameter(1) << " +- " << fit_pm->GetParError(1) << endl;


//
//	Alternate calculation of mean µ (2018-08-07 AP)
//
//cout << "Resultats Valeur moyenne" << endl;
//double value01 		= charge01_histo->GetMean()/myfit01->GetParameter(0);
//double value01_err 	= value01*sqrt(TMath::Power(charge01_histo->GetMeanError()/charge01_histo->GetMean(),2)+TMath::Power(myfit01->GetParError(0)/myfit01->GetParameter(0),2));
//cout << "µ position 0 = " << value01 << " +- " << value01_err << endl;
//
//double value02 		= charge02_histo->GetMean()/myfit02->GetParameter(0);
//double value02_err 	= value02*sqrt(TMath::Power(charge02_histo->GetMeanError()/charge02_histo->GetMean(),2)+TMath::Power(myfit02->GetParError(0)/myfit02->GetParameter(0),2));
//cout << "µ position 0 = " << value02 << " +- " << value02_err << endl;
//
//double value03 		= charge03_histo->GetMean()/myfit03->GetParameter(0);
//double value03_err 	= value03*sqrt(TMath::Power(charge03_histo->GetMeanError()/charge03_histo->GetMean(),2)+TMath::Power(myfit03->GetParError(0)/myfit03->GetParameter(0),2));
//cout << "µ position 0 = " << value03 << " +- " << value03_err << endl;

//
//	Option not to show statistics on Canvas
//
gStyle->SetOptStat(0);
//
TCanvas *c99 = new TCanvas;
c99->SetLogy();
charge03_histo->SetLineColor(kViolet);
charge03_histo->SetTitle(filename);
charge03_histo->GetXaxis()->SetTitle("charge");
charge03_histo->GetYaxis()->SetTitle("Intensite");
charge03_histo->GetYaxis()->SetRangeUser(0.1,1000.);
charge03_histo->Draw();
charge02_histo->SetLineColor(kOrange);
charge02_histo->Draw("same");
charge01_histo->SetLineColor(kSpring);
charge01_histo->Draw("same");

//
//	Adds legend to the Canvas
//
TLegend *l = new TLegend(0.6,0.7,0.9,0.9);
l->AddEntry(charge01_histo, Form("charge SiPM 0; #mu = %.02f +- %.02f",amplitude01->at(0),amplitude01->at(1)));
l->AddEntry(charge02_histo, Form("charge SiPM 1; #mu = %.02f +- %.02f",amplitude02->at(0),amplitude02->at(1)));
l->AddEntry(charge03_histo, Form("charge SiPM 2; #mu = %.02f +- %.02f",amplitude03->at(0),amplitude03->at(1)));
l->Draw("same");

//
//c1->Modified(); c1->Update();
//c99->Modified(); c99->Update();
//

//
//	Next 2 lines for closing all opened Canvas. Waits for a keyboard input. 
//
//getchar();
//delete_canvas(c1,c2,c3,c4,c99);
//


//
//	NoWash and Fibers draw routines (2018-08-10 Axel)
//		- draw_graphs_nowash()
//		- draw_graphs_fibre()
//
}
void draw_graphs_nowash()
{
  const int npts_nowash_00 = 5;
  const int npts_nowash_10 = 9;
  const int npts_nowash_15 = 10;
  const int npts_nowash_20 = 9;

  double mu_fibre_0_nowash_00[npts_nowash_00] = {6.48,4.39,6.25,5.84,6.13};
  double mu_fibre_1_nowash_00[npts_nowash_00] = {3.78,3.60,3.95,3.99,3.70};
  double mu_fibre_2_nowash_00[npts_nowash_00] = {3.41,3.45,3.31,3.39,3.37};
  double mu_fibre_0_nowash_10[npts_nowash_10] = {8.61,8.42,8.36,7.69,8.15,8.08,5.64,6.51,6.17};
  double mu_fibre_1_nowash_10[npts_nowash_10] = {3.94,3.79,3.57,3.92,3.64,3.67,3.56,3.42,3.47};
  double mu_fibre_2_nowash_10[npts_nowash_10] = {3.23,2.61,2.76,3.28,3.38,3.44,3.56,3.41,3.75};
  double mu_fibre_0_nowash_15[npts_nowash_15] = {8.94,8.61,6.41,6.19,6.48,6.15,6.11,6.90,7.04,7.13};
  double mu_fibre_1_nowash_15[npts_nowash_15] = {4.15,4.44,5.20,5.05,5.16,5.12,5.08,5.08,5.08,5.22};
  double mu_fibre_2_nowash_15[npts_nowash_15] = {3.75,3.76,3.30,3.42,3.50,3.41,3.53,3.43,3.42,3.49};
  double mu_fibre_0_nowash_20[npts_nowash_20] = {7.30,7.28,7.58,9.20,7.59,7.49,7.32,7.49};
  double mu_fibre_1_nowash_20[npts_nowash_20] = {4.73,5.13,5.22,5.20,4.49,4.79,3.75,3.01};
  double mu_fibre_2_nowash_20[npts_nowash_20] = {2.31,2.28,2.26,3.00,3.00,3.00,2.93,3.01};

  double mu_err_fibre_0_nowash_00[npts_nowash_00] = {0.36,0.20,0.30,0.28,0.31};
  double mu_err_fibre_1_nowash_00[npts_nowash_00] = {0.18,0.15,0.15,0.17,0.14};
  double mu_err_fibre_2_nowash_00[npts_nowash_00] = {0.16,0.15,0.14,0.15,0.14};
  double mu_err_fibre_0_nowash_10[npts_nowash_10] = {0.38,0.42,0.42,0.41,0.41,0.40,0.27,0.32,0.29};
  double mu_err_fibre_1_nowash_10[npts_nowash_10] = {0.15,0.15,0.13,0.17,0.15,0.14,0.15,0.13,0.15};
  double mu_err_fibre_2_nowash_10[npts_nowash_10] = {0.12,0.09,0.11,0.12,0.14,0.15,0.16,0.15,0.17};
  double mu_err_fibre_0_nowash_15[npts_nowash_15] = {0.49,0.42,0.29,0.28,0.31,0.29,0.30,0.32,0.36,0.36};
  double mu_err_fibre_1_nowash_15[npts_nowash_15] = {0.17,0.19,0.24,0.19,0.22,0.21,0.20,0.20,0.21,0.22};
  double mu_err_fibre_2_nowash_15[npts_nowash_15] = {0.16,0.16,0.13,0.14,0.14,0.14,0.16,0.15,0.15,0.15};
  double mu_err_fibre_0_nowash_20[npts_nowash_20] = {0.36,0.35,0.38,0.51,0.39,0.39,0.37,0.38};
  double mu_err_fibre_1_nowash_20[npts_nowash_20] = {0.21,0.23,0.23,0.22,0.19,0.20,0.15,0.16};
  double mu_err_fibre_2_nowash_20[npts_nowash_20] = {0.11,0.10,0.09,0.13,0.12,0.20,0.13,0.12};

  double temperature_nowash_00[npts_nowash_00] = {21.5,24.3,26.1,27.0,30.0};
  double temperature_nowash_10[npts_nowash_10] = {7.6,13.5,18.7,20.7,23.0,23.6,25.2,25.7,26.2};
  double temperature_nowash_15[npts_nowash_15] = {1.9,7.3,13.5,16.5,18.7,20.3,22.9,23.7,24.6,25.1};
  double temperature_nowash_20[npts_nowash_20] = {1.4,8.3,11.5,14.8,16.3,18.8,20.3,21.2};

  double temperature_err_nowash_00[npts_nowash_00] = {0.1,0.1,0.1,0.1,0.1};
  double temperature_err_nowash_10[npts_nowash_10] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  double temperature_err_nowash_15[npts_nowash_15] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  double temperature_err_nowash_20[npts_nowash_20] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};

  double charge_PM_nowash_00[npts_nowash_00] = {2.4729,2.4231,2.3806,2.4086,2.4004};
  double charge_PM_nowash_10[npts_nowash_10] = {1.1459,1.2980,1.3758,1.3953,1.5897,1.6549,17.116,1.7869,1.7874};
  double charge_PM_nowash_15[npts_nowash_15] = {1.1403,1.2259,1.3421,1.4201,1.4697,1.5003,1.5482,16.583,1.6520,1.6457};
  double charge_PM_nowash_20[npts_nowash_20] = {0.5549,0.5863,0.6340,0.6809,0.7093,0.7299,0.7446,0.7561};

  double charge_err_PM_nowash_00[npts_nowash_00] = {0.002,0.002,0.0011,0.002,0.002};
  double charge_err_PM_nowash_10[npts_nowash_10] = {0.0014,0.0014,0.0014,0.0015,0.0015,0.0017,0.0016,0.0017,0.0017};
  double charge_err_PM_nowash_15[npts_nowash_15] = {0.0014,0.0014,0.0015,0.0015,0.0016,0.0016,0.0016,0.0016,0.0016,0.0015};
  double charge_err_PM_nowash_20[npts_nowash_20] = {0.0009,0.0009,0.0010,0.0011,0.0010,0.0011,0.0011,0.0011};

  TGraphErrors* mu_fibre_0_nowash_00_vs_temperature = new TGraphErrors(npts_nowash_00, temperature_nowash_00, mu_fibre_0_nowash_00, temperature_err_nowash_00, mu_err_fibre_0_nowash_00);
  TGraphErrors* mu_fibre_1_nowash_00_vs_temperature = new TGraphErrors(npts_nowash_00, temperature_nowash_00, mu_fibre_1_nowash_00, temperature_err_nowash_00, mu_err_fibre_1_nowash_00);
  TGraphErrors* mu_fibre_2_nowash_00_vs_temperature = new TGraphErrors(npts_nowash_00, temperature_nowash_00, mu_fibre_2_nowash_00, temperature_err_nowash_00, mu_err_fibre_2_nowash_00);
  TGraphErrors* mu_fibre_0_nowash_10_vs_temperature = new TGraphErrors(npts_nowash_10, temperature_nowash_10, mu_fibre_0_nowash_10, temperature_err_nowash_10, mu_err_fibre_0_nowash_10);
  TGraphErrors* mu_fibre_1_nowash_10_vs_temperature = new TGraphErrors(npts_nowash_10, temperature_nowash_10, mu_fibre_1_nowash_10, temperature_err_nowash_10, mu_err_fibre_1_nowash_10);
  TGraphErrors* mu_fibre_2_nowash_10_vs_temperature = new TGraphErrors(npts_nowash_10, temperature_nowash_10, mu_fibre_2_nowash_10, temperature_err_nowash_10, mu_err_fibre_2_nowash_10);
  TGraphErrors* mu_fibre_0_nowash_15_vs_temperature = new TGraphErrors(npts_nowash_15, temperature_nowash_15, mu_fibre_0_nowash_15, temperature_err_nowash_15, mu_err_fibre_0_nowash_15);
  TGraphErrors* mu_fibre_1_nowash_15_vs_temperature = new TGraphErrors(npts_nowash_15, temperature_nowash_15, mu_fibre_1_nowash_15, temperature_err_nowash_15, mu_err_fibre_1_nowash_15);
  TGraphErrors* mu_fibre_2_nowash_15_vs_temperature = new TGraphErrors(npts_nowash_15, temperature_nowash_15, mu_fibre_2_nowash_15, temperature_err_nowash_15, mu_err_fibre_2_nowash_15);
  TGraphErrors* mu_fibre_0_nowash_20_vs_temperature = new TGraphErrors(npts_nowash_20, temperature_nowash_20, mu_fibre_0_nowash_20, temperature_err_nowash_20, mu_err_fibre_0_nowash_20);
  TGraphErrors* mu_fibre_1_nowash_20_vs_temperature = new TGraphErrors(npts_nowash_20, temperature_nowash_20, mu_fibre_1_nowash_20, temperature_err_nowash_20, mu_err_fibre_1_nowash_20);
  TGraphErrors* mu_fibre_2_nowash_20_vs_temperature = new TGraphErrors(npts_nowash_20, temperature_nowash_20, mu_fibre_2_nowash_20, temperature_err_nowash_20, mu_err_fibre_2_nowash_20);

  TGraphErrors* charge_PM_nowash_00_vs_temperature = new TGraphErrors(npts_nowash_00, temperature_nowash_00, charge_PM_nowash_00, temperature_err_nowash_00, charge_err_PM_nowash_00);
  TGraphErrors* charge_PM_nowash_10_vs_temperature = new TGraphErrors(npts_nowash_10, temperature_nowash_10, charge_PM_nowash_10, temperature_err_nowash_10, charge_err_PM_nowash_10);
  TGraphErrors* charge_PM_nowash_15_vs_temperature = new TGraphErrors(npts_nowash_15, temperature_nowash_15, charge_PM_nowash_15, temperature_err_nowash_15, charge_err_PM_nowash_15);
  TGraphErrors* charge_PM_nowash_20_vs_temperature = new TGraphErrors(npts_nowash_20, temperature_nowash_20, charge_PM_nowash_20, temperature_err_nowash_20, charge_err_PM_nowash_20);

  new TCanvas;
  mu_fibre_0_nowash_00_vs_temperature->GetYaxis()->SetRangeUser(0,10);
  mu_fibre_0_nowash_00_vs_temperature->GetYaxis()->SetTitle("#mu moyen");
  mu_fibre_0_nowash_00_vs_temperature->GetXaxis()->SetTitle("Temperature nowash (^{o}C)");
  mu_fibre_0_nowash_00_vs_temperature->SetTitle("nowash 0%");
  mu_fibre_0_nowash_00_vs_temperature->SetLineColor(kRed);
  mu_fibre_0_nowash_00_vs_temperature->SetMarkerColor(kRed);
  mu_fibre_0_nowash_00_vs_temperature->SetMarkerStyle(7);
  mu_fibre_0_nowash_00_vs_temperature->SetMarkerSize(3);
  mu_fibre_0_nowash_00_vs_temperature->Draw("ap");
  mu_fibre_1_nowash_00_vs_temperature->SetLineColor(kViolet);
  mu_fibre_1_nowash_00_vs_temperature->SetMarkerColor(kViolet);
  mu_fibre_1_nowash_00_vs_temperature->SetMarkerStyle(7);
  mu_fibre_1_nowash_00_vs_temperature->SetMarkerSize(3);
  mu_fibre_1_nowash_00_vs_temperature->Draw("psame");
  mu_fibre_2_nowash_00_vs_temperature->SetLineColor(kBlue);
  mu_fibre_2_nowash_00_vs_temperature->SetMarkerColor(kBlue);
  mu_fibre_2_nowash_00_vs_temperature->SetMarkerStyle(7);
  mu_fibre_2_nowash_00_vs_temperature->SetMarkerSize(3);
  mu_fibre_2_nowash_00_vs_temperature->Draw("psame");
  charge_PM_nowash_00_vs_temperature->SetLineColor(kBlack);
  charge_PM_nowash_00_vs_temperature->SetMarkerColor(kBlack);
  charge_PM_nowash_00_vs_temperature->SetMarkerStyle(7);
  charge_PM_nowash_00_vs_temperature->SetMarkerSize(3);
  charge_PM_nowash_00_vs_temperature->Draw("psame");
  TLegend* l_nowash_00_vs_temperature = new TLegend(0.65,0.7,0.9,0.9);
  l_nowash_00_vs_temperature->AddEntry(mu_fibre_0_nowash_00_vs_temperature, "#mu fibre 0");
  l_nowash_00_vs_temperature->AddEntry(mu_fibre_1_nowash_00_vs_temperature, "#mu fibre 1");
  l_nowash_00_vs_temperature->AddEntry(mu_fibre_2_nowash_00_vs_temperature, "#mu fibre 2");
  l_nowash_00_vs_temperature->AddEntry(charge_PM_nowash_00_vs_temperature , "charge PM (x10^{4})");
  l_nowash_00_vs_temperature->Draw("lsame");

  new TCanvas;
  mu_fibre_0_nowash_10_vs_temperature->GetYaxis()->SetRangeUser(0,10);
  mu_fibre_0_nowash_10_vs_temperature->GetYaxis()->SetTitle("#mu moyen");
  mu_fibre_0_nowash_10_vs_temperature->GetXaxis()->SetTitle("Temperature nowash (^{o}C)");
  mu_fibre_0_nowash_10_vs_temperature->SetTitle("nowash 10%");
  mu_fibre_0_nowash_10_vs_temperature->SetLineColor(kRed);
  mu_fibre_0_nowash_10_vs_temperature->SetMarkerColor(kRed);
  mu_fibre_0_nowash_10_vs_temperature->SetMarkerStyle(7);
  mu_fibre_0_nowash_10_vs_temperature->SetMarkerSize(3);
  mu_fibre_0_nowash_10_vs_temperature->Draw("ap");
  mu_fibre_1_nowash_10_vs_temperature->SetLineColor(kViolet);
  mu_fibre_1_nowash_10_vs_temperature->SetMarkerColor(kViolet);
  mu_fibre_1_nowash_10_vs_temperature->SetMarkerStyle(7);
  mu_fibre_1_nowash_10_vs_temperature->SetMarkerSize(3);
  mu_fibre_1_nowash_10_vs_temperature->Draw("psame");
  mu_fibre_2_nowash_10_vs_temperature->SetLineColor(kBlue);
  mu_fibre_2_nowash_10_vs_temperature->SetMarkerColor(kBlue);
  mu_fibre_2_nowash_10_vs_temperature->SetMarkerStyle(7);
  mu_fibre_2_nowash_10_vs_temperature->SetMarkerSize(3);
  mu_fibre_2_nowash_10_vs_temperature->Draw("psame");
  charge_PM_nowash_10_vs_temperature->SetLineColor(kBlack);
  charge_PM_nowash_10_vs_temperature->SetMarkerColor(kBlack);
  charge_PM_nowash_10_vs_temperature->SetMarkerStyle(7);
  charge_PM_nowash_10_vs_temperature->SetMarkerSize(3);
  charge_PM_nowash_10_vs_temperature->Draw("psame");
  TLegend* l_nowash_10_vs_temperature = new TLegend(0.65,0.7,0.9,0.9);
  l_nowash_10_vs_temperature->AddEntry(mu_fibre_0_nowash_10_vs_temperature, "#mu fibre 0");
  l_nowash_10_vs_temperature->AddEntry(mu_fibre_1_nowash_10_vs_temperature, "#mu fibre 1");
  l_nowash_10_vs_temperature->AddEntry(mu_fibre_2_nowash_10_vs_temperature, "#mu fibre 2");
  l_nowash_10_vs_temperature->AddEntry(charge_PM_nowash_10_vs_temperature , "charge PM (x10^{4})");
  l_nowash_10_vs_temperature->Draw("lsame");

  new TCanvas;
  mu_fibre_0_nowash_15_vs_temperature->GetYaxis()->SetRangeUser(0,10);
  mu_fibre_0_nowash_15_vs_temperature->GetYaxis()->SetTitle("#mu moyen");
  mu_fibre_0_nowash_15_vs_temperature->GetXaxis()->SetTitle("Temperature nowash (^{o}C)");
  mu_fibre_0_nowash_15_vs_temperature->SetTitle("nowash 15%");
  mu_fibre_0_nowash_15_vs_temperature->SetLineColor(kRed);
  mu_fibre_0_nowash_15_vs_temperature->SetMarkerColor(kRed);
  mu_fibre_0_nowash_15_vs_temperature->SetMarkerStyle(7);
  mu_fibre_0_nowash_15_vs_temperature->SetMarkerSize(3);
  mu_fibre_0_nowash_15_vs_temperature->Draw("ap");
  mu_fibre_1_nowash_15_vs_temperature->SetLineColor(kViolet);
  mu_fibre_1_nowash_15_vs_temperature->SetMarkerColor(kViolet);
  mu_fibre_1_nowash_15_vs_temperature->SetMarkerStyle(7);
  mu_fibre_1_nowash_15_vs_temperature->SetMarkerSize(3);
  mu_fibre_1_nowash_15_vs_temperature->Draw("psame");
  mu_fibre_2_nowash_15_vs_temperature->SetLineColor(kBlue);
  mu_fibre_2_nowash_15_vs_temperature->SetMarkerColor(kBlue);
  mu_fibre_2_nowash_15_vs_temperature->SetMarkerStyle(7);
  mu_fibre_2_nowash_15_vs_temperature->SetMarkerSize(3);
  mu_fibre_2_nowash_15_vs_temperature->Draw("psame");
  charge_PM_nowash_15_vs_temperature->SetLineColor(kBlack);
  charge_PM_nowash_15_vs_temperature->SetMarkerColor(kBlack);
  charge_PM_nowash_15_vs_temperature->SetMarkerStyle(7);
  charge_PM_nowash_15_vs_temperature->SetMarkerSize(3);
  charge_PM_nowash_15_vs_temperature->Draw("psame");
  TLegend* l_nowash_15_vs_temperature = new TLegend(0.65,0.7,0.9,0.9);
  l_nowash_15_vs_temperature->AddEntry(mu_fibre_0_nowash_15_vs_temperature, "#mu fibre 0");
  l_nowash_15_vs_temperature->AddEntry(mu_fibre_1_nowash_15_vs_temperature, "#mu fibre 1");
  l_nowash_15_vs_temperature->AddEntry(mu_fibre_2_nowash_15_vs_temperature, "#mu fibre 2");
  l_nowash_15_vs_temperature->AddEntry(charge_PM_nowash_15_vs_temperature , "charge PM (x10^{4})");
  l_nowash_15_vs_temperature->Draw("lsame");

  new TCanvas;
  mu_fibre_0_nowash_20_vs_temperature->GetYaxis()->SetRangeUser(0,10);
  mu_fibre_0_nowash_20_vs_temperature->GetYaxis()->SetTitle("#mu moyen");
  mu_fibre_0_nowash_20_vs_temperature->GetXaxis()->SetTitle("Temperature nowash (^{o}C)");
  mu_fibre_0_nowash_20_vs_temperature->SetTitle("nowash 20%");
  mu_fibre_0_nowash_20_vs_temperature->SetLineColor(kRed);
  mu_fibre_0_nowash_20_vs_temperature->SetMarkerColor(kRed);
  mu_fibre_0_nowash_20_vs_temperature->SetMarkerStyle(7);
  mu_fibre_0_nowash_20_vs_temperature->SetMarkerSize(3);
  mu_fibre_0_nowash_20_vs_temperature->Draw("ap");
  mu_fibre_1_nowash_20_vs_temperature->SetLineColor(kViolet);
  mu_fibre_1_nowash_20_vs_temperature->SetMarkerColor(kViolet);
  mu_fibre_1_nowash_20_vs_temperature->SetMarkerStyle(7);
  mu_fibre_1_nowash_20_vs_temperature->SetMarkerSize(3);
  mu_fibre_1_nowash_20_vs_temperature->Draw("psame");
  mu_fibre_2_nowash_20_vs_temperature->SetLineColor(kBlue);
  mu_fibre_2_nowash_20_vs_temperature->SetMarkerColor(kBlue);
  mu_fibre_2_nowash_20_vs_temperature->SetMarkerStyle(7);
  mu_fibre_2_nowash_20_vs_temperature->SetMarkerSize(3);
  mu_fibre_2_nowash_20_vs_temperature->Draw("psame");
  charge_PM_nowash_20_vs_temperature->SetLineColor(kBlack);
  charge_PM_nowash_20_vs_temperature->SetMarkerColor(kBlack);
  charge_PM_nowash_20_vs_temperature->SetMarkerStyle(7);
  charge_PM_nowash_20_vs_temperature->SetMarkerSize(3);
  charge_PM_nowash_20_vs_temperature->Draw("psame");
  TLegend* l_nowash_20_vs_temperature = new TLegend(0.65,0.7,0.9,0.9);
  l_nowash_20_vs_temperature->AddEntry(mu_fibre_0_nowash_20_vs_temperature, "#mu fibre 0");
  l_nowash_20_vs_temperature->AddEntry(mu_fibre_1_nowash_20_vs_temperature, "#mu fibre 1");
  l_nowash_20_vs_temperature->AddEntry(mu_fibre_2_nowash_20_vs_temperature, "#mu fibre 2");
  l_nowash_20_vs_temperature->AddEntry(charge_PM_nowash_20_vs_temperature , "charge PM (x10^{4})");
  l_nowash_20_vs_temperature->Draw("lsame");
}

void draw_graphs_fibre()
{
  const int npts_nowash_00 = 5;
  const int npts_nowash_10 = 9;
  const int npts_nowash_15 = 10;
  const int npts_nowash_20 = 9;

  double mu_fibre_0_nowash_00[npts_nowash_00] = {6.48,4.39,6.25,5.84,6.13};
  double mu_fibre_1_nowash_00[npts_nowash_00] = {3.78,3.60,3.95,3.99,3.70};
  double mu_fibre_2_nowash_00[npts_nowash_00] = {3.41,3.45,3.31,3.39,3.37};
  double mu_fibre_0_nowash_10[npts_nowash_10] = {8.61,8.42,8.36,7.69,8.15,8.08,5.64,6.51,6.17};
  double mu_fibre_1_nowash_10[npts_nowash_10] = {3.94,3.79,3.57,3.92,3.64,3.67,3.56,3.42,3.47};
  double mu_fibre_2_nowash_10[npts_nowash_10] = {3.23,2.61,2.76,3.28,3.38,3.44,3.56,3.41,3.75};
  double mu_fibre_0_nowash_15[npts_nowash_15] = {8.94,8.61,6.41,6.19,6.48,6.15,6.11,6.90,7.04,7.13};
  double mu_fibre_1_nowash_15[npts_nowash_15] = {4.15,4.44,5.20,5.05,5.16,5.12,5.08,5.08,5.08,5.22};
  double mu_fibre_2_nowash_15[npts_nowash_15] = {3.75,3.76,3.30,3.42,3.50,3.41,3.53,3.43,3.42,3.49};
  double mu_fibre_0_nowash_20[npts_nowash_20] = {7.30,7.28,7.58,9.20,7.59,7.49,7.32,7.49};
  double mu_fibre_1_nowash_20[npts_nowash_20] = {4.73,5.13,5.22,5.20,4.49,4.79,3.75,3.01};
  double mu_fibre_2_nowash_20[npts_nowash_20] = {2.31,2.28,2.26,3.00,3.00,3.00,2.93,3.01};

  double mu_err_fibre_0_nowash_00[npts_nowash_00] = {0.36,0.20,0.30,0.28,0.31};
  double mu_err_fibre_1_nowash_00[npts_nowash_00] = {0.18,0.15,0.15,0.17,0.14};
  double mu_err_fibre_2_nowash_00[npts_nowash_00] = {0.16,0.15,0.14,0.15,0.14};
  double mu_err_fibre_0_nowash_10[npts_nowash_10] = {0.38,0.42,0.42,0.41,0.41,0.40,0.27,0.32,0.29};
  double mu_err_fibre_1_nowash_10[npts_nowash_10] = {0.15,0.15,0.13,0.17,0.15,0.14,0.15,0.13,0.15};
  double mu_err_fibre_2_nowash_10[npts_nowash_10] = {0.12,0.09,0.11,0.12,0.14,0.15,0.16,0.15,0.17};
  double mu_err_fibre_0_nowash_15[npts_nowash_15] = {0.49,0.42,0.29,0.28,0.31,0.29,0.30,0.32,0.36,0.36};
  double mu_err_fibre_1_nowash_15[npts_nowash_15] = {0.17,0.19,0.24,0.19,0.22,0.21,0.20,0.20,0.21,0.22};
  double mu_err_fibre_2_nowash_15[npts_nowash_15] = {0.16,0.16,0.13,0.14,0.14,0.14,0.16,0.15,0.15,0.15};
  double mu_err_fibre_0_nowash_20[npts_nowash_20] = {0.36,0.35,0.38,0.51,0.39,0.39,0.37,0.38};
  double mu_err_fibre_1_nowash_20[npts_nowash_20] = {0.21,0.23,0.23,0.22,0.19,0.20,0.15,0.16};
  double mu_err_fibre_2_nowash_20[npts_nowash_20] = {0.11,0.10,0.09,0.13,0.12,0.20,0.13,0.12};

  double temperature_nowash_00[npts_nowash_00] = {21.5,24.3,26.1,27.0,30.0};
  double temperature_nowash_10[npts_nowash_10] = {7.6,13.5,18.7,20.7,23.0,23.6,25.2,25.7,26.2};
  double temperature_nowash_15[npts_nowash_15] = {1.9,7.3,13.5,16.5,18.7,20.3,22.9,23.7,24.6,25.1};
  double temperature_nowash_20[npts_nowash_20] = {1.4,8.3,11.5,14.8,16.3,18.8,20.3,21.2};

  double temperature_err_nowash_00[npts_nowash_00] = {0.1,0.1,0.1,0.1,0.1};
  double temperature_err_nowash_10[npts_nowash_10] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  double temperature_err_nowash_15[npts_nowash_15] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  double temperature_err_nowash_20[npts_nowash_20] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};

  double charge_PM_nowash_00[npts_nowash_00] = {2.4729,2.4231,2.3806,2.4086,2.4004};
  double charge_PM_nowash_10[npts_nowash_10] = {1.1459,1.2980,1.3758,1.3953,1.5897,1.6549,17.116,1.7869,1.7874};
  double charge_PM_nowash_15[npts_nowash_15] = {1.1403,1.2259,1.3421,1.4201,1.4697,1.5003,1.5482,16.583,1.6520,1.6457};
  double charge_PM_nowash_20[npts_nowash_20] = {0.5549,0.5863,0.6340,0.6809,0.7093,0.7299,0.7446,0.7561};

  double charge_err_PM_nowash_00[npts_nowash_00] = {0.002,0.002,0.0011,0.002,0.002};
  double charge_err_PM_nowash_10[npts_nowash_10] = {0.0014,0.0014,0.0014,0.0015,0.0015,0.0017,0.0016,0.0017,0.0017};
  double charge_err_PM_nowash_15[npts_nowash_15] = {0.0014,0.0014,0.0015,0.0015,0.0016,0.0016,0.0016,0.0016,0.0016,0.0015};
  double charge_err_PM_nowash_20[npts_nowash_20] = {0.0009,0.0009,0.0010,0.0011,0.0010,0.0011,0.0011,0.0011};

  TGraphErrors* mu_fibre_0_nowash_00_vs_temperature = new TGraphErrors(npts_nowash_00, temperature_nowash_00, mu_fibre_0_nowash_00, temperature_err_nowash_00, mu_err_fibre_0_nowash_00);
  TGraphErrors* mu_fibre_1_nowash_00_vs_temperature = new TGraphErrors(npts_nowash_00, temperature_nowash_00, mu_fibre_1_nowash_00, temperature_err_nowash_00, mu_err_fibre_1_nowash_00);
  TGraphErrors* mu_fibre_2_nowash_00_vs_temperature = new TGraphErrors(npts_nowash_00, temperature_nowash_00, mu_fibre_2_nowash_00, temperature_err_nowash_00, mu_err_fibre_2_nowash_00);
  TGraphErrors* mu_fibre_0_nowash_10_vs_temperature = new TGraphErrors(npts_nowash_10, temperature_nowash_10, mu_fibre_0_nowash_10, temperature_err_nowash_10, mu_err_fibre_0_nowash_10);
  TGraphErrors* mu_fibre_1_nowash_10_vs_temperature = new TGraphErrors(npts_nowash_10, temperature_nowash_10, mu_fibre_1_nowash_10, temperature_err_nowash_10, mu_err_fibre_1_nowash_10);
  TGraphErrors* mu_fibre_2_nowash_10_vs_temperature = new TGraphErrors(npts_nowash_10, temperature_nowash_10, mu_fibre_2_nowash_10, temperature_err_nowash_10, mu_err_fibre_2_nowash_10);
  TGraphErrors* mu_fibre_0_nowash_15_vs_temperature = new TGraphErrors(npts_nowash_15, temperature_nowash_15, mu_fibre_0_nowash_15, temperature_err_nowash_15, mu_err_fibre_0_nowash_15);
  TGraphErrors* mu_fibre_1_nowash_15_vs_temperature = new TGraphErrors(npts_nowash_15, temperature_nowash_15, mu_fibre_1_nowash_15, temperature_err_nowash_15, mu_err_fibre_1_nowash_15);
  TGraphErrors* mu_fibre_2_nowash_15_vs_temperature = new TGraphErrors(npts_nowash_15, temperature_nowash_15, mu_fibre_2_nowash_15, temperature_err_nowash_15, mu_err_fibre_2_nowash_15);
  TGraphErrors* mu_fibre_0_nowash_20_vs_temperature = new TGraphErrors(npts_nowash_20, temperature_nowash_20, mu_fibre_0_nowash_20, temperature_err_nowash_20, mu_err_fibre_0_nowash_20);
  TGraphErrors* mu_fibre_1_nowash_20_vs_temperature = new TGraphErrors(npts_nowash_20, temperature_nowash_20, mu_fibre_1_nowash_20, temperature_err_nowash_20, mu_err_fibre_1_nowash_20);
  TGraphErrors* mu_fibre_2_nowash_20_vs_temperature = new TGraphErrors(npts_nowash_20, temperature_nowash_20, mu_fibre_2_nowash_20, temperature_err_nowash_20, mu_err_fibre_2_nowash_20);

  TGraphErrors* charge_PM_nowash_00_vs_temperature = new TGraphErrors(npts_nowash_00, temperature_nowash_00, charge_PM_nowash_00, temperature_err_nowash_00, charge_err_PM_nowash_00);
  TGraphErrors* charge_PM_nowash_10_vs_temperature = new TGraphErrors(npts_nowash_10, temperature_nowash_10, charge_PM_nowash_10, temperature_err_nowash_10, charge_err_PM_nowash_10);
  TGraphErrors* charge_PM_nowash_15_vs_temperature = new TGraphErrors(npts_nowash_15, temperature_nowash_15, charge_PM_nowash_15, temperature_err_nowash_15, charge_err_PM_nowash_15);
  TGraphErrors* charge_PM_nowash_20_vs_temperature = new TGraphErrors(npts_nowash_20, temperature_nowash_20, charge_PM_nowash_20, temperature_err_nowash_20, charge_err_PM_nowash_20);


  new TCanvas;
  mu_fibre_0_nowash_00_vs_temperature->GetXaxis()->SetLimits(0,35);
  mu_fibre_0_nowash_00_vs_temperature->GetXaxis()->SetRangeUser(0,35);
  mu_fibre_0_nowash_00_vs_temperature->GetYaxis()->SetRangeUser(0,10);
  mu_fibre_0_nowash_00_vs_temperature->GetYaxis()->SetTitle("#mu moyen");
  mu_fibre_0_nowash_00_vs_temperature->GetXaxis()->SetTitle("Temperature nowash (^{o}C)");
  mu_fibre_0_nowash_00_vs_temperature->SetTitle("Fibre 0");
  mu_fibre_0_nowash_00_vs_temperature->SetLineColor(kRed);
  mu_fibre_0_nowash_00_vs_temperature->SetMarkerColor(kRed);
  mu_fibre_0_nowash_00_vs_temperature->SetMarkerStyle(7);
  mu_fibre_0_nowash_00_vs_temperature->SetMarkerSize(3);
  mu_fibre_0_nowash_00_vs_temperature->Draw("ap");
  mu_fibre_0_nowash_10_vs_temperature->SetLineColor(kMagenta);
  mu_fibre_0_nowash_10_vs_temperature->SetMarkerColor(kMagenta);
  mu_fibre_0_nowash_10_vs_temperature->SetMarkerStyle(7);
  mu_fibre_0_nowash_10_vs_temperature->SetMarkerSize(3);
  mu_fibre_0_nowash_10_vs_temperature->Draw("psame");
  mu_fibre_0_nowash_15_vs_temperature->SetLineColor(kViolet+2);
  mu_fibre_0_nowash_15_vs_temperature->SetMarkerColor(kViolet+2);
  mu_fibre_0_nowash_15_vs_temperature->SetMarkerStyle(7);
  mu_fibre_0_nowash_15_vs_temperature->SetMarkerSize(3);
  mu_fibre_0_nowash_15_vs_temperature->Draw("psame");
  mu_fibre_0_nowash_20_vs_temperature->SetLineColor(kBlue);
  mu_fibre_0_nowash_20_vs_temperature->SetMarkerColor(kBlue);
  mu_fibre_0_nowash_20_vs_temperature->SetMarkerStyle(7);
  mu_fibre_0_nowash_20_vs_temperature->SetMarkerSize(3);
  mu_fibre_0_nowash_20_vs_temperature->Draw("psame");
  TLegend* l_fibre_0_vs_temperature = new TLegend(0.7,0.7,0.9,0.9);
  l_fibre_0_vs_temperature->AddEntry(mu_fibre_0_nowash_00_vs_temperature, "nowash 0%");
  l_fibre_0_vs_temperature->AddEntry(mu_fibre_0_nowash_10_vs_temperature, "nowash 10%");
  l_fibre_0_vs_temperature->AddEntry(mu_fibre_0_nowash_15_vs_temperature, "nowash 15%");
  l_fibre_0_vs_temperature->AddEntry(mu_fibre_0_nowash_20_vs_temperature, "nowash 20%");
  l_fibre_0_vs_temperature->Draw("lsame");

  new TCanvas;
  mu_fibre_1_nowash_00_vs_temperature->GetXaxis()->SetLimits(0,35);
  mu_fibre_1_nowash_00_vs_temperature->GetXaxis()->SetRangeUser(0,35);
  mu_fibre_1_nowash_00_vs_temperature->GetYaxis()->SetRangeUser(0,10);
  mu_fibre_1_nowash_00_vs_temperature->GetYaxis()->SetTitle("#mu moyen");
  mu_fibre_1_nowash_00_vs_temperature->GetXaxis()->SetTitle("Temperature nowash (^{o}C)");
  mu_fibre_1_nowash_00_vs_temperature->SetTitle("Fibre 1");
  mu_fibre_1_nowash_00_vs_temperature->SetLineColor(kRed);
  mu_fibre_1_nowash_00_vs_temperature->SetMarkerColor(kRed);
  mu_fibre_1_nowash_00_vs_temperature->SetMarkerStyle(7);
  mu_fibre_1_nowash_00_vs_temperature->SetMarkerSize(3);
  mu_fibre_1_nowash_00_vs_temperature->Draw("ap");
  mu_fibre_1_nowash_10_vs_temperature->SetLineColor(kMagenta);
  mu_fibre_1_nowash_10_vs_temperature->SetMarkerColor(kMagenta);
  mu_fibre_1_nowash_10_vs_temperature->SetMarkerStyle(7);
  mu_fibre_1_nowash_10_vs_temperature->SetMarkerSize(3);
  mu_fibre_1_nowash_10_vs_temperature->Draw("psame");
  mu_fibre_1_nowash_15_vs_temperature->SetLineColor(kViolet+2);
  mu_fibre_1_nowash_15_vs_temperature->SetMarkerColor(kViolet+2);
  mu_fibre_1_nowash_15_vs_temperature->SetMarkerStyle(7);
  mu_fibre_1_nowash_15_vs_temperature->SetMarkerSize(3);
  mu_fibre_1_nowash_15_vs_temperature->Draw("psame");
  mu_fibre_1_nowash_20_vs_temperature->SetLineColor(kBlue);
  mu_fibre_1_nowash_20_vs_temperature->SetMarkerColor(kBlue);
  mu_fibre_1_nowash_20_vs_temperature->SetMarkerStyle(7);
  mu_fibre_1_nowash_20_vs_temperature->SetMarkerSize(3);
  mu_fibre_1_nowash_20_vs_temperature->Draw("psame");
  TLegend* l_fibre_1_vs_temperature = new TLegend(0.7,0.7,0.9,0.9);
  l_fibre_1_vs_temperature->AddEntry(mu_fibre_1_nowash_00_vs_temperature, "nowash 0%");
  l_fibre_1_vs_temperature->AddEntry(mu_fibre_1_nowash_10_vs_temperature, "nowash 10%");
  l_fibre_1_vs_temperature->AddEntry(mu_fibre_1_nowash_15_vs_temperature, "nowash 15%");
  l_fibre_1_vs_temperature->AddEntry(mu_fibre_1_nowash_20_vs_temperature, "nowash 20%");
  l_fibre_1_vs_temperature->Draw("lsame");

  new TCanvas;
  mu_fibre_2_nowash_00_vs_temperature->GetXaxis()->SetLimits(0,35);
  mu_fibre_2_nowash_00_vs_temperature->GetXaxis()->SetRangeUser(0,35);
  mu_fibre_2_nowash_00_vs_temperature->GetYaxis()->SetRangeUser(0,10);
  mu_fibre_2_nowash_00_vs_temperature->GetYaxis()->SetTitle("#mu moyen");
  mu_fibre_2_nowash_00_vs_temperature->GetXaxis()->SetTitle("Temperature nowash (^{o}C)");
  mu_fibre_2_nowash_00_vs_temperature->SetTitle("Fibre 2");
  mu_fibre_2_nowash_00_vs_temperature->SetLineColor(kRed);
  mu_fibre_2_nowash_00_vs_temperature->SetMarkerColor(kRed);
  mu_fibre_2_nowash_00_vs_temperature->SetMarkerStyle(7);
  mu_fibre_2_nowash_00_vs_temperature->SetMarkerSize(3);
  mu_fibre_2_nowash_00_vs_temperature->Draw("ap");
  mu_fibre_2_nowash_10_vs_temperature->SetLineColor(kMagenta);
  mu_fibre_2_nowash_10_vs_temperature->SetMarkerColor(kMagenta);
  mu_fibre_2_nowash_10_vs_temperature->SetMarkerStyle(7);
  mu_fibre_2_nowash_10_vs_temperature->SetMarkerSize(3);
  mu_fibre_2_nowash_10_vs_temperature->Draw("psame");
  mu_fibre_2_nowash_15_vs_temperature->SetLineColor(kViolet+2);
  mu_fibre_2_nowash_15_vs_temperature->SetMarkerColor(kViolet+2);
  mu_fibre_2_nowash_15_vs_temperature->SetMarkerStyle(7);
  mu_fibre_2_nowash_15_vs_temperature->SetMarkerSize(3);
  mu_fibre_2_nowash_15_vs_temperature->Draw("psame");
  mu_fibre_2_nowash_20_vs_temperature->SetLineColor(kBlue);
  mu_fibre_2_nowash_20_vs_temperature->SetMarkerColor(kBlue);
  mu_fibre_2_nowash_20_vs_temperature->SetMarkerStyle(7);
  mu_fibre_2_nowash_20_vs_temperature->SetMarkerSize(3);
  mu_fibre_2_nowash_20_vs_temperature->Draw("psame");
  TLegend* l_fibre_2_vs_temperature = new TLegend(0.7,0.7,0.9,0.9);
  l_fibre_2_vs_temperature->AddEntry(mu_fibre_2_nowash_00_vs_temperature, "nowash 0%");
  l_fibre_2_vs_temperature->AddEntry(mu_fibre_2_nowash_10_vs_temperature, "nowash 10%");
  l_fibre_2_vs_temperature->AddEntry(mu_fibre_2_nowash_15_vs_temperature, "nowash 15%");
  l_fibre_2_vs_temperature->AddEntry(mu_fibre_2_nowash_20_vs_temperature, "nowash 20%");
  l_fibre_2_vs_temperature->Draw("lsame");

  new TCanvas;
  charge_PM_nowash_00_vs_temperature->GetXaxis()->SetLimits(0,35);
  charge_PM_nowash_00_vs_temperature->GetXaxis()->SetRangeUser(0,35);
  charge_PM_nowash_00_vs_temperature->GetYaxis()->SetRangeUser(0,10);
  charge_PM_nowash_00_vs_temperature->GetYaxis()->SetTitle("charge PM (x10^{4})");
  charge_PM_nowash_00_vs_temperature->GetXaxis()->SetTitle("Temperature nowash (^{o}C)");
  charge_PM_nowash_00_vs_temperature->SetTitle("PM");
  charge_PM_nowash_00_vs_temperature->SetLineColor(kRed);
  charge_PM_nowash_00_vs_temperature->SetMarkerColor(kRed);
  charge_PM_nowash_00_vs_temperature->SetMarkerStyle(7);
  charge_PM_nowash_00_vs_temperature->SetMarkerSize(3);
  charge_PM_nowash_00_vs_temperature->Draw("ap");
  charge_PM_nowash_10_vs_temperature->SetLineColor(kMagenta);
  charge_PM_nowash_10_vs_temperature->SetMarkerColor(kMagenta);
  charge_PM_nowash_10_vs_temperature->SetMarkerStyle(7);
  charge_PM_nowash_10_vs_temperature->SetMarkerSize(3);
  charge_PM_nowash_10_vs_temperature->Draw("psame");
  charge_PM_nowash_15_vs_temperature->SetLineColor(kViolet+2);
  charge_PM_nowash_15_vs_temperature->SetMarkerColor(kViolet+2);
  charge_PM_nowash_15_vs_temperature->SetMarkerStyle(7);
  charge_PM_nowash_15_vs_temperature->SetMarkerSize(3);
  charge_PM_nowash_15_vs_temperature->Draw("psame");
  charge_PM_nowash_20_vs_temperature->SetLineColor(kBlue);
  charge_PM_nowash_20_vs_temperature->SetMarkerColor(kBlue);
  charge_PM_nowash_20_vs_temperature->SetMarkerStyle(7);
  charge_PM_nowash_20_vs_temperature->SetMarkerSize(3);
  charge_PM_nowash_20_vs_temperature->Draw("psame");
  TLegend* l_PM_vs_temperature = new TLegend(0.7,0.7,0.9,0.9);
  l_PM_vs_temperature->AddEntry(charge_PM_nowash_00_vs_temperature, "nowash 0%");
  l_PM_vs_temperature->AddEntry(charge_PM_nowash_10_vs_temperature, "nowash 10%");
  l_PM_vs_temperature->AddEntry(charge_PM_nowash_15_vs_temperature, "nowash 15%");
  l_PM_vs_temperature->AddEntry(charge_PM_nowash_20_vs_temperature, "nowash 20%");
  l_PM_vs_temperature->Draw("lsame");

}

//
//	charge_vs_time routine (2018-08-13 Axel)
//
void charge_vs_time(const char* filename, const int E_cut, const int npe01, const int npe02, const int npe03, const int time_step)
{
  TFile *_file0 = new TFile(filename,"READ");
    TTree* reco_tree = (TTree*)_file0->Get("reco_tree");
    TTree* event_tree = (TTree*)_file0->Get("event_tree");

  ULong64_t TDC;
  float charge01, charge02, charge03, charge04;

  event_tree->SetBranchAddress("tdc",&TDC);
    reco_tree->SetBranchAddress("ch01_charge", &charge01);
    reco_tree->SetBranchAddress("ch02_charge", &charge02);
    reco_tree->SetBranchAddress("ch03_charge", &charge03);
    reco_tree->SetBranchAddress("ch04_charge", &charge04);

  std::vector<double> time_pm;
  std::vector<double> time_diff_pm;
  std::vector<double> time_sum_pm;
  std::vector<double> charge_pm;
  std::vector<double> charge_pm_err;
  std::vector<double> charge_sipm0;
  std::vector<double> charge_sipm0_err;
  std::vector<double> charge_sipm1;
  std::vector<double> charge_sipm1_err;
  std::vector<double> charge_sipm2;
  std::vector<double> charge_sipm2_err;

  time_pm         .clear();
  time_diff_pm    .clear();
  time_sum_pm     .clear();
  charge_pm       .clear();
  charge_pm_err   .clear();
  charge_sipm0    .clear();
  charge_sipm0_err.clear();
  charge_sipm1    .clear();
  charge_sipm1_err.clear();
  charge_sipm2    .clear();
  charge_sipm2_err.clear();

  double offset_PM = 0;

  for (int entry=0; entry<event_tree->GetEntries(); entry++)
    {
      event_tree->GetEntry(entry);
      if (entry == 0)
        time_pm.push_back(TDC);
      else
        {
          if (TDC<time_pm.at(entry-1)-offset_PM)
            offset_PM = time_pm.at(entry-1);
          time_pm.push_back(TDC+offset_PM);       //vect=TDC+offset != TDC.
        }
    }

  time_diff_pm.push_back(0);
  for (unsigned int entry=0; entry<time_pm.size()-1; entry++)
    time_diff_pm.push_back((time_pm.at(entry+1)-time_pm.at(entry))/200E6);

  time_sum_pm.push_back(0);
  for (unsigned int entry=0; entry<time_diff_pm.size()-1; entry++)
    time_sum_pm.push_back(time_sum_pm.at(entry)+time_diff_pm.at(entry+1));

  // TH1F* evolution_tdc = new TH1F("tdc","tdc",time_sum_pm.size(),0, time_sum_pm.at(time_sum_pm.size()-1));
  // for (unsigned int entry=0; entry<time_sum_pm.size()-1; entry++)
  //   evolution_tdc->SetBinContent(entry+1, time_sum_pm.at(entry));
  // evolution_tdc->Draw();

  vector<double> the_time;
  the_time.clear();

  for (int entry=0; entry<reco_tree->GetEntries(); entry++)
    {
      reco_tree->GetEntry(entry);
      if (-charge04>1000)
        {
          the_time.push_back(time_sum_pm.at(entry));
          charge_pm.push_back(-charge04);
          charge_pm_err.push_back(0);
          charge_sipm0.push_back(charge01);
          charge_sipm0_err.push_back(0);
          charge_sipm1.push_back(charge02);
          charge_sipm1_err.push_back(0);
          charge_sipm2.push_back(charge03);
          charge_sipm2_err.push_back(0);
            }
    }

  new TCanvas;
  TGraphErrors* charge_PM_vs_time = new TGraphErrors(the_time.size(), &(the_time[0]), &(charge_pm[0]), 0, &(charge_pm_err[0]));
  charge_PM_vs_time->SetTitle("PM");
  charge_PM_vs_time->GetXaxis()->SetTitle("Temps (s)");
  charge_PM_vs_time->GetYaxis()->SetTitle("Charge (pC)");
  charge_PM_vs_time->Draw("ap");
  new TCanvas;
  TGraphErrors* charge_SiPM0_vs_time = new TGraphErrors(the_time.size(), &(the_time[0]), &(charge_sipm0[0]), 0, &(charge_sipm0_err[0]));
  charge_SiPM0_vs_time->SetTitle("SiPM 0");
  charge_SiPM0_vs_time->GetXaxis()->SetTitle("Temps (s)");
  charge_SiPM0_vs_time->GetYaxis()->SetTitle("Charge (pC)");
  charge_SiPM0_vs_time->Draw("ap");
  new TCanvas;
  TGraphErrors* charge_SiPM1_vs_time = new TGraphErrors(the_time.size(), &(the_time[0]), &(charge_sipm1[0]), 0, &(charge_sipm1_err[0]));
  charge_SiPM1_vs_time->SetTitle("SiPM 1");
  charge_SiPM1_vs_time->GetXaxis()->SetTitle("Temps (s)");
  charge_SiPM1_vs_time->GetYaxis()->SetTitle("Charge (pC)");
  charge_SiPM1_vs_time->Draw("ap");
  new TCanvas;
  TGraphErrors* charge_SiPM2_vs_time = new TGraphErrors(the_time.size(), &(the_time[0]), &(charge_sipm2[0]), 0, &(charge_sipm2_err[0]));
  charge_SiPM2_vs_time->SetTitle("SiPM 2");
  charge_SiPM2_vs_time->GetXaxis()->SetTitle("Temps (s)");
  charge_SiPM2_vs_time->GetYaxis()->SetTitle("Charge (pC)");
  charge_SiPM2_vs_time->Draw("ap");

  vector<double> amplitude_sipm0;
  vector<double> amplitude_sipm0_err;
  vector<double> amplitude_sipm1;
  vector<double> amplitude_sipm1_err;
  vector<double> amplitude_sipm2;
  vector<double> amplitude_sipm2_err;
  vector<double> amplitude_pm;
  vector<double> amplitude_pm_err;
  vector<double> the_time_step;
  vector<double> the_time_step_err;

  amplitude_sipm0    .clear();
  amplitude_sipm0_err.clear();
  amplitude_sipm1    .clear();
  amplitude_sipm1_err.clear();
  amplitude_sipm2    .clear();
  amplitude_sipm2_err.clear();
  amplitude_pm       .clear();
  amplitude_pm_err   .clear();
  the_time_step      .clear();
  the_time_step_err  .clear();

  int k = 1;
  int nbr_steps = 0;
  for (int i=0; i<the_time.at(the_time.size()-1); i+=time_step)
    {
      cout << "Step " << nbr_steps+1 << endl;
      the_time_step.push_back(k*time_step);
      the_time_step_err.push_back(0);
      TH1F *charge01_histo  = new TH1F("charge01_histo", "SiPM 0", 400, 0, 20000);
      TH1F *charge02_histo  = new TH1F("charge02_histo", "SiPM 1", 400, 0, 20000);
      TH1F *charge03_histo  = new TH1F("charge03_histo", "SiPM 2", 400, 0, 20000);
      TH1F *charge04_histo  = new TH1F("charge04_histo", "PM"    , 1000, 0, 20000);
      for (int entry=0; entry<the_time.size(); entry++)
        {
          if (the_time.at(entry) > time_step*(k-1) && the_time.at(entry) <= time_step*k)
            {
              charge01_histo->Fill(charge_sipm0.at(entry));
                    charge02_histo->Fill(charge_sipm1.at(entry));
                    charge03_histo->Fill(charge_sipm2.at(entry));
              charge04_histo->Fill(charge_pm.at(entry));
            }
        }
      TF1* myfit01 = fit(npe01);
      TCanvas *c1 = new TCanvas;
      c1->SetLogy();
      charge01_histo->Draw();
      charge01_histo->Fit("fitf","q");
      myfit01->Draw("same");
      vector<double>* amplitude01 = calcul_amplitude(npe01, myfit01);
      amplitude_sipm0.push_back(amplitude01->at(0));
      amplitude_sipm0_err.push_back(amplitude01->at(1));
      cout << "µ position 0 = " << amplitude01->at(0) << " +- " << amplitude01->at(1) << " (chi2 = " << myfit01->GetChisquare()/myfit01->GetNDF() << ")" << endl;

      TF1* myfit02     = fit(npe02);
      TCanvas *c2     = new TCanvas;
      c2->SetLogy();
      charge02_histo->Draw();
      charge02_histo->Fit("fitf","q");
      myfit02->Draw("same");
      vector<double>* amplitude02 = calcul_amplitude(npe02, myfit02);
      amplitude_sipm1.push_back(amplitude02->at(0));
      amplitude_sipm1_err.push_back(amplitude02->at(1));
      cout << "µ position 1 = " << amplitude02->at(0) << " +- " << amplitude02->at(1) << " (chi2 = " << myfit02->GetChisquare()/myfit02->GetNDF() << ")" << endl;

      TF1* myfit03     = fit(npe03);
      TCanvas *c3     = new TCanvas;
      c3->SetLogy();
      charge03_histo->Draw();
      charge03_histo->Fit("fitf","q");
      myfit03->Draw("same");
      vector<double>* amplitude03 = calcul_amplitude(npe03, myfit03);
      amplitude_sipm2.push_back(amplitude03->at(0));
      amplitude_sipm2_err.push_back(amplitude03->at(1));
      cout << "µ position 2 = " << amplitude03->at(0) << " +- " << amplitude03->at(1) << " (chi2 = " << myfit03->GetChisquare()/myfit03->GetNDF() << ")" << endl;

      TCanvas *c4 = new TCanvas;
      c4->SetLogy();
      charge04_histo->Draw();
      TF1* fit_pm = new TF1("fit_pm","gaus",2000,20000);
      charge04_histo->Fit("fit_pm","RQ");
      for (int i=0; i<3; i++)
        {
            double mean_gaussian  = fit_pm->GetParameter(1);
            double sigma_gaussian = fit_pm->GetParameter(2);
            fit_pm->SetRange(mean_gaussian-0.8*sigma_gaussian, mean_gaussian+2.5*sigma_gaussian );
            charge04_histo->Fit("fit_pm","RQ");
        }
      fit_pm->Draw("same");
      amplitude_pm.push_back(fit_pm->GetParameter(1));
      amplitude_pm_err.push_back(fit_pm->GetParError(1));

      k += 1;
      nbr_steps += 1;
      c1->Close();
      c2->Close();
      c3->Close();
      c4->Close();
      charge01_histo->Delete();
      charge02_histo->Delete();
      charge03_histo->Delete();
      charge04_histo->Delete();
    }

  new TCanvas;
  TGraphErrors* charge_PM_vs_time_step = new TGraphErrors(nbr_steps, &(the_time_step[0]), &(amplitude_pm[0]), &(the_time_step_err[0]), &(amplitude_pm_err[0]));
  charge_PM_vs_time_step->SetTitle("PM");
  charge_PM_vs_time_step->GetXaxis()->SetTitle("Temps (s)");
  charge_PM_vs_time_step->GetYaxis()->SetTitle("Charge (pC)");
  charge_PM_vs_time_step->Draw("ap");
  new TCanvas;
  TGraphErrors* charge_SiPM0_vs_time_step = new TGraphErrors(nbr_steps, &(the_time_step[0]), &(amplitude_sipm0[0]), &(the_time_step_err[0]), &(amplitude_sipm0_err[0]));
  charge_SiPM0_vs_time_step->SetTitle("SiPM 0");
  charge_SiPM0_vs_time_step->GetXaxis()->SetTitle("Temps (s)");
  charge_SiPM0_vs_time_step->GetYaxis()->SetTitle("mu moyen");
  charge_SiPM0_vs_time_step->Draw("ap");
  new TCanvas;
  TGraphErrors* charge_SiPM1_vs_time_step = new TGraphErrors(nbr_steps, &(the_time_step[0]), &(amplitude_sipm1[0]), &(the_time_step_err[0]), &(amplitude_sipm1_err[0]));
  charge_SiPM1_vs_time_step->SetTitle("SiPM 1");
  charge_SiPM1_vs_time_step->GetXaxis()->SetTitle("Temps (s)");
  charge_SiPM1_vs_time_step->GetYaxis()->SetTitle("mu moyen");
  charge_SiPM1_vs_time_step->Draw("ap");
  new TCanvas;
  TGraphErrors* charge_SiPM2_vs_time_step = new TGraphErrors(nbr_steps, &(the_time_step[0]), &(amplitude_sipm2[0]), &(the_time_step_err[0]), &(amplitude_sipm2_err[0]));
  charge_SiPM2_vs_time_step->SetTitle("SiPM 2");
  charge_SiPM2_vs_time_step->GetXaxis()->SetTitle("Temps (s)");
  charge_SiPM2_vs_time_step->GetYaxis()->SetTitle("mu moyen");
  charge_SiPM2_vs_time_step->Draw("ap");
}
