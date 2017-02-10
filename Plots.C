#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include <TColor.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TSpectrum.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <TLatex.h>
#include <TRandom3.h>

#include "/usr/local/include/utilities.h"

using namespace std;

//! S     I    Eg     Ece
//! Ba   133   356 - 
//! Na    22   511 -  340
//! Cs   137   662 -  478
//! Co    60  1250 - 1038
//! Co    60  1172 -  962 
//! Na    22  1274 - 1064
//! Co    60  1330 - 1117 
//! AmBe      4400 - 4196
void drawText(const char */*text*/, double /*xp*/, double /*yp*/, int /*size*/, Color_t /*icol*/);
double getEdge(TH1 */*htemp*/, int /*iChMin*/, int /*iChMax*/);

double rmsFraction = 0.68;
void adjust_range(TH1D */*h*/,double &/*fitmin*/, double &/*fitmax*/);

const int kBAR    [] = {873 , 874 , 875 , 876 , 
			877 , 878 , 879 , 880 , 
			881 , 882 , 883 , 884 , 
			885 , 886 , 887 , 888 , 
			889 , 890};

const int kVoltage[] = {1900, 1850, 1800, 1850, 
			1850, 1800, 1850, 1900, 
			1800, 1900, 1900, 1800, 
			1800, 1900, 1800, 1800, 
			1800, 1800};

Color_t kCol      [] = {
  kBlack   , kRed      , kGreen  , kBlue   , 
  kOrange  , kMagenta  , kViolet , kTeal ,
  kYellow  , kAzure+7  , kPink-2 , kSpring-7 ,
  kCyan    , kGray     , kOrange-3 , kMagenta-9,
  kGreen-6 , kRed+2
};

const int nbars = sizeof(kBAR)/sizeof(int);

const double deltparam[nbars][3]={
  {0.101076,-0.294288,1.9017},
  {0.0968632,-0.578553,2.01485},
  {0.0929069,0.54269,2.09822},
  {0.0997423,0.230336,1.96226},
  {0.0959197,-0.819494,2.01206},
  {0.0974916,-1.45718,1.98352},
  {0.0963687,0.335922,1.99265},
  {0.0950786,-0.375066,2.07102},
  {0.101955,-0.404347,1.89937},
  {0.092086,-0.274056,2.07209},
  {0.0922208,0.89881,2.13879},
  {0.101445,-1.32635,1.89283},
  {0.0989683,-0.0663093,1.96981},
  {0.0935818,-0.165555,2.0486},
  {0.0974414,-0.204958,1.96907},
  {0.0969539,-1.27048,1.98353},
  {0.0965415,2.62262,1.98178},
  {0.0943684,0.402036,2.04718}
};

TF1 *fgaus[nbars];

const char *kSide[2] = {"Near","Far"};

int plotGain(const char *source="Co60")// Na22
{

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  std::cout <<" nbars :"<< nbars << std::endl;
  //  return 0;

  int iFit=0;

  std::string indir = "/DigiDataA/BARS/SF";
  std::string outdir = "/home/user/BarPerformance/SF";

  TFile *fin = 0;
  TFile *fout = 0;

  double minCh=0;
  double maxCh=4096;


if( strcmp(source,"Cs137") == 0 ){
    minCh=0;  maxCh=1800;
  }else if( strcmp(source,"Na22") == 0 || strcmp(source,"Co60") == 0){
    minCh=0;  maxCh=2500;
  }else{
    minCh=0;  maxCh=32000;
  }

 
  TH1D *hCh    [nbars][2];
  TH1D *hDeltaT[nbars], *hGMean[nbars];
  TH1D *hDeltaTCorr[nbars];
  
  double minx=0.;
  double maxx=32000.;
  int nbins=maxx-minx;


  for(int i=0; i<nbars; i++){
    if( strcmp(source,"NoSource") == 0)fin = new TFile(Form("%s%d/Cosmic_%s_SF%d_%dV_CS1.root",indir.c_str(),kBAR[i],source,kBAR[i],kVoltage[i]),"r");
    else fin = new TFile(Form("%s%d/Couples_%s_+50cm_SF%d_%dV_CS1.root",indir.c_str(),kBar[i],source,kBar[i],kVoltage[i]),"r");

    fout = new TFile(Form("%s%d_%s_%dV_PerformancePlot.root",outdir.c_str(),kBar[i],source,kVoltage[i]),"RECREATE");
      

    for(int j=0; j<2; j++){

      hCh[i][j] = new TH1D(Form("hCh_SF%d_%s",kBAR[i],kSide[j]),"",nbins,minx-0.5,maxx-0.5);
      hCh[i][j]->SetLineColor(j+1);
      hCh[i][j]->SetLineWidth(2);
      hCh[i][j]->SetLineStyle(j+1);

    }

    hGMean [i] = new TH1D(Form("hGMean_SF%d",kBAR[i]),"",nbins,minx-0.5,maxx-0.5);
    hGMean [i]->SetLineColor(kCol[i]);
    hGMean [i]->SetLineWidth(2);

    hDeltaT[i] = new TH1D(Form("hDeltaT_SF%d",kBAR[i]),"",120,-33.,33.);
    hDeltaT[i]->SetLineColor(kCol[i]);
    hDeltaT[i]->SetLineWidth(2);

     hDeltaTCorr[i] = new TH1D(Form("hDeltaTCorr_SF%d",kBAR[i]),"",65,-6.5,6.5);
    hDeltaTCorr[i]->SetLineColor(kCol[i]);
    hDeltaTCorr[i]->SetLineWidth(2);

    //! Gaussain shift correction for deltaT for all the bars (walk correction)
    fgaus[i]= new TF1(Form("fgaus_%d",i),"gaus",-20.,20.);
    fgaus[i]->SetNpx(10000);
    for(int j=0;j<3;j++){
      fgaus[i]->SetParameter(j,deltparam[i][j]);
    }

  UShort_t brch, qlong;
  ULong64_t time;
  TTree *tr= (TTree*)fin->Get("ftree");
  tr->SetBranchAddress("fBrCh",&brch);
  tr->SetBranchAddress("fQlong",&qlong);
  tr->SetBranchAddress("fTstamp",&time);
  
  int istat=0;
  UShort_t  qno[2]={0};
  UShort_t  qch[2]={0};
  ULong64_t qst[2]={0};
  
  int nEntries = tr->GetEntries();
  for(int iev=0; iev<nEntries; iev++){
  tr->GetEntry(iev);
  if( iev % 100000 == 0 )cout << " Processing event " << iev << endl;
  
    qno[istat]=brch;
    qch[istat]=qlong;
    qst[istat]=time;
    istat++;
      
    double delt=-999;
    double gmean=-999;
    double deltcorr=-999;
    int iBar = -1;
    
    if( istat == 2 ){
      if( qno[0] == qno[1] ){
	std::cout <<" Something wrong pile up event " << std::endl;
	goto skipEvent;
      }
      
      //      if( iBar < 0 ) goto skipEvent;
      
      hCh[iBar][0]->Fill(1.0*qch[0]);
      hCh[iBar][1]->Fill(1.0*qch[1]);
      delt = (1.0*qst[0] - 1.0*qst[1])/1000.;
      deltcorr = (delt - fgaus[i]->GetParameter(1))/fgaus[i]->GetParameter(2);
      //cout <<" delt : " << delt << " deltcorr : " << deltcorr << "  val "  << fgaus[i]->Eval(delt) << endl;
      hDeltaT[iBar]->Fill(delt);
      hDeltaTCorr[i]->Fill(deltcorr);
      
      gmean = sqrt(qch[0]*qch[1]);
      hGMean[iBar]->Fill(gmean);
      
    skipEvent:
      qno[0] = qno[1] = 0;
      qst[0] = qst[1] = 0;
      qch[0] = qch[1] = 0;
      istat=0;
    }
  }

  // for(int i=0; i< nbars;  i++){
  //   std::cout << kBAR[i]<< " " << source << "  # of entries in tree : " << nEntries << std::endl
  // 	      <<" \t Ch0 : " << hGMean[0]->GetEntries() << std::endl  
  // 	      <<" \t Ch1 : " << hGMean[1]->GetEntries() << std::endl;
  // }
  //return 0;

  gStyle->SetOptStat(0);
  TLegend *leg0 = new TLegend(0.74,0.78,0.90,0.98, "","brNDC");
  leg0->SetHeader("");
  leg0->SetTextSize(0.05);
  leg0->AddEntry(hCh[0][0],"Near","l");
  leg0->AddEntry(hCh[0][1],"Far","l");

  TLegend *leg1 = new TLegend(0.6937669,0.3734729,0.9742547,0.99127, "","brNDC");
  leg1->SetHeader("");
  leg1->SetTextSize(0.03);

  for(int i=0; i<nbars; i++){
    for(int j=0; j<2; j++){
      hCh[i][j]->Scale(1./hCh[i][j]->Integral());
    }
    hGMean[i]->Scale(1./hGMean[i]->Integral());
    hDeltaT[i]->Scale(1./hDeltaT[i]->Integral());
    leg1->AddEntry(hGMean[i],Form("%dcm (-%dV)",kBAR[i],kVoltage[i]),"l");
  }  
  
  double max=-999;
  TCanvas *c1 = new TCanvas("c1","Near Far",1200,800);
  c1->Divide(3,2,1.e-9,1.e-9); c1->cd(1);//I
  gPad->SetLeftMargin(0.20);
  gPad->SetBottomMargin(0.20);
  gPad->SetLogy();
  for(int i=0; i<nbars; i++){
    max = TMath::Max(hCh[i][0]->GetMaximum(), hCh[i][1]->GetMaximum());
    hCh[i][0]->SetMaximum(10*max);
    hCh[i][0]->SetMinimum(1e-06);
    MakeHist(hCh[i][0],"Channel No.","Normalized Counts");
    if( i==0 ){
      hCh[i][0]->Draw("hist");
      leg0->Draw();
    }
    else hCh[i][0]->Draw("histsame");
    MakeHist(hCh[i][1],"Channel No.","Normalized Counts");
    hCh[i][1]->Draw("histsame");

    hCh[i][0]->GetXaxis()->SetRangeUser(minCh,maxCh);
    hCh[i][1]->GetXaxis()->SetRangeUser(minCh,maxCh);
    gPad->Update();
  }
  
  //  TCanvas *c2 = new TCanvas("c2","DeltaT",740,600);
  c1->cd(2);
  gPad->SetLeftMargin(0.20);
  gPad->SetBottomMargin(0.20);
  for(int i=0; i<nbars; i++){
    hDeltaT[i]->Scale(1./hDeltaT[i]->Integral());
    MakeHist(hDeltaT[i],"#Delta T (in n.s.)","Normalized Counts");
    if( i==0 ){
      hDeltaT[i]->Draw("hist");
      leg1->Draw();
    }
    else hDeltaT[i]->Draw("histsame");
    hDeltaT[i]->GetYaxis()->SetRangeUser(0,0.15);
  }
  //Insert part of fitting delt plots to extract params if necessary

  //   TCanvas *c22 = new TCanvas("c22","CorrDeltaT",740,600);
  c1->cd(3)
  gPad->SetLeftMargin(0.20);
  gPad->SetBottomMargin(0.20);
  for(int i=0; i<nbars; i++){
    hDeltaTCorr[i]->Scale(1./hDeltaTCorr[i]->Integral());
    MakeHist(hDeltaTCorr[i],"n#sigma(#Delta T)","Normalized Counts");
    if( i==0 ){
      hDeltaTCorr[i]->Draw("hist");
      leg1->Draw();
    }
    else hDeltaTCorr[i]->Draw("histsame");
    hDeltaTCorr[i]->GetYaxis()->SetRangeUser(0,0.095);
  }
  
  //  TCanvas *c3 = new TCanvas("c3","GM",740,600);
  c1->cd(4);
  gPad->SetLeftMargin(0.20);
  gPad->SetBottomMargin(0.20);
  gPad->SetLogy();
  for(int i=0; i<nbars; i++){
    MakeHist(hGMean[i],"Geo. Mean","Normalized Counts");
    if( i==0 ){
      hGMean[i]->Draw("hist");
      leg1->Draw();
    }
    else hGMean[i]->Draw("histsame");
    hGMean[i]->GetYaxis()->SetRangeUser(1e-06,10*max);
    hGMean[i]->GetXaxis()->SetRangeUser(minCh,maxCh);
  }
  //  fout->mkdir("Histograms");
  //  fout->cd("Histograms");
  //  fout->cd();
  c1->SaveAs(Form("%s%d_%s_%dV_PerformancePlot.pdf",outdir.c_str(),kBar[i],source,kVoltage[i]),"RECREATE");
    fout->Write();
    fout->Close();


  return 0;
}
double getEdge(TH1 *htemp, int iChMin, int iChMax)
{
  double edgeCh  = -1;
  double edgeVal = -9999;
  for(int i=iChMin; i<=iChMax; i++){
    double con = htemp->GetBinContent(i);
    double cen = htemp->GetBinCenter(i);
  }
  return edgeCh;
}
void adjust_range(TH1D *h,double &fitmin, double &fitmax)
{
  double mean = h->GetMean();
  double evFac=0;
  double intg = h->Integral();
  int meanBin = h->FindBin(mean);
  int iBin=1;
  while( evFac < rmsFraction ){
    evFac += h->GetBinContent(meanBin+iBin);
    evFac += h->GetBinContent(meanBin-iBin);
    fitmin = h->GetBinCenter(meanBin-iBin);
    fitmax = h->GetBinCenter(meanBin+iBin);
    iBin++;
  }
}
void drawText(const char *text, double xp, double yp, int size, Color_t col)
{
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(col);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
// int GetBarNo(int brch)
// {
//   int iBar=-1;
//   switch(brch){
//     //!  Board 0
//   case 0:
//   case 1:    
//     iBar = 0;
//     break;
//   case 2:
//   case 3:    
//     iBar = 1;
//     break;
//   case 4:
//   case 5:    
//     iBar = 2;
//     break;
//   case 6:
//   case 7:    
//     iBar = 3;
//     break;
//   case 8:
//   case 9:    
//     iBar = 4;
//     break;
//   case 10:
//   case 11:    
//     iBar = 5;
//     break;
//   case 12:
//   case 13:    
//     iBar = 6;
//     break;
//   case 14:
//   case 15:    
//     iBar = 7;
//     break;
//   case 16:
//   case 17:    
//     iBar = 8;
//     break;
//   case 18:
//   case 19:    
//     iBar = 9;
//     break;
//   case 20:
//   case 21:    
//     iBar = 10;
//     break;
//   case 22:
//   case 23:    
//     iBar = 11;
//     break;
//   case 24:
//   case 25:    
//     iBar = 12;
//     break;
//   case 26:
//   case 27:    
//     iBar = 13;
//     break;
//   case 28:
//   case 29:    
//     iBar = 14;
//     break;
//   case 30:
//   case 31:    
//     iBar = 6;
//     break;
//   default:
//     iBar=-1;
//     break;
//   };

//   return iBar;
// }
