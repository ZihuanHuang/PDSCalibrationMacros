// fitspe_sep.C
// Single-Gaussian SPE fit for September 2026 run 20530.
// The run uses the mixed 12.75/14.00 ADC HV settings in
// sep_voltage_settings.h. Input and all outputs live in ~/work/output/Sep/.

#include <set>
#include <vector>
#include <algorithm>
#include <fstream>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "sep_voltage_settings.h"

static const int RUN=20530;
static const double SEARCH_LO=13.0;
static const double SEARCH_HI=40.0;
static const double VALLEY_LO=5.0,VALLEY_HI=15.0;
static const double FIT_LO_REL=-4.0;
static const double FIT_HI_REL=+5.0;
static const int REBIN_FACTOR=2;
// Applied after Rebin(). 1 reproduces the original rule:
// fitLo cannot be below the low edge of bin (valleyBin + 1), i.e. the
// right-hand edge of the rebinned valley bin. Set to 0, 2, ... as needed.
static const int VALLEY_BIN_SHIFT=2;
static const double PED_GUARD=7.0;
static const double MAXERR_FRAC=0.5;

static std::set<int> readChannelsSep(TDirectory* d) {
  std::set<int> chans;
  if (!d) return chans;
  const TString prefix="amp_opchannel_";
  TIter next(d->GetListOfKeys()); TKey* key;
  while ((key=(TKey*)next())) {
    TString n=key->GetName();
    if (!n.BeginsWith(prefix)) continue;
    TString num=n(prefix.Length(),n.Length()-prefix.Length());
    if (num.IsDigit()) chans.insert(num.Atoi());
  }
  return chans;
}

// Valley = lowest local minimum in 5--15 ADC before the candidate SPE peak.
static int findValleyBinSep(TH1* h,int speBin) {
  int b1=h->GetXaxis()->FindBin(VALLEY_LO+1e-3);
  int b2=std::min(h->GetXaxis()->FindBin(VALLEY_HI-1e-3),speBin-1);
  if (b2-b1<2) return -1;
  int vbin=-1; double vmin=1e300;
  for (int b=b1+1;b<=b2-1;++b) {
    double c=h->GetBinContent(b);
    if (c<h->GetBinContent(b-1) && c<h->GetBinContent(b+1) && c<vmin) {
      vmin=c; vbin=b;
    }
  }
  return vbin;
}

void fitspe_sep() {
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.052,"XYZ");
  gStyle->SetLabelSize(0.044,"XYZ");
  gStyle->SetOptTitle(1);

  TString base=gSystem->ExpandPathName("~/work/output/Sep/");
  gSystem->mkdir(base,kTRUE);
  TString outpdf=base+"fitspe_sep.pdf";
  std::ofstream fout((base+"fitspe_sep.txt").Data());
  std::ofstream fflag((base+"fitspe_sep_flagged.txt").Data());
  fout << "# channel voltage[V] mean[ADC] mean_err[ADC] status target[ADC] hv_limit\n";
  fflag << "# channel voltage[V] reason\n";

  TFile* f=TFile::Open(base+Form("pmt_gain%d.root",RUN));
  if (!f || f->IsZombie()) { printf("ERROR opening run %d\n",RUN); return; }
  TDirectory* d=(TDirectory*)f->Get("PMTGain");
  if (!d) { printf("ERROR: PMTGain missing in run %d\n",RUN); f->Close(); return; }

  std::set<int> channelSet=readChannelsSep(d);
  std::vector<int> chans(channelSet.begin(),channelSet.end());
  std::sort(chans.begin(),chans.end());
  printf("Fitting %zu channels from run %d.\n",chans.size(),RUN);

  TCanvas* c=new TCanvas("csep","",950,700);
  c->Print(outpdf+"[");
  int nfit=0,nflag=0,noff=0,nlimit=0;

  for (int ch:chans) {
    auto sit=SEP_SETTINGS.find(ch);
    bool haveSetting=(sit!=SEP_SETTINGS.end());
    SepSetting setting={-1,0.0,false,false};
    if (haveSetting) setting=sit->second;
    if (!haveSetting) {
      fflag << ch << "  -1  no voltage entry\n";
      nflag++;
    }
    if (setting.off) noff++;
    if (setting.atLimit) nlimit++;

    c->Clear();
    c->SetLeftMargin(0.13); c->SetRightMargin(0.04);
    c->SetBottomMargin(0.13); c->SetTopMargin(0.11);
    TH1* h=(TH1*)d->Get(Form("amp_opchannel_%d",ch));
    if (!h) continue;
    h->Rebin(REBIN_FACTOR);

    int b1=h->GetXaxis()->FindBin(SEARCH_LO);
    int b2=h->GetXaxis()->FindBin(SEARCH_HI);
    int pbin=b1; double pmax=-1;
    for (int b=b1;b<=b2;++b) {
      if (h->GetBinContent(b)>pmax) { pmax=h->GetBinContent(b); pbin=b; }
    }
    double peak0=h->GetBinCenter(pbin);
    double fitLo=peak0+FIT_LO_REL;
    int vbin=findValleyBinSep(h,pbin);
    if (vbin>=0) {
      int edgeBin=std::max(1,std::min(vbin+VALLEY_BIN_SHIFT,h->GetNbinsX()+1));
      fitLo=std::max(fitLo,h->GetXaxis()->GetBinLowEdge(edgeBin));
    }
    double fitHi=peak0+FIT_HI_REL;

    TF1* g=new TF1(Form("gsep_%d",ch),"gaus",fitLo,fitHi);
    g->SetParameters(pmax,peak0,2.0);
    TFitResultPtr r=h->Fit(g,"RQS");
    int status=(r.Get()!=nullptr) ? r->Status() : -1;
    double mean=g->GetParameter(1),merr=g->GetParError(1);

    h->SetLineColor(setting.target==14.0 ? kOrange+7 : kBlue+1);
    h->SetLineWidth(2);
    h->GetXaxis()->SetRangeUser(0,40);
    TString settingLabel;
    if (!haveSetting) settingLabel="setting missing";
    else if (setting.off) settingLabel="PMT off";
    else if (setting.atLimit) settingLabel=Form("%d V limit, target %.2f ADC",setting.voltage,setting.target);
    else settingLabel=Form("%d V, target %.2f ADC",setting.voltage,setting.target);
    h->SetTitle(Form("ch %d | %s;Amplitude [ADC];Count",ch,settingLabel.Data()));
    h->Draw("HIST");

    TF1* gfull=new TF1(Form("gsep_full_%d",ch),"gaus",0,40);
    gfull->SetParameters(g->GetParameter(0),g->GetParameter(1),g->GetParameter(2));
    gfull->SetLineColor(kBlack); gfull->SetLineWidth(2);
    gfull->SetLineStyle(2); gfull->Draw("SAME");
    g->SetLineColor(kBlack); g->SetLineWidth(2);
    g->SetLineStyle(1); g->Draw("SAME");

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.041);
    tx.DrawLatex(0.63,0.80,Form("#mu = %.2f #pm %.2f",mean,merr));
    //tx.DrawLatex(0.63,0.74,Form("fit: %.2f--%.2f ADC",fitLo,fitHi));

    const char* st;
    if (!haveSetting) st="no_setting";
    else if (setting.off) st="off";
    else if (status==0) st="ok";
    else st="fail";
    fout << ch << "  " << setting.voltage << "  " << mean << "  " << merr
         << "  " << st << "  " << setting.target << "  " << (setting.atLimit ? 1 : 0) << "\n";
    nfit++;

    if (haveSetting && !setting.off) {
      TString reason;
      if (status!=0) reason += "did not converge; ";
      if (mean<PED_GUARD) reason += "mean near pedestal; ";
      if (mean>0 && merr/mean>MAXERR_FRAC) reason += "large error; ";
      if (setting.atLimit) reason += "HV at hardware limit; ";
      if (reason.Length()>0) {
        fflag << ch << "  " << setting.voltage << "  " << reason << "\n";
        nflag++;
      }
    }
    c->Print(outpdf);
  }

  c->Print(outpdf+"]");
  fout.close(); fflag.close();
  printf("Wrote %d fit points to fitspe_sep.txt (%d flagged, %d off, %d at limit).\n",
         nfit,nflag,noff,nlimit);
  printf("Diagnostic plots: %s\n",outpdf.Data());
  f->Close();
  delete c;
}
