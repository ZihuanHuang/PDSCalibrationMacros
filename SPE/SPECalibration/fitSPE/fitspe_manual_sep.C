// fitspe_manual_sep.C
// Manual refit for September 2026 run 20530.
// MANUAL[channel] = {fit_lo, fit_hi, init_mean}; init_mean=0 selects the
// highest histogram bin inside the custom range. Run fitspe_sep.C first.

#include <map>
#include <set>
#include <array>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "sep_voltage_settings.h"

static const int RUN=20530;

// Add only channels that need a manual fit, for example:
//   {140,{10.0,18.0,0.0}},
static std::map<int,std::array<double,3>> MANUAL = {
    {41,{10.0,18.0,12.0}},
    {60,{11.0,18.0,12.0}},
    {61,{10.0,18.0,12.0}},
    {62,{10.0,18.0,12.0}},
    {114,{11.0,18.0,12.0}},
    {116,{10.0,18.0,12.0}},
    {140,{10.0,18.0,12.0}},
    {142,{12.0,20.0,14.0}},
    {146,{10.0,18.0,12.0}},
    {171,{12.0,20.0,14.0}},
    {194,{10.0,18.0,12.0}},
    {195,{10.0,18.0,12.0}},
    {216,{13.0,20.0,14.0}},
    {227,{9.0,18.0,12.0}},
    {242,{10.0,18.0,12.0}},
    {249,{13.0,21.0,14.0}},
    {270,{10.0,18.0,12.0}},
    {272,{10.0,18.0,12.0}},
    {296,{10.0,18.0,12.0}},
    {298,{10.0,18.0,12.0}},
    {299,{10.0,18.0,12.0}},
    {300,{11.0,19.0,12.0}},
    {301,{10.0,18.0,12.0}},
    {305,{9.0,18.0,12.0}},
};

void fitspe_manual_sep() {
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.052,"XYZ");
  gStyle->SetLabelSize(0.044,"XYZ");
  gStyle->SetOptTitle(1);

  if (MANUAL.empty()) {
    printf("MANUAL table is empty. Review fitspe_sep.pdf and add channels to refit.\n");
    return;
  }

  TString base=gSystem->ExpandPathName("~/work/output/Sep/");
  TString txt=base+"fitspe_sep.txt";
  TString mpdf=base+"fitspe_sep_manual.pdf";

  TFile* f=TFile::Open(base+Form("pmt_gain%d.root",RUN));
  if (!f || f->IsZombie()) { printf("ERROR opening run %d\n",RUN); return; }
  TDirectory* d=(TDirectory*)f->Get("PMTGain");
  if (!d) { printf("ERROR: PMTGain missing in run %d\n",RUN); f->Close(); return; }

  std::vector<std::string> header,datalines;
  {
    std::ifstream in(txt.Data());
    if (!in) {
      printf("ERROR: cannot open %s (run fitspe_sep.C first)\n",txt.Data());
      f->Close(); return;
    }
    std::string line;
    while (std::getline(in,line)) {
      if (line.empty()) continue;
      if (line[0]=='#') header.push_back(line);
      else datalines.push_back(line);
    }
  }

  std::map<int,std::string> replacement;
  TCanvas* c=new TCanvas("cmsep","",950,700);
  c->Print(mpdf+"[");

  for (auto& item:MANUAL) {
    int ch=item.first;
    double flo=item.second[0],fhi=item.second[1],m0=item.second[2];
    auto sit=SEP_SETTINGS.find(ch);
    if (sit==SEP_SETTINGS.end()) {
      printf("WARNING: no September setting for ch %d; skipped.\n",ch);
      continue;
    }
    const SepSetting& setting=sit->second;
    TH1* h=(TH1*)d->Get(Form("amp_opchannel_%d",ch));
    if (!h) { printf("WARNING: histogram missing for ch %d; skipped.\n",ch); continue; }
    h->Rebin(2);

    int b1=h->GetXaxis()->FindBin(flo),b2=h->GetXaxis()->FindBin(fhi);
    double a0=h->GetMaximum();
    if (m0<=0) {
      double pm=-1; int pb=b1;
      for (int b=b1;b<=b2;++b) {
        if (h->GetBinContent(b)>pm) { pm=h->GetBinContent(b); pb=b; }
      }
      m0=h->GetBinCenter(pb); a0=pm;
    }

    TF1* g=new TF1(Form("gmsep_%d",ch),"gaus",flo,fhi);
    g->SetParameters(a0,m0,2.0);
    TFitResultPtr r=h->Fit(g,"RQS");
    int status=(r.Get()!=nullptr) ? r->Status() : -1;
    double mean=g->GetParameter(1),merr=g->GetParError(1);

    c->Clear();
    c->SetLeftMargin(0.13); c->SetRightMargin(0.04);
    c->SetBottomMargin(0.13); c->SetTopMargin(0.11);
    h->SetLineColor(setting.target==14.0 ? kOrange+7 : kBlue+1);
    h->SetLineWidth(2); h->GetXaxis()->SetRangeUser(0,40);
    TString settingLabel;
    if (setting.off) settingLabel="PMT off";
    else if (setting.atLimit) settingLabel=Form("%d V limit, target %.2f ADC",setting.voltage,setting.target);
    else settingLabel=Form("%d V, target %.2f ADC",setting.voltage,setting.target);
    h->SetTitle(Form("ch %d | %s [manual];Amplitude [ADC];Count",ch,settingLabel.Data()));
    h->Draw("HIST");

    TF1* gfull=new TF1(Form("gmsep_full_%d",ch),"gaus",0,40);
    gfull->SetParameters(g->GetParameter(0),g->GetParameter(1),g->GetParameter(2));
    gfull->SetLineColor(kBlack); gfull->SetLineWidth(2);
    gfull->SetLineStyle(2); gfull->Draw("SAME");
    g->SetLineColor(kBlack); g->SetLineWidth(2);
    g->SetLineStyle(1); g->Draw("SAME");

    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.041);
    tx.DrawLatex(0.63,0.80,Form("#mu = %.2f #pm %.2f",mean,merr));
    //tx.DrawLatex(0.63,0.74,Form("fit: %.2f--%.2f ADC",flo,fhi));

    const char* st;
    if (setting.off) st="off";
    else st=(status==0) ? "ok_manual" : "fail_manual";
    replacement[ch]=Form("%d  %d  %g  %g  %s  %g  %d",
                         ch,setting.voltage,mean,merr,st,setting.target,
                         setting.atLimit ? 1 : 0);
    c->Print(mpdf);
  }
  c->Print(mpdf+"]");

  std::ofstream out(txt.Data());
  for (auto& line:header) out << line << "\n";
  std::set<int> written;
  for (auto& line:datalines) {
    std::istringstream ss(line); int ch=-1; ss>>ch;
    auto it=replacement.find(ch);
    if (it!=replacement.end()) {
      out << it->second << "\n";
      written.insert(ch);
    } else out << line << "\n";
  }
  for (auto& item:replacement) {
    if (!written.count(item.first)) out << item.second << "\n";
  }
  out.close();

  printf("Refit %zu channel(s). Patched %s; wrote %s\n",
         replacement.size(),txt.Data(),mpdf.Data());
  f->Close();
  delete c;
}
