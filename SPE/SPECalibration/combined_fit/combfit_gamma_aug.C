// combfit_gamma_aug.C
//
// August 2026 combined fit with NO low-voltage pedestal reference.
// Each (channel, run) is fitted directly with
//
//   total(x) = shifted-gamma pedestal + SPE Gaussian
//
// Pedestal treatment:
//   - do NOT use a low-voltage reference pedestal;
//   - unlike the Gaussian/exponential pedestal tests, the gamma model fits the
//     WHOLE pedestal, including its rising side and peak;
//   - the left edge is chosen from the measured pedestal itself: walk left from
//     the pedestal peak until the data fall below PED_LEFT_FRAC of the peak,
//     with minimum/maximum left spans and a hard floor to avoid sparse baseline;
//   - a pedestal-only gamma prefit (left edge -> valley/right-side cap) is used
//     ONLY for initial values.  The final gamma parameters are free in the
//     simultaneous gamma + Gaussian fit.
//
// The shifted gamma is parameterised by N, k, theta and its peak position mu_ped:
//
//   x0 = mu_ped - (k-1)*theta
//   gamma(x) = N * t^(k-1) * exp(-t),  t=(x-x0)/theta,  x>x0.
//
// By default mu_ped is fixed to the measured pedestal peak.  This keeps the
// free k/theta shape fit from sliding the whole pedestal under the SPE; set
// FIX_PED_PEAK=false to allow a small +/- PED_PEAK_FLOAT motion.
//
// Inputs (same August scan as fitspe_aug.C):
//   ~/work/output/Aug_scan/pmt_gain20468.root   new nominal
//   ~/work/output/Aug_scan/pmt_gain20469.root   +10 V
//   ~/work/output/Aug_scan/pmt_gain20470.root   +20 V
//   ~/work/output/Aug_scan/pmt_gain20471.root   +30 V
//
// SPE seeds and the final single-Gaussian parameters are read from
// ~/work/output/Aug_scan/fitspe_aug.txt when available. The updated fitspe
// macros append amp, sigma, fit_lo, fit_hi and bin_width so the comparison PDF
// can reconstruct the final single-Gaussian fit without copying its algorithm.
//
// Outputs:
//   ~/work/output/combfit_aug/combfit_gamma_aug.txt
//   ~/work/output/combfit_aug/combfit_gamma_aug.pdf
//   ~/work/output/combfit_aug/combfit_gamma_aug_residual.pdf
//   ~/work/output/combfit_aug/combfit_gamma_aug_vs_single.pdf
//
// Usage:
//   root -l -b -q 'combfit_gamma_aug.C'

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TError.h>
#include <Math/MinimizerOptions.h>

#include <cstdio>
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <limits>

// ------------------------------------------------------------------ config --

static const char* FPAT  = "pmt_gain%d.root";
static const char* HDIR  = "PMTGain";
static const char* HNAME = "amp_opchannel_%d";

static const int NRUN = 4;
static const int RUNS[NRUN]   = {20468, 20469, 20470, 20471};
static const int OFFSET[NRUN] = {    0,    10,    20,    30};

// Histogram rebinning.  REBIN=1 keeps the original binning (default: no rebin).
// Set REBIN=2, 4, ... to merge that many neighbouring bins before fitting.
static const int REBIN = 1;

// pedestal peak finder: first prominent peak from the left
static const int    LOCAL_W   = 2;
static const double DROP_FRAC = 0.15;
static const double PED_LO    = 1.0;
static const double PED_HI    = 15.0;

// SPE seed/search region
static const double SPE_SEARCH_HI = 40.0;
static const double MU_GAP        = 2.0;   // SPE mu must remain >= ped_peak + MU_GAP
// Peak-to-valley (P/V) metric used for the CLEAR/BURIED label.
// P/V = observed SPE-peak height / observed valley height.
// The threshold is deliberately exposed: 1.10 is only a starting value and
// should be tuned to the criterion agreed with the supervisor.
static const double PV_RATIO_MIN        = 1.10;
static const double PV_EDGE_GAP         = 0.50;  // keep valley search away from both peaks [ADC]
static const double PV_SPE_HALF_WINDOW  = 1.00;  // search observed SPE maximum around fitted mu [ADC]

// Gamma pedestal: include the rising side as well as the peak/tail.
// Left edge is data-driven: move left from the pedestal peak until the histogram
// drops below PED_LEFT_FRAC of the pedestal height.  The span guards prevent a
// very narrow statistical crossing or a long sparse baseline from dominating.
static const double PED_LEFT_FRAC      = 0.02;  // 2% of pedestal-peak height
static const double PED_LEFT_MIN_SPAN  = 1.5;   // include at least this much left side [ADC]
static const double PED_LEFT_MAX_SPAN  = 4.0;   // never go farther left than this [ADC]
static const double PED_FLOOR          = 0.5;   // absolute lower x bound [ADC]

// Pedestal-only gamma INITIALISATION fit: full left side through the falling
// side, but stop before the SPE.  This fit supplies seeds only.
static const double PED_PREFIT_MIN_RIGHT = 1.5; // need at least this much right of peak [ADC]
static const double PED_PREFIT_MAX_RIGHT = 4.5; // right-side cap [ADC]
static const double PED_TO_SPE_GAP       = 1.0; // stop before SPE seed [ADC]

// Shifted-gamma shape.  The earlier low-voltage gamma study used k>1 and broad
// k/theta limits.  Here k and theta are fitted directly for every run/channel.
static const double GAMMA_K_SEED     = 3.0;
static const double GAMMA_K_MIN      = 1.05;
static const double GAMMA_K_MAX      = 60.0;
static const double GAMMA_THETA_MIN  = 0.05;
static const double GAMMA_THETA_MAX  = 15.0;
static const bool   FIX_PED_PEAK      = true;
static const double PED_PEAK_FLOAT    = 0.75; // only used if FIX_PED_PEAK=false

static const double SPE_SIGMA_SEED = 4.0;
static const double SPE_SIGMA_MIN  = 0.8;
static const double SPE_SIGMA_MAX  = 12.0;

// Final combined-fit right edge.  The previous version used
//   hi = speSeed + 3 * SPE_SIGMA_SEED = speSeed + 12 ADC,
// which included much more of the SPE right tail than needed here.
// 7 ADC roughly matches the upper-side span used in the August single-Gaussian
// fits (+5/+6/+6/+7 ADC for offsets 0/+10/+20/+30).
static const double SPE_RIGHT_SPAN = 7.0;   // hi ~= SPE seed + this many ADC
static const double MIN_COMB_SPAN  = 8.0;   // keep at least this much range from pedestal peak
static const double FIT_XMAX       = 50.0;
static const int    MIN_BINS       = 10;

// diagnostics / quality labels
// Residuals below are standardized residuals: (data-fit)/sqrt(data).
// residual_frac is the longest run of same-sign standardized residual bins
// (to the right of the pedestal peak), divided by the number of such bins.
static const double RESIDUALFRAC_BAD = 0.60;
static const double EDGE_EPS     = 0.05;

// ----------------------------------------------------------------- structs --

struct SeedPoint {
  double mean = std::numeric_limits<double>::quiet_NaN();
  double meanErr = std::numeric_limits<double>::quiet_NaN();

  // Optional shape parameters written by the updated fitspe_aug.C /
  // fitspe_manual_aug.C.  Keeping these in the single-fit output avoids
  // duplicating the single-Gaussian fitting algorithm here.
  double amp = std::numeric_limits<double>::quiet_NaN();
  double sigma = std::numeric_limits<double>::quiet_NaN();
  double fitLo = std::numeric_limits<double>::quiet_NaN();
  double fitHi = std::numeric_limits<double>::quiet_NaN();
  double binWidth = std::numeric_limits<double>::quiet_NaN();

  bool ok = false;
  bool off = false;
  bool hasShape = false;
  bool hasWindow = false;
};

struct RunFit {
  bool filled = false;
  bool usable = false;
  TH1* h = nullptr;

  int run = 0;
  int offset = 0;
  int fitStatus = -999;

  double pedpk = 0;
  double pedInitHi = 0;
  double lo = 0;
  double hi = 0;

  double pedN = 0;
  double pedK = 0;
  double pedTheta = 0;
  double pedMu = 0;
  double pedX0 = 0;

  double speA = 0;
  double speMu = 0;
  double speMuErr = 0;
  double speSig = 0;

  double c2 = 0;
  double residualFrac = 0;
  double sep = 0;

  double valleyX = std::numeric_limits<double>::quiet_NaN();
  double valleyY = std::numeric_limits<double>::quiet_NaN();
  double spePeakX = std::numeric_limits<double>::quiet_NaN();
  double spePeakY = std::numeric_limits<double>::quiet_NaN();
  double pvRatio = std::numeric_limits<double>::quiet_NaN();
  bool pvValid = false;

  std::vector<double> residualX;
  std::vector<double> residualY;

  TString verdict;
  TString quality;
};

// ----------------------------------------------------------------- model ----
// p0 = N_gamma, p1 = k, p2 = theta, p3 = gamma peak position mu_ped
// p4 = A_spe,   p5 = mu_spe, p6 = sigma_spe
// x0 is derived so that the gamma maximum is at mu_ped:
//   x0 = mu_ped - (k-1)*theta.

double GammaValue(double x, double N, double k, double theta, double muPed)
{
  if (!(N >= 0.0) || !(k > 1.0) || !(theta > 0.0)) return 0.0;
  const double x0 = muPed - (k - 1.0) * theta;
  const double t = (x - x0) / theta;
  if (t <= 0.0) return 0.0;
  // Evaluate in log space to avoid pow() overflow during minimizer excursions.
  const double logShape = (k - 1.0) * std::log(t) - t;
  if (logShape < -745.0) return 0.0;
  if (logShape > 700.0) return N * std::exp(700.0);
  return N * std::exp(logShape);
}

double GammaGauss(double* x, double* p)
{
  const double ped = GammaValue(x[0], p[0], p[1], p[2], p[3]);
  const double z = (x[0] - p[5]) / p[6];
  const double spe = p[4] * std::exp(-0.5 * z * z);
  return ped + spe;
}

double GammaComp(double* x, double* p)
{
  return GammaValue(x[0], p[0], p[1], p[2], p[3]);
}

double SpeComp(double* x, double* p)
{
  const double z = (x[0] - p[1]) / p[2];
  return p[0] * std::exp(-0.5 * z * z);
}

// ----------------------------------------------------------------- helpers --

TH1* GetH(TFile* f, int ch)
{
  if (!f) return nullptr;
  TDirectory* d = (TDirectory*)f->Get(HDIR);
  return d ? (TH1*)d->Get(Form(HNAME, ch)) : nullptr;
}

std::set<int> ReadChannels(TDirectory* d)
{
  std::set<int> chans;
  if (!d) return chans;

  const TString prefix = "amp_opchannel_";
  TIter next(d->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = (TKey*)next())) {
    TString n = key->GetName();
    if (!n.BeginsWith(prefix)) continue;
    TString num = n(prefix.Length(), n.Length() - prefix.Length());
    if (num.IsDigit()) chans.insert(num.Atoi());
  }
  return chans;
}

// First prominent peak from the left.  This is deliberately not the global
// maximum because the SPE peak can be taller than the pedestal.
double FirstPeak(TH1* h)
{
  const int b0 = h->FindBin(PED_LO + 1e-6);
  const int b1 = h->FindBin(PED_HI - 1e-6);

  for (int b = b0 + LOCAL_W; b <= b1 - LOCAL_W; ++b) {
    const double y = h->GetBinContent(b);
    if (y <= 0) continue;

    bool localMax = true;
    for (int j = -LOCAL_W; j <= LOCAL_W; ++j) {
      if (h->GetBinContent(b + j) > y) {
        localMax = false;
        break;
      }
    }
    if (!localMax) continue;

    double after = y;
    for (int j = 1; j <= 2 * LOCAL_W + 1 && b + j <= b1; ++j)
      after = std::min(after, h->GetBinContent(b + j));

    if (after < (1.0 - DROP_FRAC) * y)
      return h->GetBinCenter(b);
  }

  // fallback: largest bin in pedestal search region
  double best = -1.0;
  double xbest = 0.0;
  for (int b = b0; b <= b1; ++b) {
    const double y = h->GetBinContent(b);
    if (y > best) {
      best = y;
      xbest = h->GetBinCenter(b);
    }
  }
  return xbest;
}

// Search for a local minimum between the pedestal peak and the SPE seed.
// If several exist, keep the lowest one.  Returns -1 if none is usable.
double ValleyBetween(TH1* h, double pedpk, double speSeed)
{
  const double xlo = pedpk + PED_PREFIT_MIN_RIGHT;
  const double xhi = speSeed - PED_TO_SPE_GAP;
  if (xhi <= xlo) return -1.0;

  int b0 = h->FindBin(xlo + 1e-6);
  int b1 = h->FindBin(xhi - 1e-6);
  b0 = std::max(2, b0);
  b1 = std::min(h->GetNbinsX() - 1, b1);
  if (b1 <= b0) return -1.0;

  int bestBin = -1;
  double bestY = 1e300;
  for (int b = b0; b <= b1; ++b) {
    const double y = h->GetBinContent(b);
    if (y < 0) continue;
    if (y <= h->GetBinContent(b - 1) &&
        y <= h->GetBinContent(b + 1) && y < bestY) {
      bestY = y;
      bestBin = b;
    }
  }

  return (bestBin >= 0) ? h->GetBinCenter(bestBin) : -1.0;
}

// Fallback SPE seed: tallest bin to the right of the pedestal.
double SpectrumSpeSeed(TH1* h, double pedpk)
{
  const double slo = std::max(pedpk + MU_GAP, 8.0);
  int b0 = h->FindBin(slo + 1e-6);
  int b1 = h->FindBin(SPE_SEARCH_HI - 1e-6);

  double best = -1.0;
  int bestBin = b0;
  for (int b = b0; b <= b1; ++b) {
    const double y = h->GetBinContent(b);
    if (y > best) {
      best = y;
      bestBin = b;
    }
  }
  return h->GetBinCenter(bestBin);
}

// Data-driven left edge for the gamma pedestal.  This follows the same idea as
// the earlier gamma pedestal study, which defined the pedestal window by a
// fraction of the pedestal peak height, but only the LEFT edge is used here;
// the final fit must continue through the SPE.
double GammaLeftEdge(TH1* h, double pedpk)
{
  const int pb = h->FindBin(pedpk);
  const double H = h->GetBinContent(pb);
  if (H <= 0) return std::max(PED_FLOOR, pedpk - PED_LEFT_MIN_SPAN);

  const double thr = PED_LEFT_FRAC * H;
  int bl = pb;
  while (bl > 1) {
    const int bn = bl - 1;
    const double lowEdge = h->GetBinLowEdge(bn);
    if (lowEdge < PED_FLOOR) break;
    if (h->GetBinContent(bn) < thr) break;
    bl = bn;
  }

  const double fracLo = h->GetBinLowEdge(bl);
  const double farLo  = std::max(PED_FLOOR, pedpk - PED_LEFT_MAX_SPAN);
  const double nearLo = std::max(PED_FLOOR, pedpk - PED_LEFT_MIN_SPAN);

  // Do not go farther left than MAX_SPAN; do include at least MIN_SPAN.
  double lo = std::max(farLo, fracLo);
  lo = std::min(lo, nearLo);
  return lo;
}

struct GammaSeedFit {
  bool ok = false;
  double N = 0;
  double k = GAMMA_K_SEED;
  double theta = 1.0;
  double muPed = 0;
  double x0 = 0;
};

// Full-pedestal gamma prefit for INITIAL VALUES ONLY.  It starts on the rising
// side and extends to the valley/right-side cap, so k and theta see both sides
// of the pedestal instead of only its falling tail.
GammaSeedFit FullPedGammaSeed(TH1* h, double pedpk, double speSeed,
                              double lo, double& prefitHi)
{
  GammaSeedFit out;

  double hi = std::min(pedpk + PED_PREFIT_MAX_RIGHT,
                       speSeed - PED_TO_SPE_GAP);
  const double valley = ValleyBetween(h, pedpk, speSeed);
  if (valley > pedpk + PED_PREFIT_MIN_RIGHT)
    hi = std::min(hi, valley);
  if (hi < pedpk + PED_PREFIT_MIN_RIGHT)
    hi = pedpk + PED_PREFIT_MIN_RIGHT;
  prefitHi = hi;

  const double H = std::max(1.0, h->GetBinContent(h->FindBin(pedpk)));
  const double k0 = GAMMA_K_SEED;
  // Choose theta so the implied support start x0 initially sits slightly left
  // of the chosen fit edge.  This gives the prefit a visible rising side.
  const double x0seed = std::max(0.0, lo - 0.25);
  double theta0 = (pedpk - x0seed) / std::max(0.2, k0 - 1.0);
  theta0 = std::min(std::max(theta0, GAMMA_THETA_MIN * 2.0),
                    0.8 * GAMMA_THETA_MAX);
  const double tpk = k0 - 1.0;
  const double shapePk = std::pow(tpk, tpk) * std::exp(-tpk);
  const double N0 = H / std::max(1e-12, shapePk);

  TF1 fp("fgamma_seed", GammaComp, lo, hi, 4);
  fp.SetParameters(N0, k0, theta0, pedpk);
  fp.SetParLimits(0, 0.0, 1e6 * N0 + 1e3);
  fp.SetParLimits(1, GAMMA_K_MIN, GAMMA_K_MAX);
  fp.SetParLimits(2, GAMMA_THETA_MIN, GAMMA_THETA_MAX);
  fp.FixParameter(3, pedpk);
  h->Fit(&fp, "LQNR");

  const double N = fp.GetParameter(0);
  const double k = fp.GetParameter(1);
  const double th = fp.GetParameter(2);
  if (!std::isfinite(N) || !std::isfinite(k) || !std::isfinite(th) ||
      N <= 0 || k <= 1.0 || th <= 0) {
    out.N = N0; out.k = k0; out.theta = theta0; out.muPed = pedpk;
  } else {
    out.ok = true; out.N = N; out.k = k; out.theta = th; out.muPed = pedpk;
  }
  out.x0 = out.muPed - (out.k - 1.0) * out.theta;
  return out;
}

// Observed peak-to-valley ratio for judging whether the SPE is buried.
// The valley is the lowest DATA bin between the pedestal and fitted SPE peak,
// excluding a small gap near each peak.  The SPE peak height is the largest DATA
// bin within +/- PV_SPE_HALF_WINDOW of the fitted SPE mean.
// This deliberately uses the histogram values rather than the fitted component
// heights, so the metric answers the visual question: is there a real dip in
// the measured spectrum between pedestal and SPE?
bool PeakValleyMetric(TH1* h, double pedpk, double speMu,
                      double& valleyX, double& valleyY,
                      double& spePeakX, double& spePeakY,
                      double& ratio)
{
  valleyX = valleyY = spePeakX = spePeakY = ratio =
    std::numeric_limits<double>::quiet_NaN();
  if (!h || !std::isfinite(pedpk) || !std::isfinite(speMu)) return false;

  const double vlo = pedpk + PV_EDGE_GAP;
  const double vhi = speMu - PV_EDGE_GAP;
  if (vhi <= vlo) return false;

  int vb0 = h->FindBin(vlo + 1e-6);
  int vb1 = h->FindBin(vhi - 1e-6);
  vb0 = std::max(1, vb0);
  vb1 = std::min(h->GetNbinsX(), vb1);
  if (vb1 < vb0) return false;

  int vbest = -1;
  double vy = 1e300;
  for (int b = vb0; b <= vb1; ++b) {
    const double y = h->GetBinContent(b);
    if (y >= 0 && y < vy) { vy = y; vbest = b; }
  }
  if (vbest < 0 || vy <= 0) return false;

  int sb0 = h->FindBin(speMu - PV_SPE_HALF_WINDOW + 1e-6);
  int sb1 = h->FindBin(speMu + PV_SPE_HALF_WINDOW - 1e-6);
  sb0 = std::max(1, sb0);
  sb1 = std::min(h->GetNbinsX(), sb1);
  if (sb1 < sb0) return false;

  int sbest = -1;
  double sy = -1.0;
  for (int b = sb0; b <= sb1; ++b) {
    const double y = h->GetBinContent(b);
    if (y > sy) { sy = y; sbest = b; }
  }
  if (sbest < 0 || sy <= 0) return false;

  valleyX = h->GetBinCenter(vbest);
  valleyY = vy;
  spePeakX = h->GetBinCenter(sbest);
  spePeakY = sy;
  ratio = sy / vy;
  return std::isfinite(ratio) && ratio > 0;
}

double ResidualFrac(const std::vector<double>& residualX,
                const std::vector<double>& residualY,
                double pedpk)
{
  int tot = 0;
  int best = 0;
  int cur = 0;
  int prev = 0;

  for (size_t i = 0; i < residualX.size(); ++i) {
    if (residualX[i] <= pedpk) continue;
    ++tot;
    const int s = (residualY[i] >= 0) ? 1 : -1;
    if (s == prev) ++cur;
    else {
      cur = 1;
      prev = s;
    }
    best = std::max(best, cur);
  }
  return (tot > 0) ? double(best) / double(tot) : 0.0;
}

TString OffStr(int off)
{
  if (off == 0) return TString("New nominal");
  return TString(Form("New nom +%d V", off));
}

// Read fitspe_aug.txt preserving the four rows/channel in run order.
//
// Backward compatible:
//   legacy: channel voltage mean mean_err status
//   v1:     channel voltage mean mean_err status amp sigma bin_width
//   v2:     channel voltage mean mean_err status amp sigma fit_lo fit_hi bin_width
//
// The updated fitspe_aug.C and fitspe_manual_aug.C write v2.  The extra
// parameters let the comparison PDF reconstruct the FINAL single-Gaussian
// fit without duplicating the single-fit algorithm here.
void ReadSingleSeeds(const TString& path,
                     std::map<int, std::vector<SeedPoint>>& seeds)
{
  std::ifstream in(path.Data());
  if (!in) {
    printf("WARNING: cannot open %s; SPE seeds will come from spectra\n",
           path.Data());
    return;
  }

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;

    std::istringstream ss(line);
    int ch = -1;
    double volt = 0, mean = 0, err = 0;
    std::string st;
    if (!(ss >> ch >> volt >> mean >> err >> st)) continue;

    std::vector<double> extra;
    double v = 0;
    while (ss >> v) extra.push_back(v);

    SeedPoint s;
    s.ok  = (st == "ok" || st == "ok_manual");
    s.off = (st == "off");
    if (s.ok && mean > 0) {
      s.mean = mean;
      s.meanErr = err;
    }

    if (s.ok && extra.size() >= 5 &&
        std::isfinite(extra[0]) && extra[0] > 0 &&
        std::isfinite(extra[1]) && extra[1] > 0 &&
        std::isfinite(extra[4]) && extra[4] > 0) {
      s.amp      = extra[0];
      s.sigma    = extra[1];
      s.fitLo    = extra[2];
      s.fitHi    = extra[3];
      s.binWidth = extra[4];
      s.hasShape = true;
      s.hasWindow = std::isfinite(s.fitLo) && std::isfinite(s.fitHi)
                    && s.fitHi > s.fitLo;
    } else if (s.ok && extra.size() >= 3 &&
               std::isfinite(extra[0]) && extra[0] > 0 &&
               std::isfinite(extra[1]) && extra[1] > 0 &&
               std::isfinite(extra[2]) && extra[2] > 0) {
      s.amp      = extra[0];
      s.sigma    = extra[1];
      s.binWidth = extra[2];
      s.hasShape = true;
      s.hasWindow = false;
    }

    seeds[ch].push_back(s);
  }
}

bool IsOffChannel(const std::map<int, std::vector<SeedPoint>>& seeds, int ch)
{
  auto it = seeds.find(ch);
  if (it == seeds.end() || it->second.empty()) return false;
  for (const SeedPoint& s : it->second)
    if (!s.off) return false;
  return true;
}

// -------------------------------------------------------------------- main --

void combfit_gamma_aug()
{
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

  const TString base   = gSystem->ExpandPathName("~/work/output/Aug_scan/");
  const TString outDir = gSystem->ExpandPathName("~/work/output/combfit_aug/");
  gSystem->mkdir(outDir, true);

  const TString seedFile = base + "fitspe_aug.txt";
  const TString outTxt   = outDir + "combfit_gamma_aug.txt";
  const TString pdfOv      = outDir + "combfit_gamma_aug.pdf";
  const TString pdfResidual    = outDir + "combfit_gamma_aug_residual.pdf";
  const TString pdfCompare = outDir + "combfit_gamma_aug_vs_single.pdf";

  std::map<int, std::vector<SeedPoint>> seeds;
  ReadSingleSeeds(seedFile, seeds);

  std::vector<TFile*> files(NRUN, nullptr);
  std::vector<TDirectory*> dirs(NRUN, nullptr);
  for (int i = 0; i < NRUN; ++i) {
    files[i] = TFile::Open(base + Form(FPAT, RUNS[i]), "READ");
    if (!files[i] || files[i]->IsZombie()) {
      printf("ERROR: cannot open run %d\n", RUNS[i]);
      return;
    }
    dirs[i] = (TDirectory*)files[i]->Get(HDIR);
    if (!dirs[i]) {
      printf("ERROR: %s missing in run %d\n", HDIR, RUNS[i]);
      return;
    }
  }

  // Only channels that exist in all four runs.
  std::set<int> common = ReadChannels(dirs[0]);
  for (int i = 1; i < NRUN; ++i) {
    std::set<int> ci = ReadChannels(dirs[i]);
    std::set<int> tmp;
    for (int ch : common)
      if (ci.count(ch)) tmp.insert(ch);
    common.swap(tmp);
  }

  std::vector<int> channels(common.begin(), common.end());
  std::sort(channels.begin(), channels.end());
  printf("Combined gamma + Gaussian fit: %zu channels x %d runs\n",
         channels.size(), NRUN);

  FILE* fout = fopen(outTxt.Data(), "w");
  if (!fout) {
    printf("ERROR: cannot write %s\n", outTxt.Data());
    return;
  }

  fprintf(fout,
    "# channel run offset ped_peak ped_init_hi fit_lo fit_hi "
    "ped_N ped_k ped_theta ped_mu ped_x0 spe_A spe_mu spe_mu_err spe_sigma "
    "chi2ndf residual_frac valley_x valley_y spe_peak_y peak_valley_ratio "
    "separation fit_status verdict quality\n");

  TCanvas cOv("cOv", "cOv", 1200, 850);
  TCanvas cResidual("cResidual", "cResidual", 1500, 1500);
  TCanvas cCmp("cCmp", "cCmp", 1200, 850);
  cOv.Print(pdfOv + "[");
  cResidual.Print(pdfResidual + "[");
  cCmp.Print(pdfCompare + "[");

  int nGood = 0;
  int nCheck = 0;
  int nNoFit = 0;
  int nBuried = 0;
  int nOff = 0;

  for (int ch : channels) {
    if (IsOffChannel(seeds, ch)) {
      ++nOff;
      continue;
    }

    std::vector<RunFit> res(NRUN);

    for (int i = 0; i < NRUN; ++i) {
      TH1* hsrc = GetH(files[i], ch);
      if (!hsrc || hsrc->Integral() <= 0) continue;
      TH1* h = (TH1*)hsrc->Clone(Form("hcomb_%d_%d", ch, RUNS[i]));
      h->SetDirectory(nullptr);
      if (REBIN > 1) h->Rebin(REBIN);

      const double pedpk = FirstPeak(h);
      if (!(pedpk > 0)) continue;

      double speSeed = std::numeric_limits<double>::quiet_NaN();
      auto sit = seeds.find(ch);
      if (sit != seeds.end() && i < (int)sit->second.size() && sit->second[i].ok)
        speSeed = sit->second[i].mean;
      if (!std::isfinite(speSeed) || speSeed < pedpk + MU_GAP)
        speSeed = SpectrumSpeSeed(h, pedpk);

      // Gamma pedestal uses the FULL pedestal.  The fit starts on the rising
      // side at a data-driven left edge rather than at the pedestal peak.
      const double lo = GammaLeftEdge(h, pedpk);

      // Full-pedestal gamma prefit -> initial N, k, theta only.
      double pedInitHi = pedpk + PED_PREFIT_MAX_RIGHT;
      GammaSeedFit gs = FullPedGammaSeed(h, pedpk, speSeed, lo, pedInitHi);

      const double hi = std::min(FIT_XMAX,
                         std::max(speSeed + SPE_RIGHT_SPAN,
                                  pedpk + MIN_COMB_SPAN));

      const int nb = h->FindBin(hi - 1e-6) - h->FindBin(lo + 1e-6) + 1;
      if (nb < MIN_BINS || hi <= lo) continue;

      // Estimate the SPE amplitude after subtracting the gamma pedestal seed.
      const double pedAtSpe = GammaValue(speSeed, gs.N, gs.k, gs.theta, gs.muPed);
      const double ySpe = h->GetBinContent(h->FindBin(speSeed));
      const double speA0 = std::max(1.0, ySpe - pedAtSpe);

      const double muLo = pedpk + MU_GAP;
      const double muHi = std::min(SPE_SEARCH_HI, hi - 0.2);
      if (muHi <= muLo) continue;
      speSeed = std::min(std::max(speSeed, muLo + 0.1), muHi - 0.1);

      TF1 fit(Form("fcomb_%d_%d", ch, RUNS[i]), GammaGauss, lo, hi, 7);
      fit.SetParameters(gs.N, gs.k, gs.theta, pedpk,
                        speA0, speSeed, SPE_SIGMA_SEED);

      fit.SetParLimits(0, 0.0, 1e6 * std::max(1.0, gs.N) + 1e3);
      fit.SetParLimits(1, GAMMA_K_MIN, GAMMA_K_MAX);
      fit.SetParLimits(2, GAMMA_THETA_MIN, GAMMA_THETA_MAX);
      if (FIX_PED_PEAK)
        fit.FixParameter(3, pedpk);
      else
        fit.SetParLimits(3, pedpk - PED_PEAK_FLOAT, pedpk + PED_PEAK_FLOAT);
      fit.SetParLimits(4, 0.0, 50.0 * std::max(1.0, speA0) + 1e3);
      fit.SetParLimits(5, muLo, muHi);
      fit.SetParLimits(6, SPE_SIGMA_MIN, SPE_SIGMA_MAX);

      TFitResultPtr fr = h->Fit(&fit, "LQNRS");
      const int fitStatus = fr.Get() ? fr->Status() : -1;

      const double pedN     = fit.GetParameter(0);
      const double pedK     = fit.GetParameter(1);
      const double pedTheta = fit.GetParameter(2);
      const double pedMu    = fit.GetParameter(3);
      const double pedX0    = pedMu - (pedK - 1.0) * pedTheta;
      const double speA     = fit.GetParameter(4);
      const double speMu    = fit.GetParameter(5);
      const double speErr   = fit.GetParError(5);
      const double speSig   = fit.GetParameter(6);

      const bool finite = std::isfinite(pedN) && std::isfinite(pedK) &&
                          std::isfinite(pedTheta) && std::isfinite(pedMu) &&
                          std::isfinite(pedX0) && std::isfinite(speA) &&
                          std::isfinite(speMu) && std::isfinite(speErr) &&
                          std::isfinite(speSig);
      const bool usable = (fit.GetNDF() > 0) && finite;
      const double c2 = usable ? fit.GetChisquare() / fit.GetNDF() : 0.0;

      std::vector<double> residualX, residualY;
      if (usable) {
        double pp[7] = {pedN, pedK, pedTheta, pedMu, speA, speMu, speSig};
        for (int b = h->FindBin(lo + 1e-6);
             b <= h->FindBin(hi - 1e-6); ++b) {
          const double y = h->GetBinContent(b);
          if (y <= 0) continue;
          const double x = h->GetBinCenter(b);
          double xx[1] = {x};
          residualX.push_back(x);
          residualY.push_back((y - GammaGauss(xx, pp)) / std::sqrt(y));
        }
      }

      const double residualFrac = usable ? ResidualFrac(residualX, residualY, pedpk) : 0.0;
      const double sep = usable ? (speMu - pedMu) : 0.0;

      double valleyX = std::numeric_limits<double>::quiet_NaN();
      double valleyY = std::numeric_limits<double>::quiet_NaN();
      double spePeakX = std::numeric_limits<double>::quiet_NaN();
      double spePeakY = std::numeric_limits<double>::quiet_NaN();
      double pvRatio  = std::numeric_limits<double>::quiet_NaN();
      const bool pvValid = usable && PeakValleyMetric(
        h, pedpk, speMu, valleyX, valleyY, spePeakX, spePeakY, pvRatio);

      const char* verdict = "NO_PV";
      if (pvValid) verdict = (pvRatio < PV_RATIO_MIN) ? "BURIED" : "CLEAR";

      bool kEdge = false;
      bool thetaEdge = false;
      bool pedMuEdge = false;
      bool speMuEdge = false;
      bool speSigEdge = false;
      if (usable) {
        kEdge = (pedK < GAMMA_K_MIN + EDGE_EPS ||
                 pedK > GAMMA_K_MAX - EDGE_EPS);
        thetaEdge = (pedTheta < GAMMA_THETA_MIN + EDGE_EPS ||
                     pedTheta > GAMMA_THETA_MAX - EDGE_EPS);
        pedMuEdge = (!FIX_PED_PEAK) &&
                    (pedMu < pedpk - PED_PEAK_FLOAT + EDGE_EPS ||
                     pedMu > pedpk + PED_PEAK_FLOAT - EDGE_EPS);
        speMuEdge = (speMu < muLo + EDGE_EPS || speMu > muHi - EDGE_EPS);
        speSigEdge = (speSig < SPE_SIGMA_MIN + EDGE_EPS ||
                      speSig > SPE_SIGMA_MAX - EDGE_EPS);
      }

      TString quality = "nofit";
      if (usable) {
        const bool bad = (fitStatus != 0) ||
                         (residualFrac > RESIDUALFRAC_BAD) ||
                         kEdge || thetaEdge || pedMuEdge ||
                         speMuEdge || speSigEdge;
        quality = bad ? "check" : "good";
        if (bad) ++nCheck;
        else ++nGood;
        if (pvValid && pvRatio < PV_RATIO_MIN) ++nBuried;
      } else {
        ++nNoFit;
      }

      fprintf(fout,
        "%4d %6d %4d %8.3f %8.3f %8.3f %8.3f "
        "%12.4g %8.3f %8.3f %8.3f %8.3f %12.4g %8.3f %8.3f %8.3f "
        "%9.2f %8.3f %8.3f %12.4g %12.4g %8.3f %8.3f %4d %s %s\n",
        ch, RUNS[i], OFFSET[i], pedpk, pedInitHi, lo, hi,
        pedN, pedK, pedTheta, pedMu, pedX0,
        speA, speMu, speErr, speSig,
        c2, residualFrac, valleyX, valleyY, spePeakY, pvRatio, sep,
        fitStatus, verdict, quality.Data());
      fflush(fout);

      RunFit& r = res[i];
      r.filled = true;
      r.usable = usable;
      r.h = h;
      r.run = RUNS[i];
      r.offset = OFFSET[i];
      r.fitStatus = fitStatus;
      r.pedpk = pedpk;
      r.pedInitHi = pedInitHi;
      r.lo = lo;
      r.hi = hi;
      r.pedN = pedN;
      r.pedK = pedK;
      r.pedTheta = pedTheta;
      r.pedMu = pedMu;
      r.pedX0 = pedX0;
      r.speA = speA;
      r.speMu = speMu;
      r.speMuErr = speErr;
      r.speSig = speSig;
      r.c2 = c2;
      r.residualFrac = residualFrac;
      r.sep = sep;
      r.valleyX = valleyX;
      r.valleyY = valleyY;
      r.spePeakX = spePeakX;
      r.spePeakY = spePeakY;
      r.pvRatio = pvRatio;
      r.pvValid = pvValid;
      r.residualX = residualX;
      r.residualY = residualY;
      r.verdict = verdict;
      r.quality = quality;
    }

    // ------------------------------------------------ overview 2x2 per channel
    cOv.Clear();
    cOv.Divide(2, 2, 0.006, 0.012);

    for (int i = 0; i < NRUN; ++i) {
      cOv.cd(i + 1);
      gPad->SetLeftMargin(0.14);
      gPad->SetBottomMargin(0.13);
      gPad->SetTopMargin(0.12);

      RunFit& r = res[i];
      if (!r.filled) continue;

      TH1* h = r.h;
      h->SetLineColor(kBlack);
      h->SetLineWidth(1);
      h->SetTitle(Form("ch %d | %s;Amplitude [ADC];Counts",
                       ch, OffStr(r.offset).Data()));
      h->GetXaxis()->SetRangeUser(0, std::min(FIT_XMAX, r.hi + 5.0));
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitleSize(0.05);
      h->GetYaxis()->SetTitleSize(0.05);
      h->SetMinimum(0);
      h->Draw("HIST");

      if (r.usable) {
        TF1* ft = new TF1(Form("tot_%d_%d", ch, r.run),
                          GammaGauss, r.lo, r.hi, 7);
        ft->SetParameters(r.pedN, r.pedK, r.pedTheta, r.pedMu,
                          r.speA, r.speMu, r.speSig);
        ft->SetNpx(400);
        ft->SetLineColor(kRed + 1);
        ft->SetLineWidth(2);
        ft->Draw("SAME");

        TF1* fp = new TF1(Form("ped_%d_%d", ch, r.run),
                          GammaComp, r.lo, r.hi, 4);
        fp->SetParameters(r.pedN, r.pedK, r.pedTheta, r.pedMu);
        fp->SetNpx(400);
        fp->SetLineColor(kBlue + 1);
        fp->SetLineStyle(2);
        fp->SetLineWidth(2);
        fp->Draw("SAME");

        TF1* fs = new TF1(Form("spe_%d_%d", ch, r.run),
                          SpeComp, r.lo, r.hi, 3);
        fs->SetParameters(r.speA, r.speMu, r.speSig);
        fs->SetNpx(400);
        fs->SetLineColor(kGreen + 2);
        fs->SetLineStyle(2);
        fs->SetLineWidth(2);
        fs->Draw("SAME");

        TLine* lpk = new TLine(r.pedpk, 0, r.pedpk, h->GetMaximum());
        lpk->SetLineColor(kMagenta + 1);
        lpk->SetLineStyle(3);
        lpk->Draw("SAME");

        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.048);
        tx.DrawLatex(0.47, 0.82,
          Form("SPE #mu = %.2f #pm %.2f", r.speMu, r.speMuErr));
        tx.SetTextSize(0.043);
        tx.DrawLatex(0.47, 0.75,
          Form("ped k = %.2f, #theta = %.2f", r.pedK, r.pedTheta));
        tx.SetTextSize(0.048);
        tx.SetTextColor(r.verdict == "BURIED" ? kOrange + 2 :
                        (r.verdict == "CLEAR" ? kGreen + 2 : kGray + 2));
        if (r.pvValid)
          tx.DrawLatex(0.47, 0.68,
            Form("%s, P/V = %.2f", r.verdict.Data(), r.pvRatio));
        else
          tx.DrawLatex(0.47, 0.68, "P/V unavailable");
        if (r.quality == "check") {
          tx.SetTextColor(kRed + 1);
          tx.DrawLatex(0.47, 0.61, "check");
        }
      } else {
        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.06);
        tx.SetTextColor(kRed + 1);
        tx.DrawLatex(0.35, 0.5, "no fit");
      }
    }

    cOv.cd(0);
    TLatex head;
    head.SetNDC();
    head.SetTextSize(0.025);
    head.DrawLatex(0.01, 0.985,
      Form("Gamma pedestal + Gaussian SPE, channel %d", ch));
    cOv.Print(pdfOv);

    // -------------------------------- comparison with the FINAL single Gaussian
    // The final single-fit parameters are read from fitspe_aug.txt.  The manual
    // macro patches that same file, so this follows the final auto/manual result
    // without copying the single-fit algorithm.
    //
    // The comparison histogram keeps the COMBFIT binning (r.h; REBIN=1 by
    // default).  The saved single-Gaussian curve is drawn on that same histogram
    // by converting only its amplitude from its original counts/bin convention
    // to the combfit bin width.  Its fitted mu, sigma, uncertainty and original
    // fit window are unchanged.
    cCmp.Clear();
    cCmp.Divide(2, 2, 0.006, 0.012);

    for (int i = 0; i < NRUN; ++i) {
      cCmp.cd(i + 1);
      gPad->SetLeftMargin(0.14);
      gPad->SetBottomMargin(0.13);
      gPad->SetTopMargin(0.10);

      RunFit& r = res[i];
      if (!r.filled) continue;

      const SeedPoint* sp = nullptr;
      auto sit = seeds.find(ch);
      if (sit != seeds.end() && i < (int)sit->second.size() &&
          sit->second[i].ok)
        sp = &sit->second[i];

      // Comparison plot deliberately uses the SAME histogram binning as the
      // combined fit (r.h).  With the default REBIN=1 this is the original
      // unrebinned August histogram.  The saved single-Gaussian amplitude is
      // converted from its own counts/bin convention to this display bin width
      // below; its mu, sigma and original fit window are never changed.
      TH1* hCmp = (TH1*)r.h->Clone(Form("hcmp_%d_%d", ch, r.run));
      hCmp->SetDirectory(nullptr);

      hCmp->SetLineColor(kBlack);
      hCmp->SetLineWidth(1);
      hCmp->SetTitle(Form("ch %d | %s;Amplitude [ADC];Counts",
                          ch, OffStr(r.offset).Data()));

      double drawHi = std::min(FIT_XMAX, r.hi + 5.0);
      if (sp && sp->hasWindow)
        drawHi = std::min(FIT_XMAX, std::max(drawHi, sp->fitHi + 2.0));

      hCmp->GetXaxis()->SetRangeUser(0, drawHi);
      hCmp->GetXaxis()->SetLabelSize(0.05);
      hCmp->GetYaxis()->SetLabelSize(0.05);
      hCmp->GetXaxis()->SetTitleSize(0.05);
      hCmp->GetYaxis()->SetTitleSize(0.05);
      hCmp->SetMinimum(0);
      hCmp->Draw("HIST");

      const int bDisp = hCmp->FindBin(
        (sp && std::isfinite(sp->mean)) ? sp->mean : r.speMu);
      const double displayBW = hCmp->GetBinWidth(bDisp);
      const int bComb = r.h->FindBin(r.speMu);
      const double combBW = r.h->GetBinWidth(bComb);

      TF1* fCombDraw = nullptr;
      if (r.usable && displayBW > 0 && combBW > 0) {
        const double combScale = displayBW / combBW;
        fCombDraw = new TF1(Form("cmp_comb_%d_%d", ch, r.run),
                            GammaGauss, r.lo, r.hi, 7);
        fCombDraw->SetParameters(r.pedN * combScale, r.pedK, r.pedTheta, r.pedMu,
                                 r.speA * combScale, r.speMu, r.speSig);
        fCombDraw->SetNpx(400);
        fCombDraw->SetLineColor(kRed + 1);
        fCombDraw->SetLineStyle(1);
        fCombDraw->SetLineWidth(2);
        fCombDraw->Draw("SAME");
      }

      TF1* fSingleDraw = nullptr;
      if (sp && sp->hasShape && displayBW > 0) {
        // ROOT histogram amplitudes are counts/bin.  If the original single fit
        // used a different bin width (e.g. Rebin(2)), convert only the vertical
        // normalization to the current combfit binning.  This does NOT alter
        // the fitted single-Gaussian mu, sigma, error, or fit window.
        double singleAmpDraw = sp->amp;
        if (sp->binWidth > 0)
          singleAmpDraw *= displayBW / sp->binWidth;

        const double singleLo = sp->hasWindow ? sp->fitLo
                                              : std::max(0.0, sp->mean - 5.0);
        const double singleHi = sp->hasWindow ? sp->fitHi
                                              : std::min(drawHi, sp->mean + 5.0);

        fSingleDraw = new TF1(Form("cmp_single_%d_%d", ch, r.run),
                              "gaus", singleLo, singleHi);
        fSingleDraw->SetParameters(singleAmpDraw, sp->mean, sp->sigma);
        fSingleDraw->SetNpx(400);
        fSingleDraw->SetLineColor(kBlue + 1);
        fSingleDraw->SetLineStyle(1);
        fSingleDraw->SetLineWidth(2);
        fSingleDraw->Draw("SAME");
      }

      const double txtX = 0.45;
      const double txtY0 = 0.84;
      const double txtDY = 0.060;

      TLegend* leg = new TLegend(0.54, 0.54, 0.88, 0.65);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.032);
      if (fCombDraw) leg->AddEntry(fCombDraw, "combined fit", "l");
      if (fSingleDraw) leg->AddEntry(fSingleDraw, "single Gaussian fit", "l");
      leg->Draw();

      TLatex tx;
      tx.SetNDC();
      tx.SetTextAlign(13);
      tx.SetTextSize(0.040);

      if (r.usable) {
        tx.SetTextColor(kRed + 1);
        tx.DrawLatex(txtX, txtY0,
          Form("comb: #mu = %.2f #pm %.2f", r.speMu, r.speMuErr));
      }

      if (sp) {
        tx.SetTextColor(kBlue + 1);
        tx.DrawLatex(txtX, txtY0 - txtDY,
          Form("single: #mu = %.2f #pm %.2f", sp->mean, sp->meanErr));

        if (r.usable && std::isfinite(sp->mean)) {
          tx.SetTextColor(kBlack);
          tx.DrawLatex(txtX, txtY0 - 2.0 * txtDY,
            Form("#Delta#mu (comb #minus single) = %.2f ADC",
                 r.speMu - sp->mean));
        }

        if (!sp->hasShape) {
          tx.SetTextColor(kGray + 2);
          tx.SetTextSize(0.028);
          tx.DrawLatex(txtX, txtY0 - 3.0 * txtDY,
            "single curve unavailable: rerun updated fitspe macros");
        } else if (!sp->hasWindow) {
          tx.SetTextColor(kGray + 2);
          tx.SetTextSize(0.028);
          tx.DrawLatex(txtX, txtY0 - 3.0 * txtDY,
            "single window unavailable: using fallback display range");
        }
      } else {
        tx.SetTextColor(kGray + 2);
        tx.SetTextSize(0.032);
        tx.DrawLatex(txtX, txtY0 - txtDY, "single fit unavailable");
      }
    }

    cCmp.cd(0);
    TLatex cmpHead;
    cmpHead.SetNDC();
    cmpHead.SetTextSize(0.025);
    cmpHead.DrawLatex(0.01, 0.985,
      Form("Gamma+Gaussian fit vs final single-Gaussian fit, channel %d", ch));
    cCmp.Print(pdfCompare);

    // ---------------------------------------- all fits with residual, 4 rows x 2 columns
    // One page per channel.  Each row is one August voltage point:
    // left = spectrum + fit components, right = residual.
    cResidual.Clear();
    cResidual.Divide(2, 4, 0.004, 0.006);
    bool anyResidual = false;

    for (int i = 0; i < NRUN; ++i) {
      RunFit& r = res[i];
      if (!r.filled || !r.usable) continue;
      anyResidual = true;

      cResidual.cd(2 * i + 1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.13);
      gPad->SetTopMargin(0.10);
      r.h->SetLineColor(kBlack);
      r.h->SetLineWidth(1);
      r.h->SetTitle(Form("ch %d | %s;Amplitude [ADC];Counts",
                         ch, OffStr(r.offset).Data()));
      r.h->GetXaxis()->SetRangeUser(0, std::min(FIT_XMAX, r.hi + 5.0));
      r.h->SetMinimum(0);
      r.h->Draw("HIST");

      TF1* ftAll = new TF1(Form("ptot_%d_%d", ch, r.run),
                           GammaGauss, r.lo, r.hi, 7);
      ftAll->SetParameters(r.pedN, r.pedK, r.pedTheta, r.pedMu,
                           r.speA, r.speMu, r.speSig);
      ftAll->SetNpx(400);
      ftAll->SetLineColor(kRed + 1);
      ftAll->SetLineWidth(2);
      ftAll->Draw("SAME");

      TF1* fpAll = new TF1(Form("pped_%d_%d", ch, r.run),
                           GammaComp, r.lo, r.hi, 4);
      fpAll->SetParameters(r.pedN, r.pedK, r.pedTheta, r.pedMu);
      fpAll->SetNpx(400);
      fpAll->SetLineColor(kBlue + 1);
      fpAll->SetLineStyle(2);
      fpAll->SetLineWidth(1);
      fpAll->Draw("SAME");

      TF1* fsAll = new TF1(Form("pspe_%d_%d", ch, r.run),
                           SpeComp, r.lo, r.hi, 3);
      fsAll->SetParameters(r.speA, r.speMu, r.speSig);
      fsAll->SetNpx(400);
      fsAll->SetLineColor(kGreen + 2);
      fsAll->SetLineStyle(2);
      fsAll->SetLineWidth(1);
      fsAll->Draw("SAME");

      TLatex ptxt;
      ptxt.SetNDC();
      ptxt.SetTextSize(0.050);
      ptxt.DrawLatex(0.50, 0.83,
        Form("#mu_{SPE}=%.2f #pm %.2f", r.speMu, r.speMuErr));
      ptxt.SetTextSize(0.045);
      if (r.pvValid)
        ptxt.DrawLatex(0.50, 0.75,
          Form("P/V=%.2f  (%s)", r.pvRatio, r.verdict.Data()));
      ptxt.SetTextColor(r.quality == "check" ? kRed + 1 : kGray + 2);
      ptxt.DrawLatex(0.50, 0.67, r.quality.Data());

      cResidual.cd(2 * i + 2);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.13);
      gPad->SetTopMargin(0.10);
      if (r.residualX.size() > 1) {
        TGraph* grAll = new TGraph((int)r.residualX.size(), &r.residualX[0], &r.residualY[0]);
        grAll->SetTitle(Form("Residual | %s;Amplitude [ADC];(data-fit)/#sqrt{data}",
                             OffStr(r.offset).Data()));
        grAll->SetMarkerStyle(20);
        grAll->SetMarkerSize(0.45);
        grAll->SetMarkerColor(kBlue + 1);
        grAll->Draw("AP");
        gPad->Update();

        const double uy0 = gPad->GetUymin();
        const double uy1 = gPad->GetUymax();
        TLine* l0All = new TLine(r.lo, 0, r.hi, 0);
        l0All->SetLineColor(kGray + 2);
        l0All->SetLineStyle(2);
        l0All->Draw();

        TLine* lp1 = new TLine(r.lo, +1, r.hi, +1);
        lp1->SetLineColor(kGray + 1); lp1->SetLineStyle(3); lp1->Draw();
        TLine* lm1 = new TLine(r.lo, -1, r.hi, -1);
        lm1->SetLineColor(kGray + 1); lm1->SetLineStyle(3); lm1->Draw();

        if (r.pvValid) {
          TLine* lvAll = new TLine(r.valleyX, uy0, r.valleyX, uy1);
          lvAll->SetLineColor(kOrange + 7);
          lvAll->SetLineStyle(3);
          lvAll->Draw();
        }

        TLine* lsAll = new TLine(r.speMu, uy0, r.speMu, uy1);
        lsAll->SetLineColor(kGreen + 2);
        lsAll->SetLineStyle(2);
        lsAll->Draw();
      }
    }
    if (anyResidual) cResidual.Print(pdfResidual);

  }

  cOv.Print(pdfOv + "]");
  cResidual.Print(pdfResidual + "]");
  cCmp.Print(pdfCompare + "]");

  fclose(fout);
  for (TFile* f : files)
    if (f) f->Close();

  printf("Done: %d good, %d check, %d nofit, %d BURIED points, %d off channels skipped\n",
         nGood, nCheck, nNoFit, nBuried, nOff);
  printf("  overview -> %s\n", pdfOv.Data());
  printf("  residuals -> %s\n", pdfResidual.Data());
  printf("  compare  -> %s\n", pdfCompare.Data());
  printf("  table    -> %s\n", outTxt.Data());
}

