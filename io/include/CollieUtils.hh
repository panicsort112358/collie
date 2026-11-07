#ifndef CollieUtils_hh_included
#define CollieUtils_hh_included

#include "CollieHistogram.hh"
#include "CollieQuickKEYS.hh"
#include "CollieDistribution.hh"
#include "CollieChannel.hh"
#include "CollieLoader.hh"
#include <iostream>

class CollieUtils {

public:
  CollieUtils();
  ~CollieUtils();

  inline void setCutoffs(float low, float high) {
    m_cutLow = low;
    m_cutHigh = high;
    return;
  }
  inline void setInputHist(float min, float max, int bins) {
    m_histMin = min;
    m_histMax = max;
    m_histBins = bins;
    return;
  }
  inline void setRebin(int rebin) {
    m_rebin = rebin;
    return;
  }
  inline void clear() {
    m_rebin = 1;
    m_cutLow = -1234;
    m_cutHigh = -1234;
    m_histMin = -1234;
    m_histMax = -1234;
    m_histBins = -1234;
    return;
  }
  bool check();
  bool getSmoothedROOT(float alpha, int iter, TH1F* in, TH1F& out);
  bool getSmoothedROOT(float alpha, int iter, CollieHistogram& in, TH1F& out);
  bool getSmoothedCollie(float alpha, int iter, TH1F* in, CollieHistogram& out);
  bool getSmoothedCollie(float alpha, int iter, CollieHistogram& in, CollieHistogram& out);
  bool getROOT(CollieHistogram& in, TH1F& out);
  bool getCollie(TH1F* hist, CollieHistogram& hist);

  void fillDist(CollieHistogram* hist, CollieDistribution* dist);
  void massPoint(int mass, CollieChannel* chan, CollieHistogram* sig,
                 std::vector<CollieHistogram*>& bkgd, CollieHistogram* data);

  void interpolate(CollieHistogram& h1, const double p1, const double xsec1, CollieHistogram& h2,
                   const double p2, const double xsec2, CollieHistogram& hf, const double pf,
                   const double xsecF);

  CollieLoader getLoader(const char* filename, const char* options, bool ok);

  TH1D* h1FMT(TH1D* hist);

private:
  float m_cutLow;
  float m_cutHigh;
  int m_rebin;
  float m_histMin;
  float m_histMax;
  int m_histBins;
  int m_count;
};

#endif // CollieUtils_hh_included
