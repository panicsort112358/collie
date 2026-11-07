// CollieHistogram.cpp: implementation of the CollieHistogram class.
//
//////////////////////////////////////////////////////////////////////

#include "CollieHistogram.hh"

#ifndef MIN
/// Returns the smaller of its two arguments
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef MAX
/// Returns the larger of its two arguments
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifndef PI
#define PI 3.141592654
#endif

#ifndef SQR
/// Returns the square of its argument
#define SQR(x) ((x) * (x))
#endif

#ifndef ABS
/// Returns the fast absolute value of its argument
#define ABS(x) ((x) < 0 ? (-(x)) : (x))
#endif

#define NORTH2D 0
#define SOUTH2D 1
#define WEST2D 2
#define EAST2D 3
#define NORTHEAST2D 4
#define SOUTHEAST2D 5
#define NORTHWEST2D 6
#define SOUTHWEST2D 7

/*
All CollieHistograms share one list...
*/
CollieHistogramContainer CollieHistogram::DefaultContainer;

///
// Histogram container
///

bool CollieHistogramContainer::AddHistogram(CollieHistogram* p_hist) {
  p_histograms.push_back(p_hist);
  return true;
}

bool CollieHistogramContainer::AddROOT(TH1* hist) {
  if (hist == NULL)
    return false;
  p_extern.push_back(hist);
  return true;
}

void CollieHistogramContainer::ClearHistograms() {
  p_extern.clear();
  p_histograms.clear();
}

bool CollieHistogramContainer::RemoveHistogram(CollieHistogram* p_hist) {
  vector<CollieHistogram*>::iterator i;
  for (i = p_histograms.begin(); i != p_histograms.end() && (*i) != p_hist; i++)
    ;
  if (i != p_histograms.end()) {
    p_histograms.erase(i);
    return true;
  } else
    return false;
}

int CollieHistogramContainer::GetUniqueId(int id) {
  // need lock here...
  unsigned int i;

  if (id > 0) {
    for (i = 0; i < p_histograms.size(); i++)
      if (p_histograms[i] != NULL && p_histograms[i]->Id() == id)
        break;
    if (i == p_histograms.size())
      return id;
  } else {
    for (int testid = 100000; testid < 2000000; testid++) {
      for (i = 0; i < p_histograms.size(); i++)
        if (p_histograms[i] != NULL && p_histograms[i]->Id() == testid)
          break;
      if (i == p_histograms.size())
        return testid;
    }
  }
  return 0;
}

// returns true if x preceeds y in the list
static bool StrictWeak(CollieHistogram* x, CollieHistogram* y) {
  if (x == y)
    return false;
  if (x == NULL)
    return false;
  if (y == NULL)
    return true;
  if (x->Id() < y->Id())
    return true;
  return false;
}

void CollieHistogramContainer::StoreHistograms(const char* filename, const char* filetype) {
  if (p_histograms.size() == 0)
    return;
  sort(p_histograms.begin(), p_histograms.end(), StrictWeak);
  CollieHistogramStorer* storer = CollieHistogramStorer::GetStorerFor(filetype);
  if (storer == NULL)
    return;
  storer->StoreHistograms(filename, p_histograms, p_extern);
}

class CollieHistogramStorerROOT : public CollieHistogramStorer {
public:
  virtual void StoreHistograms(const char* filename, vector<CollieHistogram*>& hists) {
    vector<CollieHistogram*>::iterator i;
    TFile* hfile;
    hfile = new TFile(filename, "RECREATE", "Collie Root Histogram File");
    hfile->cd();
    for (i = hists.begin(); i != hists.end(); i++)
      if (*i != NULL)
        StoreROOT(*i);
    // Save all objects in this file
    hfile->Write();
    hfile->Close();
  }

  virtual void StoreHistograms(const char* filename, vector<CollieHistogram*> hists,
                               vector<TH1*> rHists) {

    vector<CollieHistogram*>::iterator i;

    vector<TH1*>::iterator j;

    TFile* hfile;

    hfile = new TFile(filename, "RECREATE", "Collie Root Histogram File");
    if (!hfile) {
      printf("StoreHistograms: Error, could not create input file!\n");
      return;
    }
    hfile->cd();

    for (i = hists.begin(); i != hists.end(); i++)
      if (*i != NULL)
        StoreROOT(*i);

    // Save all objects in this file

    for (j = rHists.begin(); j != rHists.end(); j++) {
      if ((*j) == NULL)
        continue;
      (*j)->Write();
    }

    hfile->Write();
    hfile->Close();
  }

  void StoreROOT(CollieHistogram* hist) {
    if (!hist->Valid())
      return;
    if (!strcmp(hist->getType(), "CollieHistogram")) {
      TH1D* myHist;
      if (hist->GetName() == NULL)
        myHist = new TH1D("", "", hist->nBins(), hist->MinEdge(), hist->MaxEdge());
      else
        myHist = new TH1D(hist->GetName(), hist->GetName(), hist->nBins(), hist->MinEdge(),
                          hist->MaxEdge());
      myHist->Sumw2();
      for (int i = 0; i < hist->nBins(); i++) {
        double pos =
            hist->MinEdge() + (hist->MaxEdge() - hist->MinEdge()) / (hist->nBins()) * (i + 0.5);
        int fills = hist->FillsInBin(i);
        //  double val=(fills>0)?hist->InBin(i)/fills:0;
        //  for (int j=0; j<fills; j++) myHist->Fill(pos,val);
        int iBin = myHist->FindBin(pos);
        double val = (fills > 0) ? hist->InBin(i) : 0;
        iBin = i + 1;
        myHist->SetBinContent(iBin, val);
        if (hist->BinErr(i) > 0)
          myHist->SetBinError(iBin, hist->BinErr(i));
      }
      myHist->SetBinContent(0, hist->Underflow());
      myHist->SetBinContent(hist->nBins() + 1, hist->Overflow());

    } else if (!strcmp(hist->getType(), "CollieHistogram2d")) {
      CollieHistogram2d* h2d = (CollieHistogram2d*)hist;
      TH2D* myHist;

      if (h2d->GetName() == NULL)
        myHist = new TH2D("", "", h2d->nxBins(), h2d->LeftEdge(), h2d->RightEdge(), h2d->nyBins(),
                          h2d->BottomEdge(), h2d->TopEdge());
      else
        myHist = new TH2D(h2d->GetName(), h2d->GetName(), h2d->nxBins(), h2d->LeftEdge(),
                          h2d->RightEdge(), h2d->nyBins(), h2d->BottomEdge(), h2d->TopEdge());

      // store data...
      double ystep = double((h2d->TopEdge() - h2d->BottomEdge()) / (h2d->nyBins()));
      double xstep = double((h2d->RightEdge() - h2d->LeftEdge()) / (h2d->nxBins()));
      for (int i = 0; i < h2d->nxBins(); i++)
        for (int j = 0; j < h2d->nyBins(); j++) {
          int loc = i * h2d->nyBins() + j;
          int fills = hist->FillsInBin(loc);
          double val = (fills > 0) ? hist->InBin(loc) / fills : 0;
          double xpos = double(i) * xstep + h2d->LeftEdge() + xstep / 2.0;
          double ypos = double(j) * ystep + h2d->BottomEdge() + ystep / 2.0;
          for (int k = 0; k < fills; k++)
            myHist->Fill(xpos, ypos, val);

          int iBin = myHist->FindBin(xpos, ypos);
          if (h2d->BinErr(i, j) > 0)
            myHist->SetBinError(iBin, h2d->BinErr(i, j));
        }
      ystep = double((h2d->TopEdge() - h2d->BottomEdge()) / (4.0 * h2d->nyBins()));
      xstep = double((h2d->RightEdge() - h2d->LeftEdge()) / (4.0 * h2d->nxBins()));

      myHist->Fill(double(h2d->LeftEdge()) - xstep, double(h2d->BottomEdge()) + ystep,
                   double(h2d->Overflow(WEST2D)));
      myHist->Fill(double(h2d->RightEdge()) + xstep, double(h2d->TopEdge()) - ystep,
                   double(h2d->Overflow(EAST2D)));

      myHist->Fill(double(h2d->LeftEdge()) + xstep, double(h2d->BottomEdge()) - ystep,
                   double(h2d->Overflow(SOUTH2D)));
      myHist->Fill(double(h2d->RightEdge()) - xstep, double(h2d->TopEdge()) + ystep,
                   double(h2d->Overflow(NORTH2D)));

      myHist->Fill(double(h2d->LeftEdge()) - xstep, double(h2d->BottomEdge()) - ystep,
                   double(h2d->Overflow(SOUTHWEST2D)));
      myHist->Fill(double(h2d->RightEdge()) + xstep, double(h2d->TopEdge()) + ystep,
                   double(h2d->Overflow(NORTHEAST2D)));

      myHist->Fill(double(h2d->LeftEdge()) - xstep, double(h2d->TopEdge()) + ystep,
                   double(h2d->Overflow(NORTHWEST2D)));
      myHist->Fill(double(h2d->RightEdge()) + xstep, double(h2d->BottomEdge()) - ystep,
                   double(h2d->Overflow(SOUTHEAST2D)));
    }
  }
} gl_ROOTStorer;

map<string, CollieHistogramStorer*> CollieHistogramStorer::gl_stores;
CollieHistogramStorer* CollieHistogramStorer::GetStorerFor(const char* filetype) {
  if (filetype == NULL) { // default
    return GetStorerFor("ROOT");
    return NULL;
  }
  map<string, CollieHistogramStorer*>::iterator i = gl_stores.find((string)filetype);
  if (i == gl_stores.end())
    return NULL;
  return i->second;
}
void CollieHistogramStorer::AddStorer(const char* filetype, CollieHistogramStorer* storer) {
  map<string, CollieHistogramStorer*>::iterator i = gl_stores.find((string)filetype);
  if (i != gl_stores.end())
    return;
  gl_stores[filetype] = storer;
}

//////////////////////////////////////////////////////////////////////
// Construction/Destruction/Booking
//////////////////////////////////////////////////////////////////////

static void Setup() {
  if (CollieHistogramStorer::GetStorerFor("ROOT") == NULL)
    CollieHistogramStorer::AddStorer("ROOT", &gl_ROOTStorer);
}

CollieHistogram::CollieHistogram() {
  Setup();
  p_Histogram = NULL;
  p_Fills = NULL;
  p_Errs = NULL;
  p_container = NULL;
  m_end = -1;
  m_start = 0;
  n_bins = 0;
  m_hist = -1;
  m_name = NULL;

  n_entries = 0;
  m_offleft = 0.0;
  m_offright = 0.0;
}

CollieHistogram::CollieHistogram(const CollieHistogram& h) {
  CollieHistogram();
  p_Histogram = NULL;
  p_Fills = NULL;
  p_Errs = NULL;
  p_container = NULL;
  m_end = -1;
  m_start = 0;
  n_bins = 0;
  m_hist = -1;
  m_name = NULL;

  n_entries = 0;
  m_offleft = 0.0;
  m_offright = 0.0;
  *this = h;
}

/*
CollieHistogram::CollieHistogram(int id, const char* name, int bins, double left, double right, const char* dir)
{
  CollieHistogram();
  Book(id,name,bins,left,right,dir);
}

CollieHistogram::CollieHistogram(const char* name, int bins, double left, double right, const char* dir)
{
  CollieHistogram();
  Book(name,bins,left,right,dir);
}
*/

CollieHistogram::~CollieHistogram() {
  // NEED LOCK HERE
  if (p_container != NULL)
    p_container->RemoveHistogram(this);
  // NEED UNLOCK HERE
  if (p_Histogram != NULL) {
    delete[] p_Histogram;
    p_Histogram = NULL;
  }
  if (p_Fills != NULL) {
    delete[] p_Fills;
    p_Fills = NULL;
  }
  if (p_Errs != NULL) {
    delete[] p_Errs;
    p_Errs = NULL;
  }
  if (m_name != NULL) {
    delete[] m_name;
    m_name = NULL;
  }
}

bool CollieHistogram::BookROOT(TH1* inhist, CollieHistogramContainer* container) {
  if (p_container != NULL)
    p_container->RemoveHistogram(this);
  p_container = (container == NULL) ? (&DefaultContainer) : container;
  return p_container->AddROOT(inhist);
}

int CollieHistogram::Book(int id, const char* name, int bins, double left, double right,
                          CollieHistogramContainer* container) {
  if (p_container != NULL)
    p_container->RemoveHistogram(this);
  p_container = (container == NULL) ? (&DefaultContainer) : container;
  p_container->AddHistogram(this);

  if (p_Histogram != NULL)
    delete[] p_Histogram;
  if (p_Fills != NULL)
    delete[] p_Fills;
  if (p_Errs != NULL)
    delete[] p_Errs;
  if (m_name != NULL)
    delete[] m_name;

  n_entries = 0;
  m_offleft = 0.0;
  m_offright = 0.0;
  m_start = left;
  m_end = right;
  n_bins = bins;
  assert(n_bins > 0);
  p_Histogram = new double[n_bins];
  p_Fills = new int[n_bins];
  p_Errs = new double[n_bins];
  assert(p_Histogram != NULL && p_Fills != NULL);
  for (int i = 0; i < n_bins; i++) {
    p_Histogram[i] = 0.0;
    p_Fills[i] = 0;
    p_Errs[i] = 0.0;
  }

  if (name != NULL) {
    m_name = new char[strlen(name) + 1];
    assert(m_name != NULL);
    strcpy(m_name, name);
  } else
    m_name = NULL;

  m_hist = p_container->GetUniqueId(id);

  if (m_hist == 0) {
    char errmsg[400];
    sprintf(errmsg, "CollieHistogram : Histogram %d already exists, cannot allocate %s", id,
            m_name);
    fprintf(stderr, "%s\n", errmsg);
  }
  assert(m_hist != 0);

  return true;
}

int CollieHistogram::Book(const char* name, int bins, double left, double right,
                          CollieHistogramContainer* container) {
  return Book(0, name, bins, left, right, container);
}

//////////////////////////////////////////////////////////////////////
// Operations
//////////////////////////////////////////////////////////////////////

// definition : ranges are [a,b)

int CollieHistogram::FindBin(double x) {
  if (p_Histogram == NULL)
    return -1;

  if (x < m_start)
    return -1;
  if (x > m_end)
    return n_bins;

  int bin = int((x - m_start) / (m_end - m_start) * double(n_bins));
  if (bin == n_bins)
    bin--;
  return bin;
}

void CollieHistogram::Fill(double x, double w) {
  if (p_Histogram == NULL)
    return;

  n_entries++;

  if (x < m_start) {
    m_offleft += w;
    return;
  }
  // may add flag for this...  p_Histogram[0]++;
  if (x > m_end) {
    m_offright += w;
    return;
  }
  // may add flag for this... p_Histogram[n_bins-1]++;

  int bin = int((x - m_start) / (m_end - m_start) * double(n_bins));
  if (bin == n_bins)
    bin--;

  p_Histogram[bin] += w;
  p_Fills[bin]++;
  return;
}

void CollieHistogram::SetBinError(int i, double w) {
  if (p_Histogram == NULL)
    return;
  if (w < 0)
    return;
  if (i < 0 || i >= n_bins)
    return;
  p_Errs[i] = w;
  return;
}

void CollieHistogram::SetBinContent(int i, double w) {
  if (p_Histogram == NULL)
    return;

  if (i < 0 || i >= n_bins)
    return;
  p_Histogram[i] = w;
  p_Fills[i] = 1;
  return;
}

void CollieHistogram::FillBin(int i, double w) {
  if (p_Histogram == NULL)
    return;
  ;

  n_entries++;

  if (i < 0) {
    m_offleft += w;
    return;
  }
  // may add flag for this...  p_Histogram[0]++;
  if (i >= n_bins) {
    m_offright += w;
    return;
  }
  // may add flag for this... p_Histogram[n_bins-1]++;

  p_Histogram[i] += w;
  p_Fills[i]++; // sqrt(n)
}

void CollieHistogram::Scale(double s) {
  if (p_Histogram == NULL)
    return;
  for (int i = 0; i < n_bins; i++)
    p_Histogram[i] *= s;
  //  if (p_Fills!=NULL)
  //    for (int i=0; i<n_bins; i++)
  //      p_Fills[i]*=s;
  m_offleft *= s;
  m_offright *= s;
}
void CollieHistogram::ScaleBin(int bin, double s) {
  if (p_Histogram == NULL || bin < -1 || bin > n_bins)
    return;
  if (bin == -1)
    m_offleft *= s;
  else if (bin == n_bins)
    m_offright *= s;
  else
    p_Histogram[bin] *= s;
}

void CollieHistogram::Clear() {
  if (p_Histogram == NULL)
    return;
  for (int i = 0; i < n_bins; i++)
    p_Histogram[i] = 0;
  if (p_Fills != NULL)
    for (int i = 0; i < n_bins; i++)
      p_Fills[i] = 0;
  m_offleft = 0;
  m_offright = 0;
  n_entries = 0;
}

void CollieHistogram::Add(const CollieHistogram& h) {
  if (p_Histogram == NULL)
    return;
  ;
  if (h.n_bins != n_bins)
    printf("%d != %d!\n", h.n_bins, n_bins);
  if (h.m_start != m_start)
    printf("%f != %f!\n", h.m_start, m_start);
  if (h.m_end != m_end)
    printf("%f != %f!\n", h.m_end, m_end);
  if (!(p_Histogram != NULL && h.n_bins == n_bins && h.m_start == m_start && h.m_end == m_end))
    return;
  assert(p_Histogram != NULL && h.n_bins == n_bins && h.m_start == m_start && h.m_end == m_end);
  for (int i = 0; i < n_bins; i++) {
    double perr = p_Histogram[i] * p_Errs[i] * p_Errs[i];
    perr += h.p_Histogram[i] * h.p_Errs[i] * h.p_Errs[i];
    double tot = p_Histogram[i] + h.p_Histogram[i];
    if (tot > 0)
      perr /= tot;
    p_Errs[i] = sqrt(perr);
    p_Histogram[i] += h.p_Histogram[i];
    p_Fills[i] += h.p_Fills[i];
  }
  m_offleft += h.m_offleft;
  m_offright += h.m_offright;
  n_entries += h.n_entries;
}

void CollieHistogram::Subtract(const CollieHistogram& h) {
  if (p_Histogram == NULL)
    return;
  ;
  assert(p_Histogram != NULL && h.n_bins == n_bins && h.m_start == m_start && h.m_end == m_end);
  for (int i = 0; i < n_bins; i++) {
    p_Histogram[i] -= h.p_Histogram[i];
    p_Fills[i] = p_Fills[i] + h.p_Fills[i];
  }
  m_offleft -= h.m_offleft;
  m_offright -= h.m_offright;
  n_entries -= h.n_entries;
}

/// Returns a Poisson chi2 compared to the reference histogram
//  Return value is -2ln(L(r;r)/L(r;t)) where
//  L(a;b) is the Poisson likelihood for a events with a mean
//  of b events, were r=reference points and t=test points.
//  Treat this histogram as data
double CollieHistogram::Chi2Test(CollieHistogram* testHist) {
  if (!Valid())
    return -1;
  if (n_bins != testHist->nBins()) {
    printf("CollieHistogram::Chi2Test, histogram binning doesn't match!\n");
    return -1;
  }
  if ((m_start != testHist->MinEdge()) || (m_end != testHist->MaxEdge())) {
    printf("CollieHistogram::Chi2Test, histogram ranges don't match!\n");
    return -1;
  }

  double llrtotal = 0;
  double test;
  double ref;
  double A = 0;
  double B = 0;

  //Construct 2*sum_bins[ t - r - r*ln(t/r) ]
  for (int i = 0; i < n_bins; i++) {
    test = testHist->InBin(i);
    ref = p_Histogram[i];
    if (test == 0 && ref == 0)
      continue;

    A = test - ref;
    B = 0;
    if (ref > 0 && test > 0)
      B = ref * log(test / ref);

    llrtotal += (A - B);
  }
  return 2.0 * llrtotal;
}

/*
 * based on Interpolation routine from Alex Read.  Converted to C++ by Jeremiah Mans.

*.=========================================================================
      SUBROUTINE D_PVMORPH(NB1,XMIN1,XMAX1,DIST1
     &                  ,NB2,XMIN2,XMAX2,DIST2
     &                  ,NBN,XMINN,XMAXN,DISTN
     &                  ,PAR1,PAR2,PARN)
*--------------------------------------------------------------------------
* Author           : Alex Read 11-FEB-1998
*/
CollieHistogram& CollieHistogram::Interpolate(const CollieHistogram& h1, const double p1,
                                              const CollieHistogram& h2, const double p2,
                                              const double pf) {
  const int nybins = 100000; // number of "ybins"  increase this number if spiky-ness occurs
  const int nb_max = 1000;

  int ix, iy, iy0, iy1 = 0;
  //  int ix0,ix1;
  double x0, x, x1, y, y0, y1, yold, dx1, dx2, dxn;
  double dyinv;
  //  double x1old,y1old;
  double slope, b;
  double* sigdis[3];
  double* ydis[3];
  double wta, wtb;
  double inth1, inth2;

  bool allzero;
  //
  //......Debugging on=1.
  //
  const int idebug = 0;
  Clear();
  //
  //......Check that input histograms are not too large.
  //
  if (h1.n_bins > nb_max || h2.n_bins > nb_max) {

    fprintf(stderr, "%s\n", "CollieHistogram::Interpolate :  Too many bins in input histograms.");
    n_bins = -1;
    return *this;
  }
  //
  // From Tom Junk:
  // special case -- zero histogram contents -- sometimes comes up in 2D
  // interpolations
  //
  allzero = true;
  for (ix = 0; ix < h1.n_bins && allzero; ix++)
    if (h1.p_Histogram[ix] != 0.0)
      allzero = false;
  for (ix = 0; ix < h2.n_bins && allzero; ix++)
    if (h2.p_Histogram[ix] != 0.0)
      allzero = false;
  //
  if (allzero) {
    n_bins = MAX(h1.n_bins, h2.n_bins);
    m_start = MIN(h1.m_start, h2.m_start);
    m_end = MAX(h1.m_end, h2.m_end);
    if (p_Histogram != NULL)
      delete[] p_Histogram;
    if (p_Fills != NULL)
      delete[] p_Fills;
    p_Histogram = new double[n_bins];
    p_Fills = new int[n_bins];
    Clear();
    return *this;
  }

  //
  //.....The weights (a,b) are the "distances" between the values of the
  //      parameters at the histograms and the desired interpolation point.
  //      Check that they make sense. If par1=par2 then we can choose any
  //      valid set of wta,wtb so why not take the average?
  //
  if (p2 != p1) {
    wta = 1.0 - (pf - p1) / (p2 - p1);
    wtb = 1.0 + (pf - p2) / (p2 - p1);
  } else {
    wta = 0.5;
    wtb = 0.5;
  }
  //  if (wta<0  || wta>1.0 || wtb<0 || wtb > 1.0 ||
  //    fabs(1.0-(wta+wtb))>1.0e-4) {
  //    fprintf(stderr,"pvmorph.cc :: Invalid weights: %f + %f = %f\n",wta,wtb, wta+wtb);
  //    return;
  //  }
  //
  //......Determine the resulting binning constants.
  //
  if (MIN(wta, wtb) > 0.0) {
    if (h1.m_start == h2.m_start)
      m_start = h1.m_start;
    else
      m_start = wta * h1.m_start + wtb * h2.m_start;
    if (h1.m_end == h2.m_end)
      m_end = h1.m_end;
    else
      m_end = wta * h1.m_end + wtb * h2.m_end;
    if (h1.n_bins == h2.n_bins)
      n_bins = h1.n_bins;
    else
      n_bins = int(wta * h1.n_bins + wtb * h2.n_bins);
  } else {
    m_start = MIN(h1.m_start, h2.m_start);
    m_end = MAX(h1.m_end, h2.m_end);
    n_bins = MAX(h1.n_bins, h2.n_bins);
  }
  dxn = (m_end - m_start) / double(n_bins);

  //
  //......Determine the sizes of the bins used in the
  //      probability inversions.
  //
  dyinv = 1.0 / double(nybins);
  dx1 = (h1.m_end - h1.m_start) / double(h1.n_bins);
  dx2 = (h2.m_end - h2.m_start) / double(h2.n_bins);
  //  if(idebug>=1) {
  //    printf("pvmorph.cc: (debug) %d %f %f %d \n",nybins,dx1,wta,h1.n_bins);
  //    printf("pvmorph.cc: (debug) %d %f %f %d \n",nybins,dx2,wtb,h2.n_bins);
  //  }
  //
  //......Copy the two histograms to temporary locations and
  //      overwrite the contents with the normalized, cummulative
  //      probability distribution.
  //
  // JMM -> allocate memory
  sigdis[0] = new double[h1.n_bins + 10];
  sigdis[1] = new double[h2.n_bins + 10];
  sigdis[2] = new double[n_bins + 10];
  ydis[0] = new double[nybins + 10];
  ydis[1] = new double[nybins + 10];
  ydis[2] = new double[nybins + 10];

  sigdis[0][0] = h1.p_Histogram[0];
  for (ix = 1; ix < h1.n_bins; ix++)
    sigdis[0][ix] = h1.p_Histogram[ix] + sigdis[0][ix - 1];
  inth1 = sigdis[0][h1.n_bins - 1];
  for (ix = 0; ix < h1.n_bins; ix++)
    sigdis[0][ix] /= inth1;

  sigdis[1][0] = h2.p_Histogram[0];
  for (ix = 1; ix < h2.n_bins; ix++)
    sigdis[1][ix] = h2.p_Histogram[ix] + sigdis[1][ix - 1];
  inth2 = sigdis[1][h2.n_bins - 1];
  for (ix = 0; ix < h2.n_bins; ix++)
    sigdis[1][ix] /= inth2;

  if (idebug >= 1) {
    fprintf(stderr, "%s\n",
            "CollieHistogram::Interpolate(): (debug) temp cumm. dists. computed.\n");
  }

  //
  //......Invert the cummulative distributions by scanning in the y-direction
  //      and looking for the x that corresponds to this y. The complicated
  //      business is due to the digitization of the probability. Ugh.
  //
  yold = 0.0;
  for (ix = 1; ix <= h1.n_bins && iy1 < nybins; ix++) {
    //    x0 = (ix-1)*dx1 + h1.m_start;
    x0 = (ix - 1) * dx1 + h1.m_start;
    x1 = x0 + dx1;
    y0 = yold;
    y1 = sigdis[0][ix - 1];
    if (y1 >= (1.0 - dyinv))
      y1 = 1.0;
    slope = (y1 - y0) / (x1 - x0);
    b = y1 - slope * x1;
    iy0 = int(y0 * double(nybins)); // /dyinv);
    iy1 = int(y1 * double(nybins)); ///dyinv);
    if (idebug >= 2)
      printf("pvmorph.cc: (debug) %d %f %f %f %f %d %d\n", ix, x0, x1, y0, y1, iy0 + 1, iy1 + 1);
    if (iy1 > 0) {
      if (idebug >= 2)
        printf("pvmorph.cc: (debug) iy0=%d,iy1=%d,x0=%f,x1=%f\n", iy0, iy1, x0, x1);
      for (iy = iy0; iy <= iy1; iy++) {
        if (y0 < (1.0 - dyinv) && slope != 0.0) {
          //                  y = iy*dyinv
          y = iy * dyinv;
          //          ydis(iy+1,1) = (y-b)/slope
          ydis[0][iy + 1] = (y - b) / slope;
        } else if (y0 >= (1.0 - dyinv)) {
          //                  ydis(iy+1,1) = ydis(iy0-1,1)
          ydis[0][iy + 1] = ydis[0][iy0 - 1];
        }
      }
    }

    yold = y1;
  }
  //
  yold = 0.0;
  iy1 = 0;
  for (ix = 1; ix <= h2.n_bins && iy1 < nybins; ix++) {
    x0 = (ix - 1) * dx2 + h2.m_start;
    x1 = x0 + dx2;
    y0 = yold;
    y1 = sigdis[1][ix - 1];
    if (y1 >= (1.0 - dyinv))
      y1 = 1.0;
    slope = (y1 - y0) / (x1 - x0);
    b = y1 - slope * x1;
    iy0 = int(y0 * double(nybins)); ///dyinv);
    iy1 = int(y1 * double(nybins)); ///dyinv);
    if (idebug >= 2)
      printf("pvmorph.cc: (debug) iy0=%d,iy1=%d,x0=%f,x1=%f\n", iy0, iy1, x0, x1);
    if (iy1 > 0) {
      for (iy = iy0; iy <= iy1; iy++) {
        if (y0 < (1.0 - dyinv) && slope != 0.0) {
          y = iy * dyinv;
          //        ydis[1][iy+1] = (y-b)/slope;
          ydis[1][iy + 1] = (y - b) / slope;
        } else if (y0 >= (1.0 - dyinv)) {
          ydis[1][iy + 1] = ydis[1][iy0 - 1];
        }
      }
    }
    yold = y1;
  }
  //
  //......Compute interpolated mean x(y)
  //
  for (iy = 1; iy <= nybins; iy++) {
    ydis[2][iy] = wta * ydis[0][iy] + wtb * ydis[1][iy];
    //    ydis[2][iy]=(ydis[0][iy]-ydis[1][iy])/(p1-p2)*pf+(p2*ydis[0][iy]/(p2-p1)-p1*ydis[1][iy]/(p2-p1));

    if (idebug >= 2 && (iy < 12 || iy >= nybins - 12))
      printf("pvmorph.cc: (debug) %d, %f, %f %f\n", iy, ydis[0][iy], ydis[1][iy], ydis[2][iy]);
  }
  if (idebug >= 1)
    printf("pvmorph.cc: (debug) solving x_interp(y) complete.\n");

  //
  //......Now find y(x) for the x's at the upper edges of the new
  //      coarse-grained histogram. This looks much more complicated
  //      than one might expect because (1) extrapolation needs
  //      special treatment (one of the weights is negative) and
  //      (2) we have to insure that users who send two identical
  //      histogram binnings in get out EXACTLY the same binning.

  yold = 0.0;
  iy = 1;
  for (ix = 1; ix <= n_bins; ix++) {
    x = ix * dxn + m_start;
    //         IF(idebug.GE.2) PRINT *,'PVMORPH ',ix0,ix1,ix,x1old,y1old
    while (ydis[2][iy] < x && iy <= nybins)
      iy = iy + 1;
    y = (iy - 1) * dyinv;
    sigdis[2][ix - 1] = y;
  }
  if (idebug >= 1)
    printf("pvmorph.cc: (debug) solving y(x_interp) complete.\n");

  //
  //.....Differentiate to get distribution from cummulative distribution.
  //
  for (ix = n_bins - 1; ix > 0; ix--)
    sigdis[2][ix] = sigdis[2][ix] - sigdis[2][ix - 1];
  //
  //......Copy the result to the output distribution.
  //
  if (p_Histogram != NULL)
    delete[] p_Histogram;
  if (p_Fills != NULL)
    delete[] p_Fills;
  p_Histogram = new double[n_bins];
  p_Fills = new int[n_bins];

  double scale;

  if (p1 == p2)
    scale = inth1;
  else {
    //    scale=(inth1-inth2)/(p1-p2)*pf+(p2*inth1/(p2-p1)-p1*inth2/(p2-p1));
    scale = (pf - p1) * (inth2 - inth1) / (p2 - p1) + inth1;
  }

  for (int q = 0; q < n_bins; q++) {
    p_Histogram[q] = sigdis[2][q] * scale;
    //    p_Fills[q]=(int)sqrt(ABS(p_Histogram[q]));
    p_Fills[q] = 1;
    m_offleft = (h2.m_offleft - h1.m_offleft) / (p2 - p1) * (pf - p1) + h1.m_offleft;
    m_offright = (h2.m_offright - h1.m_offright) / (p2 - p1) * (pf - p1) + h1.m_offright;
  }

  for (int q = 0; q < 3; q++) {
    if (sigdis[q])
      delete[] (sigdis[q]);
    if (ydis[q])
      delete[] (ydis[q]);
  }

  return *this;
}

//////////////////////////////////////////////////////////////////////
// Information Methods
//////////////////////////////////////////////////////////////////////

double CollieHistogram::Integral() const {
  return Integral(-1, n_bins);
}
double CollieHistogram::Integral(int from, int to) const {
  double retval = 0.0;
  if (p_Histogram == NULL)
    return 0.0;
  for (int i = MAX(0, from); i <= (MIN(n_bins - 1, to)); i++)
    retval += p_Histogram[i];

  if (from < 0)
    retval += m_offleft;
  if (to == n_bins)
    retval += m_offright;

  return retval;
}
double CollieHistogram::MathInt() const {
  return MathInt(-1, n_bins);
}
double CollieHistogram::MathInt(int from, int to) const {
  return Integral(from, to) * (m_end - m_start) / double(n_bins);
}
double* CollieHistogram::CopyOf() const {
  if (p_Histogram == NULL)
    return NULL;
  double* copyd = new double[n_bins];
  if (copyd == NULL)
    return NULL;
  memcpy(copyd, p_Histogram, sizeof(double) * n_bins);
  return copyd;
}

double CollieHistogram::InBin(int nbin) const {
  if (p_Histogram == NULL || nbin < 0 || nbin >= n_bins) {
    char s[1024];
    sprintf(s, "CollieHistogram::InBin(%d) invalid access! (name='%s',max=%d)", nbin, m_name,
            n_bins);
    fprintf(stderr, "%s\n", s);
    return 0.0;
  }
  return p_Histogram[nbin];
}

double CollieHistogram::BinErr(int nbin) const {
  if (p_Histogram == NULL || nbin < 0 || nbin >= n_bins) {
    char s[1024];
    sprintf(s, "CollieHistogram::BinErr(%d) invalid access! (name='%s',max=%d)", nbin, m_name,
            n_bins);
    fprintf(stderr, "%s\n", s);
    return 0.0;
  }
  return p_Errs[nbin];
}

int CollieHistogram::FillsInBin(int nbin) const {
  if (p_Histogram == NULL || nbin < 0 || nbin >= n_bins) {
    char s[1024];
    sprintf(s, "CollieHistogram::FillsInBin(%d) invalid access! (name='%s',max=%d)", nbin, m_name,
            n_bins);
    fprintf(stderr, "%s\n", s);
    return 0;
  }
  return p_Fills[nbin];
}

double CollieHistogram::ValueIn(double v) const {
  if (p_Histogram == NULL) {
    fprintf(stderr, "%s\n", "CollieHistogram::ValueIn() invalid access!");
    return 0.0;
  }
  if (v < m_start)
    return m_offleft;
  if (v >= m_end)
    return m_offright;

  int bin = int((v - m_start) / (m_end - m_start) * double(n_bins));
  if (bin == n_bins)
    bin--;

  return p_Histogram[bin];
}

double CollieHistogram::GetMaximum() const {
  if (p_Histogram == NULL)
    return 0.0;
  double max = 0.0;
  for (int i = 0; i < n_bins; i++)
    if (p_Histogram[i] > max)
      max = p_Histogram[i];

  return max;
}

double CollieHistogram::Mean() const {
  double retval = 0.0;
  if (p_Histogram == NULL)
    return 0.0;
  double x;
  double sum = 0.0;
  double mean = 0.0;
  for (int i = 0; i < n_bins; i++) {
    x = (double(i) + 0.5) * (m_end - m_start) / double(n_bins) + m_start;
    sum += p_Histogram[i];
    mean += x * p_Histogram[i];
  }
  if (sum > 0.0)
    retval = mean / sum;
  return retval;
}

double CollieHistogram::RMS() const {
  double retval = 0.0;
  if (p_Histogram == NULL)
    return 0.0;
  double mean, x;
  double sum = 0.0;
  double rms = 0.0;
  mean = Mean();
  for (int i = 0; i < n_bins; i++) {
    x = (double(i) + 0.5) * (m_end - m_start) / double(n_bins) + m_start;
    sum += p_Histogram[i];
    rms += (x - mean) * (x - mean) * p_Histogram[i];
  }
  if (sum > 0.0 && rms > 0.0) {
    retval = sqrt(rms / sum);
  } else
    retval = 0.0;
  return retval;
}

//////////////////////////////////////////////////////////////////////
// Storage Methods
//////////////////////////////////////////////////////////////////////
void CollieHistogram::StoreHistograms(const char* filename, const char* filetype,
                                      CollieHistogramContainer* container) {
  if (container == NULL)
    DefaultContainer.StoreHistograms(filename, filetype);
  else
    container->StoreHistograms(filename, filetype);
}

void CollieHistogram::StoreHistograms(const char* filename, const char* filetype,
                                      vector<TH1*> externHists) {

  for (vector<TH1*>::iterator ir = externHists.begin(); ir != externHists.end(); ir++) {
    if ((*ir) != NULL) {
      if ((*ir)->GetNbinsY() <= 1) {
        CollieHistogram* inHist = new CollieHistogram();
        inHist->BookROOT(*ir);
        delete inHist;
      } else {
        CollieHistogram2d* inHist = new CollieHistogram2d();
        inHist->BookROOT(*ir);
        delete inHist;
      }
    }
  }
  DefaultContainer.StoreHistograms(filename, filetype);
}

void CollieHistogram::ClearHistograms(CollieHistogramContainer* container) {
  p_container = (container == NULL) ? (&DefaultContainer) : container;
  if (p_container != NULL)
    p_container->ClearHistograms();
}

//////////////////////////////////////////////////////////////////////
// Operator
////////////////////////////////////////////////////////////////////////

CollieHistogram& CollieHistogram::operator=(const CollieHistogram& h) {
  if (this == &h)
    return *this;
  n_bins = h.n_bins;
  n_entries = h.n_entries;
  m_start = h.m_start;
  m_end = h.m_end;
  m_offleft = h.m_offleft;
  m_offright = h.m_offright;

  if (p_Histogram != NULL)
    delete[] p_Histogram;
  if (p_Fills != NULL)
    delete[] p_Fills;
  if (m_name != NULL)
    delete[] m_name;
  if (p_container != NULL)
    p_container->RemoveHistogram(this);
  p_container = h.p_container;
  p_container->AddHistogram(this);

  if (h.p_Histogram == NULL)
    p_Histogram = NULL;
  else {
    p_Histogram = new double[n_bins];
    assert(p_Histogram != NULL);
    for (int i = 0; i < n_bins; i++)
      p_Histogram[i] = h.p_Histogram[i];
  }
  if (h.p_Fills != NULL) {
    p_Fills = new int[n_bins];
    assert(p_Fills != NULL);
    for (int i = 0; i < n_bins; i++)
      p_Fills[i] = h.p_Fills[i];
  } else
    p_Fills = NULL;

  if (h.m_name != NULL) {
    m_name = new char[strlen(h.m_name) + 1];
    assert(m_name != NULL);
    strcpy(m_name, h.m_name);
  } else
    m_name = NULL;

  return *this;
}

CollieHistogram& CollieHistogram::operator+=(const CollieHistogram& h) {
  Add(h);
  return *this;
}
CollieHistogram& CollieHistogram::operator-=(const CollieHistogram& h) {
  Subtract(h);
  return *this;
}
CollieHistogram& CollieHistogram::operator*=(double s) {
  Scale(s);
  return *this;
}

//class  : public CollieHistogram
//{
//protected:
//  int n_x,n_y;
//  double m_lowx,m_lowy,m_highx,m_highy;
//public:
// constructor/book
////////////////////////////////////////////

CollieHistogram2d::CollieHistogram2d()
    : CollieHistogram() {
  n_x = 0;
  n_y = 0;
}

CollieHistogram2d::CollieHistogram2d(const CollieHistogram2d& h) {
  (*this) = h;
}

CollieHistogram2d::~CollieHistogram2d() {}

int CollieHistogram2d::Book(int id, const char* name, int xbins, int ybins, double xmin,
                            double xmax, double ymin, double ymax,
                            CollieHistogramContainer* container) {

  CollieHistogram::Book(id, name, xbins * ybins, 0, 1, container);

  n_x = xbins;
  n_y = ybins;

  m_lowx = xmin;
  m_highx = xmax;
  m_lowy = ymin;
  m_highy = ymax;

  for (int i = 0; i < 8; i++)
    off[i] = 0.0;

  assert(m_hist != 0);

  return true;
}

int CollieHistogram2d::Book(const char* name, int xbins, int ybins, double xmin, double xmax,
                            double ymin, double ymax, CollieHistogramContainer* container) {
  return Book(0, name, xbins, ybins, xmin, xmax, ymin, ymax, container);
}
// operations
//////////////////////////////////
void CollieHistogram2d::Fill(double x, double y, double w) {
  if (p_Histogram == NULL)
    return;

  int i, j;
  i = int((x - m_lowx) / (m_highx - m_lowx) * double(n_x));

  j = int((y - m_lowy) / (m_highy - m_lowy) * double(n_y));

  Fill(i, j, w);
}

void CollieHistogram2d::Fill(int i, int j, double w) {
  if (p_Histogram == NULL)
    return;

  n_entries++;

  if (i < 0) {
    if (j < 0)
      off[SOUTHWEST2D] += w;
    if (j >= n_y)
      off[NORTHWEST2D] += w;
    else
      off[WEST2D] += w;
    return;
  }
  if (i >= n_x) {
    if (j < 0)
      off[SOUTHEAST2D] += w;
    if (j >= n_y)
      off[NORTHEAST2D] += w;
    else
      off[EAST2D] += w;
    return;
  }

  if (j < 0 || j >= n_y) {
    if (j < 0)
      off[SOUTH2D] += w;
    else
      off[NORTH2D] += w;
    return;
  }

  p_Histogram[i * n_y + j] += w;
  p_Fills[i * n_y + j]++;
}

void CollieHistogram2d::SetBinError(int i, int j, double w) {
  if (p_Histogram == NULL)
    return;
  if (w < 0)
    return;
  if (i < 0 || i >= n_x)
    return;
  if (j < 0 || j >= n_y)
    return;
  p_Errs[i * n_y + j] = w;
  return;
}

void CollieHistogram2d::SetBinContent(int i, int j, double w) {
  if (p_Histogram == NULL)
    return;
  if (i < 0 || i >= n_x)
    return;
  if (j < 0 || j >= n_y)
    return;
  p_Histogram[i * n_y + j] = w;
  p_Fills[i * n_y + j] = 1;
  return;
}

void CollieHistogram2d::Add(const CollieHistogram2d& h) {
  if (p_Histogram == NULL)
    return;
  assert(!(n_x == 0 || n_x != h.n_x || n_y != h.n_y || m_lowx != h.m_lowx || m_lowy != h.m_lowy));
  for (int i = 0; i < n_x * n_y; i++) {
    double perr = p_Histogram[i] * p_Errs[i] * p_Errs[i];
    perr += h.p_Histogram[i] * h.p_Errs[i] * h.p_Errs[i];
    double tot = p_Histogram[i] + h.p_Histogram[i];
    if (tot > 0)
      perr /= tot;
    p_Errs[i] = sqrt(perr);
    p_Histogram[i] += h.p_Histogram[i];
    p_Fills[i] += h.p_Fills[i];
  }
  for (int i = 0; i < 8; i++)
    off[i] += h.off[i];
}

void CollieHistogram2d::Subtract(const CollieHistogram2d& h) {
  if (p_Histogram == NULL)
    return;
  ;
  assert(!(n_x == 0 || n_x != h.n_x || n_y != h.n_y || m_lowx != h.m_lowx || m_lowy != h.m_lowy));
  for (int i = 0; i < n_x * n_y; i++) {
    p_Histogram[i] -= h.p_Histogram[i];
    p_Fills[i] += h.p_Fills[i];
  }
  for (int i = 0; i < 8; i++)
    off[i] -= h.off[i];
}

//////////////////////////////////////////////////////////////////////
// Information Methods
//////////////////////////////////////////////////////////////////////

double CollieHistogram2d::Integral() const {
  return CollieHistogram2d::Integral(-1, n_x, -1, n_y);
}
double CollieHistogram2d::Integral(int fromX, int toX, int fromY, int toY) const {
  if (p_Histogram == NULL)
    return 0.0;

  double retval = 0.0;

  for (int x = MAX(0, fromX); x <= (MIN(n_x - 1, toX)); x++) {
    for (int y = MAX(0, fromY); y <= (MIN(n_y - 1, toY)); y++) {
      retval += p_Histogram[x * n_y + y];
    }
  }

  if (fromX < 0)
    retval += off[WEST2D];
  if (fromY < 0)
    retval += off[SOUTH2D];
  if (fromX < 0 && fromY < 0)
    retval += off[SOUTHWEST2D];

  if (toX == n_x)
    retval += off[EAST2D];
  if (toY == n_y)
    retval += off[NORTH2D];
  if (toX == n_x && toY == n_y)
    retval += off[NORTHEAST2D];

  if (fromX < 0 && toY == n_y)
    retval += off[NORTHWEST2D];
  if (fromY < 0 && toX == n_x)
    retval += off[SOUTHEAST2D];

  return retval;
}

/*
double* CollieHistogram2d::CopyOf() const {
  if (p_Histogram==NULL) return NULL;
  double* copyd=new double[n_bins];
  if (copyd==NULL) return NULL;
  memcpy(copyd,p_Histogram,sizeof(double)*n_bins);
  return copyd;
}
*/
double CollieHistogram2d::InBin(int nbinx, int nbiny) const {
  if (p_Histogram == NULL || nbinx < 0 || nbinx >= n_x || nbiny < 0 || nbiny >= n_y) {
    char s[1024];
    sprintf(s, "CollieHistogram2d::InBin(%d,%d) invalid access! (name='%s',max=%d,%d)", nbinx,
            nbiny, m_name, n_x, n_y);
    fprintf(stderr, "%s\n", s);
    return 0.0;
  }
  int loc = nbinx * n_y + nbiny;
  return p_Histogram[loc];
}

double CollieHistogram2d::BinErr(int nbinx, int nbiny) const {
  if (p_Histogram == NULL || nbinx < 0 || nbinx >= n_x || nbiny < 0 || nbiny >= n_y) {
    char s[1024];
    sprintf(s, "CollieHistogram2d::InBin(%d,%d) invalid access! (name='%s',max=%d,%d)", nbinx,
            nbiny, m_name, n_x, n_y);
    fprintf(stderr, "%s\n", s);
    return 0.0;
  }
  int loc = nbinx * n_y + nbiny;
  return p_Errs[loc];
}

int CollieHistogram2d::FillsInBin(int nbinx, int nbiny) const {
  if (p_Histogram == NULL || nbinx < 0 || nbinx >= n_x || nbiny < 0 || nbiny >= n_y) {
    char s[1024];
    sprintf(s, "CollieHistogram2d::FillsInBin(%d,%d) invalid access! (name='%s',max=%d,%d)", nbinx,
            nbiny, m_name, n_x, n_y);
    fprintf(stderr, "%s\n", s);
    return 0;
  }
  int loc = nbinx * n_y + nbiny;
  return p_Fills[loc];
}
/*
double CollieHistogram2d::ValueIn(double v) const {
  if (p_Histogram==NULL) {
    fprintf(stderr,"%s\n","CollieHistogram2d::ValueIn() invalid access!");
    return 0.0;
  }
  if (v<m_start) return m_offleft;
  if (v>=m_end) return m_offright;

  int bin=int((v-m_start)/(m_end-m_start)*double(n_bins));
  if (bin==n_bins) bin--;

  return p_Histogram[bin];
}

double CollieHistogram2d::Mean() const {
  double retval=0.0;
  if (p_Histogram==NULL) return 0.0;
  double x;
  double sum=0.0;
  double mean=0.0;
  for (int i=0; i<n_bins; i++) {
    x = (double(i)+0.5)*(m_end-m_start)/double(n_bins)+m_start;
    sum+=p_Histogram[i];
    mean+=x*p_Histogram[i];
  }
  if(sum>0.0) retval=mean/sum;
  return retval;
}

double CollieHistogram2d::RMS() const {
  double retval=0.0;
  if (p_Histogram==NULL) return 0.0;
  double mean,x;
  double sum=0.0;
  double rms=0.0;
  mean = Mean();
  for (int i=0; i<n_bins; i++) {
    x = (double(i)+0.5)*(m_end-m_start)/double(n_bins)+m_start;
    sum+=p_Histogram[i];
    rms+=(x-mean)*(x-mean)*p_Histogram[i];
  }
  if(sum>0.0 && rms>0.0) {
    retval=sqrt(rms/sum);
  } else retval=0.0;
  return retval;
}
*/
// storage methods
//////////////////////////////////

CollieHistogram2d& CollieHistogram2d::operator=(const CollieHistogram2d& h) {
  CollieHistogram(*this) = (CollieHistogram)h;
  n_x = h.n_x;
  n_y = h.n_y;
  m_lowx = h.m_lowx;
  m_highx = h.m_highx;
  m_lowy = h.m_lowy;
  m_highy = h.m_highy;
  for (int i = 0; i < 8; i++)
    off[i] = h.off[i];
  return *this;
}
