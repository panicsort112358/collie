/* -*- C++ -*- */
#ifndef CollieHISTOGRAM_H_INCLUDED
#define CollieHISTOGRAM_H_INCLUDED

class CollieHistogram;

#include <vector>
#include <map>
#include <string>
#include <math.h>
#include <algorithm>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"

using namespace std;

/** \brief Container for groups of histograms
 ** \ingroup IO
 **/
class CollieHistogramContainer {
public:
  /// add a histogram to this container
  bool AddHistogram(CollieHistogram* p_hist);
  bool AddROOT(TH1* p_hist);
  void ClearHistograms();

  /// remove a histogram from this container
  bool RemoveHistogram(CollieHistogram* p_hist);
  /// get a unique histogram id for this container
  int GetUniqueId(int id);
  /// store the histograms in this container to a file
  void StoreHistograms(const char* filename, const char* filetype);

private:
  std::vector<TH1*> p_extern;

  std::vector<CollieHistogram*> p_histograms;
};

/** \brief Histogram storage interface
    This interface can be used to define storage methods for histograms.
    As shipped, Histograms can be stored in JArchive files, PAW HBOOK files,
    or ROOT files.  Additional storage mechanisms can be registered by
    calling AddStorer().
*/
class CollieHistogramStorer {
private:
  static std::map<std::string, CollieHistogramStorer*> gl_stores;

public:
  /// Obtain the storer object for a given file type string (NULL gets the default)
  static CollieHistogramStorer* GetStorerFor(const char* filetype);
  /// Add a storer for a given file type
  static void AddStorer(const char* filetype, CollieHistogramStorer* storer);
  /// Store the given histograms into the given file
  virtual void StoreHistograms(const char* filename, vector<CollieHistogram*>& hists) = 0;
  virtual void StoreHistograms(const char* filename, vector<CollieHistogram*> hists,
                               vector<TH1*> rHists) = 0;
};

#include <assert.h>

/** \brief 1D Histogramming class
  \ingroup Core
CollieHistogram is a histogramming class much like HBOOK, but written in C++.  The histograms can be
written to ZEBRA PAW files, JArchive direct-access files, or (soon) ROOT files.
*/
class CollieHistogram {
public:
protected:
  /// pointer to the container of this histogram
  CollieHistogramContainer* p_container;
  /// pointer to the histogram data
  double* p_Histogram;
  /// pointer to the histogram fill count
  int* p_Fills;
  /// pointer to the histogram fill count
  double* p_Errs;
  /// number of bins in the histogram
  int n_bins;
  /// index of this histogram
  int m_hist;
  /// number of entries (calls to fill)
  int n_entries;
  /// lower edge of the histogram
  double m_start;
  /// upper edge of the histogram
  double m_end;
  /// lower overflow bin
  double m_offleft;
  /// upper overflow bin
  double m_offright;
  /// name of the histogram container
  char* m_name;

public:
  /// default container for all histograms which don't get assigned a container
  static CollieHistogramContainer DefaultContainer;
  // constructor/book
  /// default constuctor requires eventual use of Book
  CollieHistogram();
  /// copy constructor
  CollieHistogram(const CollieHistogram& h);
  /// Destructor (frees id)
  virtual ~CollieHistogram();

  /**  Creates a CollieHistogram with \c bins bins, starting at \c left and reaching to
  \c right.  The histogram is booked with the \c name and \c id given in the
  directory \c dir.  If \c dir is NULL, the default directory is used.  If \c id conflicts with a previously booked histogram,
  this condition triggers an assert.
  */
  virtual int Book(int id, const char* name, int bins, double left, double right,
                   CollieHistogramContainer* container = NULL);
  /** This method works exactly as the Book method with specified id, except
      the id is automatically allocated from those not yet in use.
  */
  virtual int Book(const char* name, int bins, double left, double right,
                   CollieHistogramContainer* container = NULL);

  bool BookROOT(TH1* inhist, CollieHistogramContainer* container = NULL);

  // operations
  /// Determines the bin which contains \c x , \c y (-1 if underflow, nbins if overflow)
  int FindBin(double x);
  /// Adds \c w to the bin which contains \c x.
  virtual void Fill(double x, double w = 1.0);
  /// Adds \c w to bin number \c i.
  inline void Fill(int i, double w = 1.0) {
    FillBin(i, w);
  }
  /// Adds \c w to bin number \c i.
  virtual void FillBin(int i, double w = 1.0);
  /// Sets err = \c w for bin number \c i.
  virtual void SetBinError(int i, double err);
  /// Sets err = \c w for bin number \c i.
  virtual void SetBinContent(int i, double val);
  /// Multiplies every bin by \c s.  Alternatively explained, Scale() scales the integrated histogram by \c s.
  virtual void Scale(double s);
  /// Multiplies the given bin by \c s
  virtual void ScaleBin(int bin, double s);
  /// Empties the histogram, but does not unbook it.
  virtual void Clear();
  /** Adds \c h to the current histogram bin by bin.
      If the \c left, \c right, or \c bins parameters are different, Add() dies with
      an assert.
  */
  virtual void Add(const CollieHistogram& h);
  /** Subtracts \c h from the current histogram bin by bin.
      If the \c left, \c right, or \c bins parameters are different, Subtract() dies with
      an assert.
  */
  virtual void Subtract(const CollieHistogram& h);

  inline bool Valid() const {
    return p_Histogram != NULL;
  }

  // information methods
  ///  Returns the number of bins in the histogram
  inline int nBins() const {
    return n_bins;
  }
  /// Returns the lower edge of the histogram
  inline double MinEdge() const {
    return m_start;
  }
  /// Returns the upper edge of the histogram
  inline double MaxEdge() const {
    return m_end;
  }
  /// Returns the amount off the left side of the histogram
  inline double Underflow() const {
    return m_offleft;
  }
  /// Returns the amount off the right side of the histogram
  inline double Overflow() const {
    return m_offright;
  }
  /// Return the title of this Histogram
  inline const char* GetName() const {
    return m_name;
  }
  /** Returns the sum of the bin contents.
      (Total number of events in most cases)  Adds all bins, including the underflow and overflow bins.
  */
  double Integral() const;
  /** Returns the sum of the bin contents.
      (Number of events in most cases)  Includes the
      underflow bin if \c from=-1 and the overflow bin if \c to = \c bins.
  */
  double Integral(int from, int to) const;
  /** works like Integral() except
      it returns the sum of the bin contents times the bin width.  This is the mathematical integration
      of the histogram.
  */
  double MathInt() const;
  /** works like Integral(int from, int to) except
      it returns the sum of the bin contents times the bin width.  This is the mathematical integration
      of the histogram.
  */
  double MathInt(int from, int to) const;
  ///  Returns an array of double containing a copy of the histogram bins.
  double* CopyOf() const;
  /// Returns the value of the n_th bin
  double InBin(int nbin) const;
  /// Returns the error of the n_th bin
  double BinErr(int nbin) const;
  /// Returns the value of the n_th bin
  int FillsInBin(int nbin) const;

  /// Returns the value of the bin containing this point
  double ValueIn(double v) const;
  /// returns the id of this histogram
  inline int Id() const {
    return m_hist;
  }
  /// returns computed mean of this histogram
  double Mean() const;
  /// returns computed rms of this histogram
  double RMS() const;

  /// returns a Poisson chi2 compared to the reference histogram
  double Chi2Test(CollieHistogram* test);

  /// return the maximum value in the histogram
  double GetMaximum() const;

  /// return the bin width
  inline double GetBinWidth() const {
    return (m_end - m_start) / double(n_bins);
  }

  /// return the bin center
  inline double GetBinCenter(int ib) const {
    return m_start + GetBinWidth() * (ib + 0.5);
  }

  /** \brief interpolation method
   **
   ** This method replaces the current histogram with an interpolated histogram between
   ** the two passed histograms with the given interpolation parameters.  These parameters
   ** might be Higgs masses, or just about anything...
   **/
  CollieHistogram& Interpolate(const CollieHistogram& h1, const double p1,
                               const CollieHistogram& h2, const double p2, const double pf);

  // storage methods
  /** \brief Global output method
   **
   ** This method outputs all histograms to the file \c filename of type \c filetype.
   ** Eventually this method will be replaced by a HistogramFileContext but not now.
   **/
  void ClearHistograms(CollieHistogramContainer* container = NULL);
  static void
  StoreHistograms(const char* filename, const char* filetype = NULL,
                  CollieHistogramContainer* container = NULL); // Store all Histograms...
  static void StoreHistograms(const char* filename, const char* filetype,
                              vector<TH1*> externHists); // Store all Histograms...

  // operators
  /// Sets this histogram equal to \c h
  CollieHistogram& operator=(const CollieHistogram& h);
  /// Adds \c h to this histogram. (See Add())
  CollieHistogram& operator+=(const CollieHistogram& h);
  /// Subtracts \c h from this histogram. (See Subtract())
  CollieHistogram& operator-=(const CollieHistogram& h);
  /// Scales this histogram by \c s. (See Scale())
  CollieHistogram& operator*=(double s);
  virtual const char* getType() {
    return "CollieHistogram";
  }
};

/** Two-dimensional histogram class based on CollieHistogram
 ** \ingroup Core
 **/
class CollieHistogram2d : public CollieHistogram {
private:
  int n_x, n_y;
  double m_lowx, m_lowy, m_highx, m_highy;
  double off[8];

protected:
public:
  /// Creates blank 2d Histogram.  Requires use of Book()
  CollieHistogram2d();
  /// Copy constructor
  CollieHistogram2d(const CollieHistogram2d& h);
  /// Destructor
  virtual ~CollieHistogram2d();

  /**  \brief Main 2d Book method
   **
   **  Creates a CollieHistogram with \c xbins across, starting at \c xmin at the left
   **  and reaching to  \c xmax on the right.  The histogram is booked with \c ybins in the
   **  the second dimension, from \c ymin to \c ymax.  The histogram is given the
   **  \c name and \c id passed and stored in the (internal) directory \c dir.
   **  If \c dir is NULL, the default directory is used.  If \c id conflicts
   **  with a previously booked histogram,  this condition triggers an assert.
   */
  int Book(int id, const char* name, int xbins, int ybins, double xmin, double xmax, double ymin,
           double ymax, CollieHistogramContainer* container = NULL);
  /** This method works exactly as the Book method with specified id, except
      the id is automatically allocated from those not yet in use.
  */
  int Book(const char* name, int xbins, int ybins, double xmin, double xmax, double ymin,
           double ymax, CollieHistogramContainer* container = NULL);

  // operations
  /// Adds \c w to the bin which contains \c x , \c y.
  void Fill(double x, double y, double w = 1.0);
  /// Adds \c w to bin  [ \c i , \c j ]
  void Fill(int i, int j, double w = 1.0);
  void SetBinError(int i, int j, double w);
  void SetBinContent(int i, int j, double w);

  /** Adds \c h to the current histogram bin by bin.
      If the \c xmin, \c xmax, \c xbins, \c ymin, \c ymax, or \c ybins parameters are different, Add() dies with
      an assert.
  */
  void Add(const CollieHistogram2d& h);
  /** Subtracts \c h from the current histogram bin by bin.
      If the \c xmin, \c xmax, \c xbins, \c ymin, \c ymax, or \c ybins parameters are different, Subtract() dies with
      an assert.
  */
  void Subtract(const CollieHistogram2d& h);

  // information methods
  /// returns \c xbins
  inline int nxBins() const {
    return n_x;
  }
  /// returns \c ybins
  inline int nyBins() const {
    return n_y;
  }
  /// returns the left edge value
  inline double LeftEdge() const {
    return m_lowx;
  }
  /// returns the right edge value
  inline double RightEdge() const {
    return m_highx;
  }
  /// returns the bottom edge value
  inline double BottomEdge() const {
    return m_lowy;
  }
  /// returns the top edge value
  inline double TopEdge() const {
    return m_highy;
  }
  /// returns an overflow bin
  inline double Overflow(int bin) const {
    return off[bin % 8];
  }
  //  double MinEdge() const { return m_start; }
  //  double MaxEdge() const { return m_end; }
  //  double Underflow() const { return m_offleft; }
  //  double Overflow() const { return m_offright; }

  /// return the bin width
  inline double GetBinWidthX() const {
    return (m_highx - m_lowx) / double(n_x);
  }
  inline double GetBinWidthY() const {
    return (m_highy - m_lowy) / double(n_y);
  }

  /// return the bin center
  inline double GetBinCenterX(int ib) const {
    return m_lowx + GetBinWidthX() * (ib + 0.5);
  }
  inline double GetBinCenterY(int ib) const {
    return m_lowy + GetBinWidthY() * (ib + 0.5);
  }

  double Integral() const;
  double Integral(int fromx, int tox, int fromy, int toy) const;
  //  double* CopyOf() const;
  double InBin(int nbinx, int nbiny) const;
  double BinErr(int nbinx, int nbiny) const;
  int FillsInBin(int nbinx, int nbiny) const;

  // operators
  /// Sets this 2d histogram equal to \c h
  CollieHistogram2d& operator=(const CollieHistogram2d& h);
  /// Adds \c h to this histogram. (See Add())
  CollieHistogram2d& operator+=(const CollieHistogram2d& h) {
    Add(h);
    return *this;
  }
  /// Subtracts \c h from this histogram. (See Subtract())
  CollieHistogram2d& operator-=(const CollieHistogram2d& h) {
    Subtract(h);
    return *this;
  }
  /// Scales this histogram by \c s. (See CollieHistogram::Scale())
  CollieHistogram2d& operator*=(double s) {
    Scale(s);
    return *this;
  }

  virtual const char* getType() {
    return "CollieHistogram2d";
  }
};

#endif
