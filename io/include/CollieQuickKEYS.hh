#ifndef CollieQuickKEYS_h_included
#define CollieQuickKEYS_h_included

#include "CollieHistogram.hh"
#include <math.h>

/** \brief Smoothing tool based on work of Kyle Cranmer
    \ingroup util
    This class is a smoothing tool based on the work of Kyle Cranmer.  It implements
    the same basic algorithm as the KEYS package, but uses C and C++ instead
    of FORTRAN, and has no dependency on PAW.  The class is designed to be
    mostly independant of CollieAnalLib, and should be removable as a single
    package.

  Features:
    1) reflection and rescaling for boundary purposes
    2) creation of histogram from smoothed distribution

This is a version of the same algorithm as in the original CollieKeys (now
obsolete).  However, it runs much quicker in most cases, for a small cost
in accuracy.
*/
class CollieQuickKEYS : public CollieHistogram {
public:
  // overrides
  virtual int Book(int id, const char* name, int bins, double left, double right,
                   CollieHistogramContainer* container = NULL);
  virtual int Book(const char* name, int bins, double left, double right,
                   CollieHistogramContainer* container = NULL);

  virtual void Fill(double x, double w = 1.0);
  virtual void Fill(int i, double w = 1.0);
  virtual void Clear();

  /// smooth with KEYS...
  void Smooth(double m_alpha = 1.0, bool refl_f0 = false);
  /// copy into array of doubles, optionally rebinning
  void Copy(double*, int bins_per_bin = 1);
  /// constructor
  CollieQuickKEYS();
  /// destructor
  virtual ~CollieQuickKEYS();

private:
  double* m_centers;
  int n_fills;
};

#endif // CollieQuickKEYS_h_included
