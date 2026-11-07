#ifndef FileSetLoader_hh_included
#define FileSetLoader_hh_included

#include <Loader.hh>

/** This class manages a group of Loaders and produces just one
    SBD out of the entire set.

    Future versions should correctly manage systematic errors.
*/
class FileSetLoader : public Loader {
public:
  /// constructor
  FileSetLoader();
  /// destructor
  virtual ~FileSetLoader();
  /// add a loader to the list to manage
  void addLoader(Loader* aLoader);
  /// enable the use of log10(signal/background) rebinning which can greatly speed confidence level calculation
  void rebin(double min_l10_sob, double max_l10_sob, int bins);

  virtual SigBkgdDist* get(int var1);
  virtual SigBkgdDist* get(int var1, int var2);
  virtual SigBkgdDist* get(int var1, int var2, int var3);
  virtual int getMasspointList(int nPoints, int* v1, int* v2, int* v3);

  /// obtain the number of mass points
  virtual int getNMasspoints();

private:
  Loader** m_loaders;
  int n_loaders;
  bool f_rebin;
  double min_sob, max_sob;
  int bins_sob;
  SigBkgdDist* p_sbd;
  void* m_set;
};

#endif // FileSetLoader_hh_included
