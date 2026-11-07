#ifndef CollieLoader_HH_INCLUDED
#define CollieLoader_HH_INCLUDED

#include "Loader.hh"

class CollieChannel;

/** Loader for CollieIO files. */
class CollieLoader : public Loader {
public:
  /// constructor
  CollieLoader();

  /** Required option:
      name='(name of the channel from the file)'
  */
  virtual bool open(const char* filename, const char* options = 0);
  virtual bool close();

  virtual void setbkgdscales(int nscales, double* scales);
  virtual SigBkgdDist* MassPoint2SBD(CollieMasspoint* mp);
  virtual SigBkgdDist* get(int var1);
  virtual SigBkgdDist* get(int var1, int var2);
  virtual SigBkgdDist* get(int var1, int var2, int var3);

  virtual int getMasspointList(int nPoints, int* v1, int* v2, int* v3);

  virtual int getNMasspoints();

  virtual bool ok() {
    return m_channel != NULL;
  }

  /// provides access to contained CollieChannel object
  inline CollieChannel* accessChannel() {
    return m_channel;
  }

private:
  CollieChannel* m_channel;
  TFile* m_file;
};

#endif // CollieLoader_HH_INCLUDED
