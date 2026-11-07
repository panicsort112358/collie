#ifndef CDF_DZero_IOpoint_hh_included
#define CDF_DZero_IOpoint_hh_included

#include "TObjArray.h"
#include "TH1.h"
#include <TString.h>
#include "TRef.h"
#include "TNamed.h"
#include "string"
#include "CDF_DZero_Distribution.hh"
#include <cmath>

class CDF_DZero_IOfile;

class CDF_DZero_IOpoint : public TNamed {

public:
  /// constructor for streaming
  CDF_DZero_IOpoint(); // for streaming...

  /// destructor
  virtual ~CDF_DZero_IOpoint();

  void Book(const char* name, int nx, double minX, double maxX);
  void Book(const char* name, int nx, double minX, double maxX, int ny, double minY, double maxY);

  /// get the number of the bins in the X direction in each distribution
  inline int getNXbins() const {
    return fNXbins;
  }
  /// get the number of the bins in the Y direction in each distribution
  inline int getNYbins() const {
    return fNYbins;
  }
  /// get the maximum edge value in the X direction in each distribution
  inline double getMaxX() const {
    return fMaxX;
  }
  /// get the minimum edge value in the X direction in each distribution
  inline double getMinX() const {
    return fMinX;
  }
  /// get the maximum edge value in the Y direction in each distribution
  inline double getMaxY() const {
    return fMaxY;
  }
  /// get the minimum edge value in the Y direction in each distribution
  inline double getMinY() const {
    return fMinY;
  }

  /// get the value of the first independent variable(coordinate) for this mass point.  Usually this is a mass.
  inline int getVar1() const {
    return fIndepVar1;
  }
  /// get the value of the second independent variable(coordinate) for this mass point, if one is used
  inline int getVar2() const {
    return fIndepVar2;
  }
  /// get the value of the third independent variable(coordinate) for this mass point, if one is used
  inline int getVar3() const {
    return fIndepVar3;
  }

  /// add a signal distribution
  void addSignalDist(CDF_DZero_Distribution* dist);
  /// access a signal distribution by index number
  CDF_DZero_Distribution* getSignalDist(int i = 0);
  /// get the number of signal distributions
  inline int getNSignalDists() const {
    return fSignals.GetEntriesFast();
  }

  /// add a signal distribution
  void addBkgdDist(CDF_DZero_Distribution* dist);
  /// access a background distribution by index number
  CDF_DZero_Distribution* getBkgdDist(int i = 0);
  /// get the number of background distributions
  inline int getNBkgdDists() const {
    return fBackgrounds.GetEntries();
  }

  /// enter a data distribution
  void addDataDist(TH1* data);
  void addDataDist2D(TH2* data);

  /// access the data
  TH1* getDataDist() const {
    return (TH1*)fDataDistribution;
  }
  TH2* getDataDist2D() const {
    return (TH2*)fDataDistribution;
  }

  void checkPoint();

  void print(map<string, int>& nSyst, map<string, double>& posSyst, map<string, double>& negSyst);

protected:
  /// constructor for creation
  static std::string buildName(const char* channelname, int var1, int var2, int var3, int nvar);

  CDF_DZero_IOpoint(int var1, int var2, int var3);

  friend class CDF_DZero_IOfile;

private:
  bool checkDist(CDF_DZero_Distribution* dist);
  bool checkHist(TH1* hist);
  bool checkHist(TH2* hist);

  // distributions for the signals
  TObjArray fSignals;
  // distributions for the backgrounds
  TObjArray fBackgrounds;

  TH1* fDataDistribution; // distribution for the data (if stored as a distribution)

  Int_t fIndepVar1; // first independent variable for this mass point (usually a mass)
  Int_t fIndepVar2; // second independent variable for this mass point (if used)
  Int_t fIndepVar3; // third independent variable for this mass point (if used)

  Int_t fNXbins; // number of bins in X (also stored in each distribution)
  Int_t fNYbins; // number of bins in Y (also stored in each distribution)

  Double_t fMinX; // minimum value of distribution in X
  Double_t fMaxX; // maximum value of distribution in X
  Double_t fMinY; // minimum value of distribution in Y (if relevant)
  Double_t fMaxY; // maximum value of distribution in Y (if relevant)

  ClassDef(CDF_DZero_IOpoint, 1)
};

#endif // CDF_DZero_IOpoint_hh_included
