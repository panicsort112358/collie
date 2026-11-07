#ifndef CDF_DZero_Distribution_hh_included
#define CDF_DZero_Distribution_hh_included

#include "TObjArray.h"
#include "TObjString.h"
#include "TH1.h"
#include "TH2.h"
#include "TRef.h"
#include "TNamed.h"
#include <map>

using namespace std;

class CDF_DZero_Distribution : public TNamed {

public:
  /// constructor for streaming
  CDF_DZero_Distribution(); // for streaming...
  /// destructor
  virtual ~CDF_DZero_Distribution();

  /// add a distribution
  void addDistribution(TH1* dist);
  void addDistribution2D(TH2* dist);

  /// get the number of the bins in the X direction in each distribution
  int getNXbins() const;
  /// get the number of the bins in the Y direction in each distribution
  int getNYbins() const;
  /// get the maximum edge value in the X direction in each distribution
  double getMaxX() const;
  /// get the minimum edge value in the X direction in each distribution
  double getMinX() const;
  /// get the maximum edge value in the Y direction in each distribution
  double getMaxY() const;
  /// get the minimum edge value in the Y direction in each distribution
  double getMinY() const;

  // Is this a 2D distribution?
  inline bool get2Dflag() const {
    return f2Ddistribution;
  }

  /// access a distribution
  TH1* getDistribution() const {
    return (TH1D*)fDistribution;
  }
  TH2* getDistribution2D() const {
    return (TH2D*)fDistribution;
  }

  /// add a systematic distribution
  void addSystematic(const char* name, TH1* pos, TH1* neg);
  void addSystematic2D(const char* name, TH2* pos, TH2* neg);

  /// access a systematic distribution by index number
  TH1* getPositiveSystematic(int i);
  TH2* getPositiveSystematic2D(int i);

  /// access a systematic distribution by index number
  TH1* getNegativeSystematic(int i);
  TH2* getNegativeSystematic2D(int i);
  /// access a systematic name by index number
  const char* getSystematicName(int i);

  /// get the number of systematics
  inline int getNSystematics() const {
    return fSystematicsNames.GetEntriesFast();
  }

  double inBin(int i, int j = 0) const;

  void setDistName(const char* name);
  const char* getDistName();

  void print(map<string, int>& nSyst, map<string, double>& posSyst, map<string, double>& negSyst);

private:
  TH1* fDistribution; // distribution
  TObjArray fSystematicsPos;
  TObjArray fSystematicsNeg;
  TObjArray fSystematicsNames;
  bool f2Ddistribution;

  Int_t fNXbins; // number of bins in X (also stored in each distribution)
  Int_t fNYbins; // number of bins in Y (also stored in each distribution)

  Double_t fMinX; // minimum value of distribution in X
  Double_t fMaxX; // maximum value of distribution in X
  Double_t fMinY; // minimum value of distribution in Y (if relevant)
  Double_t fMaxY; // maximum value of distribution in Y (if relevant)

  ClassDef(CDF_DZero_Distribution, 1)
};

#endif // CDF_DZero_Distribution_hh_included
