#ifndef CDF_DZero_IOfile_hh_included
#define CDF_DZero_IOfile_hh_included

#include "TObjArray.h"
#include "TNamed.h"
#include "TFile.h"
#include <TDatime.h>
#include <TArrayI.h>
#include "CDF_DZero_IOpoint.hh"

class CDF_DZero_IOfile : public TNamed {

public:
  /// constructor for streaming
  CDF_DZero_IOfile(); // for streaming...

  static CDF_DZero_IOfile* createChannel(const char* channelName, TFile* rootFile,
                                         const char* options = NULL);

  static CDF_DZero_IOfile* loadChannel(const char* channelName, TFile* rootFile);

  void Store();

  void logPoint(CDF_DZero_IOpoint* point);

  CDF_DZero_IOpoint* createPoint(int var1, int var2 = 0, int var3 = 0);
  CDF_DZero_IOpoint* getPoint(int var1, int var2 = 0, int var3 = 0);
  inline int getNPoints() const {
    return fPoints_var1.GetSize();
  }

  void setNIndepVariables(int nVariables);
  inline int getNIndepVariables() const {
    return fNindepVars;
  }

  void getIndependentVariables(int idx, int& var1, int& var2, int& var3);

  inline const char* getIndepVariableName(int i) const {
    return fIndepVarName[i % 3].Data();
  }
  inline void setIndepVarName(int nvar, const char* name) {
    if (!fMutable)
      return;
    if (nvar == 0)
      fIndepVarName[0] = name;
    else if (nvar == 1)
      fIndepVarName[1] = name;
    else if (nvar == 2)
      fIndepVarName[2] = name;
  }

  /// get the name of the channel
  inline const char* getChannelName() const {
    return fChannelName.Data();
  }
  /// get the name of the channel
  inline TString getChannelNameT() const {
    return fChannelName;
  }

  /// get the host name recorded during the initial creation of the channel object
  inline const char* getCreationComputer() const {
    return fCreationComputer.Data();
  }
  /// get the user login name recorded during the initial creation of the channel object
  inline const char* getCreationUser() const {
    return fCreationUser.Data();
  }
  /// get the time and date recorded during the initial creation of the channel object
  inline const TDatime& getCreationTime() const {
    return fCreationTime;
  }

  // dump contents of file to screen
  void print();

  /// utility routine
  inline void doCD() {
    fDirectory->cd();
  }

private:
  CDF_DZero_IOfile(const char* name);
  void associateFile(TFile* tf);
  void storePoint(CDF_DZero_IOpoint* iop);
  void storePoints();

  bool fMutable;             //! whether this object can be changed
  TDirectory* fDirectory;    //! pointer to the TDirectory
  TString fChannelName;      // the name of the channel
  TString fCreationComputer; // computer on which this channel was created
  TString fCreationUser;     // user who created this channel (may be bogus)
  TDatime fCreationTime;     // when this channel was created (not final modification!)

  Int_t fNindepVars;        // number of independent variables ( up to three )
  TString fIndepVarName[3]; // names of the (up to three) independent variables

  TArrayI fPoints_var1; // array/index of mass points (variable 1)
  TArrayI fPoints_var2; // array/index of mass points (variable 2 if used)
  TArrayI fPoints_var3; // array/index of mass points (variable 3 if used)

  ClassDef(CDF_DZero_IOfile, 1)
};

#endif // CDF_DZero_IOfile_hh_included
