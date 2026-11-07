#ifndef CollieChannel_hh_included
#define CollieChannel_hh_included

#include <TNamed.h>
#include <TString.h>
#include <TFile.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TDatime.h>
#include <TArrayI.h>
#include <string>
#include <vector>
#include "CollieMasspoint.hh"

class CollieIterator;

/** \brief Header class for final-variable distributions for a specific search channel

    Structurally, the data for a channel are stored together in a Directory of
    a ROOT/IO file.  The CollieChannel and all associated CollieMasspoints are
    stored with names which depend on the name given to the channel.

 */
class CollieChannel : public TNamed {
public:
  /// constuctor for streaming
  CollieChannel();

  /// destructor
  ~CollieChannel();

  /** Create a new channel.  The channel pointer created here is
      read/write(mutable).  \param channelName Name to associate with
      this channel.  This name should not contain any spaces, tabs, or
      any of the following characters: '(',')','"',''',':'.  Any name
      containing these letters will be rejected.  \param rootFile
      Pointer to a TFile to create the channel in.  \param options
      String containing options.  Currently none are defined.  \return
      A new channel if the creation was successful, or NULL if it was
      not.  \warning After setting up the structure of the channel,
      you must call Store() to save the CollieChannel to the file!
  */
  static CollieChannel* createChannel(const char* channelName, TFile* rootFile,
                                      const char* options = NULL);

  /** Load a channel from a file.  The channel pointer obtained from this
      call is read-only and the contents cannot be changed.
      \param channelName Name given to the channel when it was created.
      \param rootFile Pointer to a TFile containing the channel data.
      \return Pointer to the channel if successful, NULL if the channel was not found.
  */
  static CollieChannel* loadChannel(const char* channelName, TFile* rootFile);

  CollieChannel* CopyChannel(TFile* rootFile);

  /// If this channel is read/write (mutable), set the comments string associated with the channel.
  inline void setComments(const char* comments) {
    if (fMutable)
      fComments = comments;
  }
  /** \brief If this channel is read/write (mutable), set the luminosity value associated with the channel.
      \warning The units of the luminosity should be consistent with the units
      of the cross-sections in the CollieDistributions.
  */
  inline void setLuminosity(double lumi) {
    if (fMutable)
      fLuminosity = lumi;
  }

  /** \brief If this channel is read/write (mutable), set the number of independent
      variables which will define the set of Masspoints.

      This is <b>not</b>
      the dimensionality of the final variable distributions, but rather the
      number of coordinates which define the scan line, plane, or volume.
      For example, in a Standard Model Higgs search, there is one
      variable: mH; while in an MSSM Higgs search, there are
      two: mH and tan(beta).
      The maximum number of independent variables is three.  Independent
      variables are stored as integers.
  */
  void setNIndepVariables(int nVariables);
  /** \brief If this channel is read/write (mutable), set the number of different
      types of signal distribution which can
      be added together as the "signal" for this search.
  */
  void setNSignals(int nSignals);
  /** \brief If this channel is read/write (mutable), set the number of different
      types of background distribution which should
      be added together as the "background" for this search.
  */
  void setNBackgrounds(int nBackgrounds);

  /** \brief If this channel is read/write (mutable), set the number of different
      models (sets of cross-sections) which this channel may use.
  */
  void setNmodels(int nModels);

  /// store the channel into the file (if mutable)
  void Store();

  /** \brief Create a new mass point indexed by the given independent variables.
      \warning The CollieMasspoint object must be Booked before accessing any distribution.
      \return The new mass point object if successful, NULL if not.
  */
  CollieMasspoint* createMasspoint(int var1, int var2 = 0, int var3 = 0);

  //  CollieMasspoint* createMasspoint(CollieMasspoint*);

  /** \brief Look up the mass point indexed by the given independent variables.
      \return The mass point object if successful, NULL if not.
  */
  CollieMasspoint* getMasspoint(int var1, int var2 = 0, int var3 = 0);
  /** \brief Get an iterator over the mass points in the channel */
  CollieIterator* getIterator();
  /// get a count of the mass points in this channel
  inline int getNMasspoints() const {
    return fMasspoints_var1.GetSize();
  }

  /// if this channel is read/write (mutable), set the name associated with the indexed independent variable
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

  /// if this channel is read/write (mutable), set the name associated with the indexed signal distribution
  inline void setSignalName(int nSignal, const char* name) {
    if (fMutable && nSignal >= 0 && nSignal < fNsignals)
      ((TObjString*)fSignalNames[nSignal])->SetString((char*)name);
  }
  /// if this channel is read/write (mutable), set the name associated with the indexed background distribution
  inline void setBackgroundName(int nBkgd, const char* name) {
    if (fMutable && nBkgd >= 0 && nBkgd < fNbackgrounds)
      ((TObjString*)fBkgdNames[nBkgd])->SetString((char*)name);
  }

  /// if this channel is read/write (mutable), set the name associated with the indexed (cross-section) model
  inline void setModelName(int nModel, const char* name) {
    if (fMutable && nModel >= 0 && nModel < fNmodels)
      ((TObjString*)fModelNames[nModel])->SetString((char*)name);
  }

  /** \brief get the name associated with the indexed variable.
      \see setNIndepVariables()
  */
  inline const char* getIndepVariableName(int i) const {
    return fIndepVarName[i % 3].Data();
  }
  /** \brief get the number of independent variables (max 3)
      \see setNIndepVariables()
  */
  inline int getNIndepVariables() const {
    return fNindepVars;
  }

  /// get the name of the channel
  inline const char* getChannelName() const {
    return fChannelName.Data();
  }
  /// get the name of the channel
  inline TString getChannelNameT() const {
    return fChannelName;
  }
  /// get the comments associated with this channel
  inline const char* getComments() const {
    return fComments.Data();
  }
  /** \brief get the luminosity totaled in this channel
      \warning The units of the luminosity should be consistent with the units
      of the cross-sections in the CollieDistributions.
  */
  inline double getLuminosity() const {
    return fLuminosity;
  }
  /// get the number of signal distributions
  inline int getNSignals() const {
    return fNsignals;
  }
  /// get the name associated with a given indexed signal distribution
  inline const char* getSignalName(int i = 0) const {
    return (i < 0 || i >= fNsignals) ? NULL : ((TObjString*)fSignalNames[i])->GetString().Data();
  }
  /// get the number of background distributions
  inline int getNBackgrounds() const {
    return fNbackgrounds;
  }
  /// get the name associated with a given indexed background distribution
  inline const char* getBackgroundName(int i = 0) const {
    return (i < 0 || i >= fNbackgrounds) ? NULL : ((TObjString*)fBkgdNames[i])->GetString().Data();
  }

  void scaleSignal(double sf);
  void scaleBackground(int idx, double sf);
  void scaleData(double sf);

  /// get the number of (cross-section) models
  inline int getNModels() const {
    return fNmodels;
  }
  /// get the name associated with a (cross-section) model
  inline const char* getModelName(int i = 0) const {
    return (i < 0 || i >= fNmodels) ? NULL : ((TObjString*)fModelNames[i])->GetString().Data();
  }

  /// get the host name recorded during the initial creation of the channel object
  inline const char* getCreationComputer() const {
    return fCreationComputer.Data();
  }
  /// get the user login name recorded during the initial creation of the channel object
  inline const char* getCreationUser() const {
    return fCreationUser.Data();
  }
  /// get the Collie CVS tag version used to make this file
  //  inline const char* getCollieVersion() const { return fCollieVersion.Data(); }
  /// get the time and date recorded during the initial creation of the channel object
  inline const TDatime& getCreationTime() const {
    return fCreationTime;
  }

  /// utility routine
  inline void doCD() {
    fDirectory->cd();
  }

  // print all channel information
  void print();

protected:
  /// store a mass point and record it in the lists...
  void storePoints();
  void logPoint(CollieMasspoint* cmp);

  friend class CollieMasspoint;

private:
  CollieChannel(const char* name);
  void storePoint(CollieMasspoint* cmp);

  bool fMutable;          //! whether this CollieChannel object can be changed
  TDirectory* fDirectory; //! pointer to the TDirectory

  TString fChannelName;   // the name of the channel
  TString fComments;      // any comments about this channel
  TString fCollieVersion; // the collie version used to make the file
  Double_t fLuminosity;   // the luminosity of the channel

  Int_t fNindepVars;        // number of independent variables ( up to three )
  TString fIndepVarName[3]; // names of the (up to three) independent variables

  Int_t fNsignals;        // the number of signal distributions
  TObjArray fSignalNames; // names of the signals
  Int_t fNbackgrounds;    // the number of background distributions
  TObjArray fBkgdNames;   // names of the backgrounds

  Int_t fNmodels;        // the number of models for production, etc
  TObjArray fModelNames; // names of the models

  TString fCreationComputer; // computer on which this channel was created
  TString fCreationUser;     // user who created this channel (may be bogus)
  TDatime fCreationTime;     // when this channel was created (not final modification!)
  //  TString fCollieVersion; // CVS tag version of collie package used to create the channel

  TArrayI fMasspoints_var1; // array/index of mass points (variable 1)
  TArrayI fMasspoints_var2; // array/index of mass points (variable 2 if used)
  TArrayI fMasspoints_var3; // array/index of mass points (variable 3 if used)

  void associateFile(TFile* tf);

  ClassDef(CollieChannel, 3)
};

#endif // CollieChannel_hh_included
