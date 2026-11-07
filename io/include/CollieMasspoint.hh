#ifndef CollieMasspoint_hh_included
#define CollieMasspoint_hh_included

#include <TObjArray.h>
#include <TRef.h>
#include <TNamed.h>
#include <CollieDistribution.hh>
#include <CollieEventList.hh>
#include <string>

class CollieChannel;

/** \brief Distribution collection for a specific mass/variable hypothesis.

    A CollieMasspoint represents the distributions for a given set of
    inpendent variables.  The structure of all mass points in a
    channel is specified by the associated CollieChannel.  The channel
    specifies the number (and name) of the signal and background distributions.
    The mass point provides access each of these distributions through
    CollieDistribution objects.
*/
class CollieMasspoint : public TNamed {
public:
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

  /// access a signal distribution by index number (read-only)
  const CollieDistribution* getSignalDist(int i = 0) const;
  /// access a signal distribution by name (read-only)
  inline const CollieDistribution* getSignalDist(const char* name) const {
    return getSignalDist(lookupSignal(name));
  }
  /// get the number of signal distributions
  inline int getNSignalDists() const {
    return fSignals.GetEntriesFast();
  }

  /// access a background distribution by index number (read-only)
  const CollieDistribution* getBkgdDist(int i = 0) const;
  /// access a background distribution by name (read-only)
  inline const CollieDistribution* getBkgdDist(const char* name) const {
    return getBkgdDist(lookupBkgd(name));
  }
  /// get the number of background distributions
  inline int getNBkgdDists() const {
    return fBackgrounds.GetEntries();
  }

  /// access the data in the form of a distribution (read-only).  If the data is in the form of an event list, this method will convert it into an appropriate distribution
  const CollieDistribution* getDataDist() const;
  /// if possible, access the event list (read-only)
  const CollieEventList* getDataEventList() const;

  /// if possible, access a signal distribution by index number in read/write mode
  CollieDistribution* getSignalDistMutable(int i = 0);
  /// if possible, access a background distribution by index number in read/write mode
  CollieDistribution* getBkgdDistMutable(int i = 0);
  /// if possible, access the data distribution in read/write mode
  CollieDistribution* getDataDistMutable();
  /// if possible, access the event list in read/write mode
  CollieEventList* getDataEventListMutable();

  void scaleSignal(double sf);
  void scaleBackground(int idx, double sf);
  void scaleData(double sf);

  /// store the modified mass point into the current ROOT file
  void Store();
  void logPoint();

  /** if the mass point is new, specify the structure of the distributions
      including the number of bins and the limits of the distribution.  This
      version of the Book command handles only 1d distributions.

      \warning A Book method must be called after a mass point is created before accessing any of the distributions.
  */
  void Book(int nx, double minX, double maxX);
  /** if the mass point is new, specify the structure of the distributions
      including the number of bins and the limits of the distribution.  This
      version of the Book command allows 2d distributions.

      \warning A Book method must be called after a mass point is created before accessing any of the distributions.
 */
  void Book(int nx, double minX, double maxX, int ny, double minY, double maxY);
  /// if the mass point is new, create a data distribution of the same format
  void BookData();
  /** if the mass point is new, specify that the data be stored in the form
      of an event list with the given format.

      \param ndoubleVars Number of additional doubleing-point variables to provide in the event list.
      \param nintVats Number of additional integer variables to provide in the vent list.
      \warning A Book method must be called before calling this method (to establish the dimensionality of the distributions)
  */
  void ConfigEventList(int ndoubleVars, int nintVars);

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

  /// access the CollieChannel associated with this mass point
  inline const CollieChannel* getChannel() const {
    return pChannel;
  }

  inline void setMutable(bool mut) {
    fMutable = mut;
  }

  // print information about this masspoint
  void print(map<string, int>& count, map<string, double>& posSyst, map<string, double>& negSyst);

  /// constructor for streaming
  CollieMasspoint(); // for streaming...
  /// destructor
  virtual ~CollieMasspoint();

protected:
  /// assemble the name used for a mass point from a channel name and variable values
  static std::string buildName(const char* channelname, int var1, int var2, int var3, int nvar);
  /// constructor used by the CollieChannel when creating mass points
  CollieMasspoint(CollieChannel* theChannel, int var1, int var2, int var3);

  friend class CollieChannel;

private:
  int lookupBkgd(const char* name) const;
  int lookupSignal(const char* name) const;
  CollieChannel* pChannel; //! the Channel...

  bool fMutable; //! am I changable?

  Int_t fIndepVar1; // first independent variable for this mass point (usually a mass)
  Int_t fIndepVar2; // second independent variable for this mass point (if used)
  Int_t fIndepVar3; // third independent variable for this mass point (if used)

  Int_t fNXbins; // number of bins in X (also stored in each distribution)
  Int_t fNYbins; // number of bins in Y (also stored in each distribution)

  Int_t fSigSyst; // number of signal systematics
  Int_t fBkgSyst; // number of background systematics

  Double_t fMinX; // minimum value of distribution in X
  Double_t fMaxX; // maximum value of distribution in X
  Double_t fMinY; // minimum value of distribution in Y (if relevant)
  Double_t fMaxY; // maximum value of distribution in Y (if relevant)

  // distributions for the signals
  TObjArray fSignals;
  // distributions for the backgrounds
  TObjArray fBackgrounds;

  CollieDistribution* fDataDistribution; // distribution for the data (if stored as a distribution)
  CollieEventList* fDataEventList;       // event list for the data

  CollieDistribution* BookDataIntl() const;

  ClassDef(CollieMasspoint, 2)
};

#endif // CollieMasspoint_hh_included
