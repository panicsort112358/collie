#ifndef CollieEventList_hh_included
#define CollieEventList_hh_included

#include "CollieDistribution.hh"
#include "TNamed.h"
#include "TArrayI.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"

class CollieMasspoint;

/** \brief Final variable distribution on an event-by-event basis.

    This class represents a list of events with X (and optionally Y) values
    which map to the values used in the CollieDistribution objects in the
    same CollieMasspoint.

    Additional doubleing-point and integer values may also be stored for each
    selected data event in the list.

*/
class CollieEventList : public TNamed {
public:
  /// constuctor for streaming
  CollieEventList(); // for streaming
  /// copy constructor
  //  CollieEventList(const CollieEventList& d);
  /// copy assignment operator
  //  CollieEventList& operator=(const CollieEventList& d);

  /// creates a new Root Histogram (1D) and fills it with the contents of this distribution
  TH1D* draw(const char* title, int nbins, double xlow, double xhi);
  /// creates a new Root Histogram (2D) and fills it with the contents of this distribution
  TH2D* draw2D(const char* title, int nbinsx, double xlow, double xhi, int nbinsy, double ylow,
               double yhi);

  /// add an event to the list (single dimensional version)
  void addEvent(double x, double* dvars = NULL, int* ivars = NULL);
  /// add an event to the list (two dimensional version)
  void addEvent(double x, double y, double* dvars = NULL, int* ivars = NULL);
  /// set the name of a doubleing point variable
  void setDVarName(int nvar, const char* name);
  /// set the name of an integer variable
  void setIVarName(int nvar, const char* name);
  /// get the name of a doubleing point variable
  const char* getDVarName(int nvar) const {
    return (nvar < 0 || nvar >= fNdoubles) ? NULL
                                           : (((TObjString*)fValueNames[nvar])->GetString().Data());
  }
  /// get the name of an integer variable
  const char* getIVarName(int nvar) const {
    return (nvar < 0 || nvar >= fNints)
               ? NULL
               : (((TObjString*)fValueNames[fNdoubles + nvar])->GetString().Data());
  }

  /// get the number of events in the list
  inline int getNEvents() const {
    return fNevents;
  }
  /// get the X value for the given event
  inline double getX(int nevent) const {
    return getDVar(nevent, 0);
  }
  /// get the Y value for the given event
  inline double getY(int nevent) const {
    return getDVar(nevent, 1);
  }

  /// get a doubleing point variable by number (X is 0, Y [if used] is 1)
  inline double getDVar(int nevent, int nvar) const {
    return (nvar < 0 || nvar >= fNdoubles || nevent < 0 || nevent >= fNevents)
               ? 0
               : fDoubles[nvar * fNevents + nevent];
  }
  /// get a doubleing point variable by name
  inline double getDVar(int nevent, const char* name) const {
    return getDVar(nevent, lookupDVarName(name));
  }

  /// get an integer variable by number
  inline int getIVar(int nevent, int nvar) const {
    return (nvar < 0 || nvar >= fNints || nevent < 0 || nevent >= fNevents)
               ? 0
               : fIntegers[nvar * fNevents + nevent];
  }
  /// get an integer variable by name
  inline int getIVar(int nevent, const char* name) const {
    return getIVar(nevent, lookupIVarName(name));
  }

  /// convert into a distribution
  void fillDistribution(CollieDistribution* cd, int overflow = 1) const;

  /// access the mass point object which contains this event list
  const CollieMasspoint* getMasspoint() const {
    return p_MassPoint;
  }

  inline int nIVars(void) const {
    return fNints;
  }
  inline int nDVars(void) const {
    return fNdoubles;
  }
  inline int nDimensions(void) const {
    return fDimensionality;
  }

  /// destructor.
  ~CollieEventList();

private:
  CollieEventList(const char* name, int dims, int nF, int nI);
  friend class CollieMasspoint;
  bool fMutable;                      //! whether this distribution can be changed
  const CollieMasspoint* p_MassPoint; //! pointer to the mass point (set by CollieMasspoint)
  TString fName;                      // name of the distribution
  // Name is based on CollieMasspoint name, distribution type, and more
  Int_t fDimensionality; // number of doubleing point values which are final variable values
  Int_t fNdoubles;       // number of doubleing point values per event
  Int_t fNints;          // number of integer values per event
  Int_t fNevents;        // number of events

  TArrayD fDoubles;      // doubleing point values
  TArrayI fIntegers;     // integer values
  TObjArray fValueNames; // names of the various values

  int lookupIVarName(const char* name) const;
  int lookupDVarName(const char* name) const;

  void BookDataIntl();

  ClassDef(CollieEventList, 1)
};

#endif // CollieEventList_hh_included
