#ifndef Loader_HH_INCLUDED
#define Loader_HH_INCLUDED

#include "CollieMasspoint.hh"
#include "SigBkgdDist.hh"

/** Abstract base class defining the required behavior of a Loader for a file
type.
*/
class Loader {
public:
  /// constructor
  Loader();
  /// destructor
  virtual ~Loader();

  /// open a file
  virtual bool open(const char* filename, const char* options = 0);

  virtual void setbkgdscales(int ns, double* sc) = 0;

  virtual SigBkgdDist* MassPoint2SBD(CollieMasspoint* mp) = 0;

  /// get the nearest distribution to this variable (mass) value
  virtual SigBkgdDist* get(int var1) = 0;

  /// get the nearest distribution to these variable values (mass, tan beta)
  virtual SigBkgdDist* get(int var1, int var2) = 0;

  /// get the nearest distribution to these variable values
  virtual SigBkgdDist* get(int var1, int var2, int var3) = 0;

  /** \brief get a list of the mass points contained in this file into user-supplied arrays.
      \param nPoints number of entries in the arrays
      \param v1 array for the first independent variable
      \param v2 array for the second independent variable or NULL if not needed
      \param v3 array for the third independent variable of NULL if not needed
      \return False if unsupported, or nPoints is too small.
  */
  virtual int getMasspointList(int nPoints, int* v1, int* v2, int* v3) {
    return false;
  }

  /// obtain the number of mass points
  virtual int getNMasspoints() {
    return -1;
  }

  /// get the name of the file associated with this Loader
  inline const char* getFilename() const {
    return m_filename;
  }
  /// get the options string used by this loader
  inline const char* getOptions() const {
    return m_options;
  }
  /// test the happiness of the loader (implementation dependent)
  virtual bool ok() {
    return true;
  }

  int nscales;
  double scales[100];

private:
  char* m_filename;
  char* m_options;
};

#endif // Loader_HH_INCLUDED
