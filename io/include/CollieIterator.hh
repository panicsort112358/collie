#ifndef CollieIterator_hh_included
#define CollieIterator_hh_included

#include <CollieMasspoint.hh>
#include <CollieChannel.hh>

/** Utility class for stepping through all the CollieMasspoint entries
    in a CollieChannel.
*/
class CollieIterator {
public:
  /// returns true if there are more mass points
  inline bool hasNext() const {
    return n_pos < n_MassPoints;
  }
  /// advances the iterator and loads the Masspoint
  inline const CollieMasspoint* next() {
    if (!hasNext())
      return NULL;
    n_pos++;
    //return NULL;
    return p_Channel->getMasspoint(p_var1[n_pos - 1], (p_var2 == NULL) ? -1 : p_var2[n_pos - 1],
                                   (p_var3 == NULL) ? -1 : p_var3[n_pos - 1]);
  }
  /// advances the iterator and does not load the Masspoint
  inline bool nextNoMasspoint() {
    if (!hasNext())
      return false;
    n_pos++;
    return true;
  }

  /// get the current value of first independent variable
  inline int getVar1() const {
    return (n_pos <= 0) ? -1 : p_var1[n_pos - 1];
  }
  /// get the current value of second independent variable
  inline int getVar2() const {
    return (n_pos <= 0 || p_var2 == NULL) ? -1 : p_var2[n_pos - 1];
  }
  /// get the current value of third independent variable
  inline int getVar3() const {
    return (n_pos <= 0 || p_var3 == NULL) ? -1 : p_var3[n_pos - 1];
  }

  /// constructor (should be called by CollieChannel)
  CollieIterator(CollieChannel* chan, int n, const int* v1, const int* v2, const int* v3);

private:
  int n_MassPoints;
  int n_pos;
  const int* p_var1;
  const int* p_var2;
  const int* p_var3;
  CollieChannel* p_Channel;
};

#endif // CollieIterator_hh_included
