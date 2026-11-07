#ifndef TrialPoint_HH_INCLUDED
#define TrialPoint_HH_INCLUDED

#include "findSoleDifference.hh"
#include <vector>

class TrialPoint {
  typedef std::vector<double> seq;

public:
  TrialPoint()
      : last_point() {}

  bool is_clear() const {
    return last_point.empty();
  }

  void clear() {
    last_point.clear();
  }

  seq::const_iterator diff(seq const& this_point) {
    seq::const_iterator result = this_point.end();

    if (this_point.size() == last_point.size())
      result = findSoleDifference<seq::const_iterator>(this_point.begin(), this_point.end(),
                                                       last_point.begin());

    if (result == this_point.end())
      last_point = this_point;

    return result;
  }

private:
  seq last_point;

}; // TrialPoint

#endif // TrialPoint_HH_INCLUDED
