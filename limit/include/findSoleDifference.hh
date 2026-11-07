#ifndef findSoleDifference_HH_INCLUDED
#define findSoleDifference_HH_INCLUDED

#include <cmath>

#ifdef ALTERNATIVE1
template <class T> bool areDifferent(T a, T b) {
  static const double epsilon = 1.0e-13;
  T absa = (a > 0 ? a : -a);
  T absb = (b > 0 ? b : -b);
  T minab = (absa < absb ? absa : absb);
  if (minab == 0)
    return (a != b);
  T absdiff = (a > b ? a - b : b - a);
  return (absdiff > epsilon * minab);
}

template <class It> It findSoleDifference(It b1, It e1, It b2) {
  It result = b1;
  unsigned ndiffs = 0u;

  for (; b1 != e1; ++b1, ++b2)
    if (areDifferent(*b1, *b2)) {
      if (++ndiffs > 1u)
        return e1;
      result = b1;
    }
  return result;
}
#endif

template <class T> bool areDifferent(T a, T b) {
  return std::abs(a - b) > 1.0e-13 * std::abs(a);
}

template <class It> It findSoleDifference(It b1, It e1, It b2) {
  It result = b1;
  unsigned ndiffs = 0u;

  for (; b1 != e1; ++b1, ++b2)
    if (std::abs(*b1 - *b2) > 1.0e-13 * std::abs(*b1)) {
      if (++ndiffs > 1u)
        return e1;
      result = b1;
    }
  return result;
}

#endif
