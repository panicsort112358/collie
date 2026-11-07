#include "CollieIterator.hh"

CollieIterator::CollieIterator(CollieChannel* chan, int n, const int* v1, const int* v2,
                               const int* v3)
    : p_var1(v1),
      p_var2(v2),
      p_var3(v3) {
  p_Channel = chan;
  n_pos = 0;
  n_MassPoints = n;
  //  printf("%p %i %p \n",v1,v1[0],chan);
}
