#include <stdlib.h>
#include "FileSetLoader.hh"
#include <set>

struct triple_point {
  int v1, v2, v3;
  bool operator()(const triple_point& a, const triple_point& b) const {
    return (a.v1 < b.v1 || a.v2 < b.v2 || a.v3 < b.v3);
  }
};

inline static std::set<triple_point, triple_point>* toSet(void* v) {
  return (std::set<triple_point, triple_point>*)v;
}

FileSetLoader::FileSetLoader() {
  m_loaders = NULL;
  n_loaders = 0;
  f_rebin = false;
  p_sbd = NULL;
  m_set = new std::set<triple_point, triple_point>;
}

FileSetLoader::~FileSetLoader() {
  if (m_loaders != NULL)
    delete[] m_loaders;
  m_loaders = NULL;
  if (p_sbd != NULL)
    delete p_sbd;
  p_sbd = NULL;
  if (m_set != NULL) {
    std::set<triple_point, triple_point>* a = toSet(m_set);
    delete a;
    m_set = NULL;
  }
}

void FileSetLoader::addLoader(Loader* aLoader) {
  Loader** temp = new Loader*[n_loaders + 1];
  for (int i = 0; i < n_loaders; i++)
    temp[i] = m_loaders[i];
  if (m_loaders != NULL)
    delete[] m_loaders;
  m_loaders = temp;
  temp[n_loaders] = aLoader;
  n_loaders++;

  int len = aLoader->getNMasspoints();
  if (len <= 0)
    return;

  int v1[len], v2[len], v3[len];
  aLoader->getMasspointList(len, v1, v2, v3);
  for (int i = 0; i < len; i++) {
    triple_point a;
    a.v1 = v1[i];
    a.v2 = v2[i];
    a.v3 = v3[i];
    toSet(m_set)->insert(a);
  }
}

int FileSetLoader::getMasspointList(int nPoints, int* v1, int* v2, int* v3) {
  if (nPoints < (int)toSet(m_set)->size())
    return false;
  int j = 0;
  for (std::set<triple_point, triple_point>::iterator i = toSet(m_set)->begin();
       i != toSet(m_set)->end(); i++) {
    v1[j] = i->v1;
    if (v2 != NULL)
      v2[j] = i->v2;
    if (v3 != NULL)
      v3[j] = i->v3;
    j++;
  }
  return true;
}

int FileSetLoader::getNMasspoints() {
  return toSet(m_set)->size();
}

void FileSetLoader::rebin(double min_l10_sob, double max_l10_sob, int bins) {
  min_sob = min_l10_sob;
  max_sob = max_l10_sob;
  bins_sob = bins;
  f_rebin = true;
}

SigBkgdDist* FileSetLoader::get(int var1) {
  if (p_sbd != NULL && p_sbd->var1() == var1)
    return p_sbd;

  if (p_sbd != NULL)
    delete p_sbd;
  p_sbd = new SigBkgdDist(var1);
  if (f_rebin)
    p_sbd->rebin_sob(min_sob, max_sob, bins_sob);

  for (int i = 0; i < n_loaders; i++) {
    SigBkgdDist* sbd = m_loaders[i]->get(var1);
    if (sbd == NULL)
      continue;
    //    if (f_rebin) p_sbd->sob_append(*sbd);
    else
      p_sbd->append(*sbd);
    delete sbd;
  }

  return p_sbd;
}

SigBkgdDist* FileSetLoader::get(int var1, int var2) {
  if (p_sbd != NULL && p_sbd->var1() == var1 && p_sbd->var2() == var2)
    return p_sbd;

  if (p_sbd != NULL)
    delete p_sbd;
  p_sbd = new SigBkgdDist(var1, var2);
  if (f_rebin)
    p_sbd->rebin_sob(min_sob, max_sob, bins_sob);

  for (int i = 0; i < n_loaders; i++) {
    SigBkgdDist* sbd = m_loaders[i]->get(var1, var2);
    if (sbd == NULL)
      continue;
    //    if (f_rebin) p_sbd->sob_append(*sbd);
    else
      p_sbd->append(*sbd);
    delete sbd;
  }

  return p_sbd;
}

SigBkgdDist* FileSetLoader::get(int var1, int var2, int var3) {
  if (p_sbd != NULL && p_sbd->var1() == var1 && p_sbd->var2() == var2 && p_sbd->var3() == var3)
    return p_sbd;

  if (p_sbd != NULL)
    delete p_sbd;
  p_sbd = new SigBkgdDist(var1, var2, var3);
  if (f_rebin)
    p_sbd->rebin_sob(min_sob, max_sob, bins_sob);

  for (int i = 0; i < n_loaders; i++) {
    SigBkgdDist* sbd = m_loaders[i]->get(var1, var2, var3);
    if (sbd == NULL)
      continue;
    //    if (f_rebin) p_sbd->sob_append(*sbd);
    else
      p_sbd->append(*sbd);
    delete sbd;
  }

  return p_sbd;
}
