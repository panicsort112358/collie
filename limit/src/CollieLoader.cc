#include "CollieIterator.hh"
#include "CollieChannel.hh"
#include "CollieLoader.hh"
#include <string>

CollieLoader::CollieLoader()
    : m_channel(NULL),
      m_file(NULL) {}

bool CollieLoader::open(const char* filename, const char* options) {

  // find the name of the channel
  if (options == NULL || strstr(options, "name=") == NULL)
    return false;

  std::string name = strstr(options, "name=");
  std::string::size_type pos = name.find("'");
  if (pos == std::string::npos)
    return false;
  name.erase(0, pos + 1);
  pos = name.find("'");
  if (pos == std::string::npos || pos == 0)
    return false;
  name.erase(pos);
  if (!Loader::open(filename, options))
    return false;
  if (m_file != NULL)
    close();
  m_file = NULL;
  m_channel = NULL;
  m_file = TFile::Open(filename, "");
  if (m_file == NULL || m_file->IsZombie()) {
    if (m_file != NULL) {
      m_file->Close();
      delete m_file;
      m_file = NULL;
    }
    return false;
  }

  m_channel = CollieChannel::loadChannel(name.c_str(), m_file);
  if (m_channel == NULL)
    return false;
  return true;
}

bool CollieLoader::close() {
  if (m_file != NULL) {
    m_file->Close();
    delete m_file;
    m_file = NULL;
    m_channel = NULL;
    return true;
  } else
    return false;
}

void CollieLoader::setbkgdscales(int ns, double* sc) {
  nscales = ns;
  for (int in = 0; in < ns; ++in) {
    scales[in] = sc[in];
  }
}

SigBkgdDist* CollieLoader::MassPoint2SBD(CollieMasspoint* mp) {

  if (mp == NULL)
    return NULL;
  // get the shape information from the distribution
  int nbins;
  double min = 0, max = 1;

  if (mp->getNYbins() > 1) {
    nbins = mp->getNYbins() * mp->getNXbins();
    min = 0.0;
    max = 1.0;
  } else {
    nbins = mp->getNXbins();
    min = mp->getMinX();
    max = mp->getMaxX();
  }

  SigBkgdDist* sbd = new SigBkgdDist(m_channel->getChannelName(), (int)nbins, min, max,
                                     mp->getVar1(), mp->getVar2(), mp->getVar3());

  // important questions: what model rate do we apply?
  //Should be set by options, but for now we just take the first rate for everything...
  //  double lumi=m_channel->getLuminosity();

  for (int n = 0; n < mp->getNSignalDists(); ++n) {
    CollieDistribution* dist = new CollieDistribution(*(mp->getSignalDist(n)));
    if (dist == NULL)
      return NULL;
    sbd->addSigDist(m_channel->getChannelName(), m_channel->getSignalName(n), dist);
  }

  double tot = 0.0;
  for (int n = 0; n < mp->getNBkgdDists(); ++n) {
    CollieDistribution* dist = new CollieDistribution(*(mp->getBkgdDist(n)));
    if (dist == NULL)
      return NULL;
    tot += dist->sumEfficiency();
    sbd->addBkgdDist(m_channel->getChannelName(), dist);
  }

  CollieDistribution* dist = new CollieDistribution(*(mp->getDataDist()));
  if (dist == NULL)
    return NULL;
  sbd->addDataDist(m_channel->getChannelName(), dist);

  sbd->setBaselineModel();
  sbd->fillArrays(true);

  return sbd;
}

SigBkgdDist* CollieLoader::get(int var1) {
  if (m_channel == NULL)
    return NULL;
  CollieMasspoint* mp = m_channel->getMasspoint(var1);
  if (mp == NULL)
    return NULL;
  return MassPoint2SBD(mp);
}

SigBkgdDist* CollieLoader::get(int var1, int var2) {
  if (m_channel == NULL)
    return NULL;
  CollieMasspoint* mp = m_channel->getMasspoint(var1, var2);
  if (mp == NULL)
    return NULL;
  return MassPoint2SBD(mp);
}

SigBkgdDist* CollieLoader::get(int var1, int var2, int var3) {
  if (m_channel == NULL)
    return NULL;
  CollieMasspoint* mp = m_channel->getMasspoint(var1, var2, var3);
  if (mp == NULL)
    return NULL;
  return MassPoint2SBD(mp);
}

int CollieLoader::getMasspointList(int nPoints, int* v1, int* v2, int* v3) {
  if (m_channel == NULL || v1 == NULL || nPoints < 0)
    return false;
  int n = getNMasspoints();
  if (nPoints < n || n == -1)
    return false;
  if (v2 == NULL && m_channel->getNIndepVariables() > 1)
    return false;
  if (v3 == NULL && m_channel->getNIndepVariables() > 2)
    return false;

  // all checked.  Now we must assume the arrays to be otherwise ok.
  CollieIterator* i = m_channel->getIterator();
  int j = 0;
  while (i->hasNext()) {
    i->nextNoMasspoint();
    v1[j] = i->getVar1();
    if (v2 != NULL)
      v2[j] = i->getVar2();
    if (v3 != NULL)
      v3[j] = i->getVar3();
    ++j;
  }
  delete i;
  return true;
}

int CollieLoader::getNMasspoints() {
  return (m_channel == NULL) ? -1 : m_channel->getNMasspoints();
}
