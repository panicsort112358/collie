#ifndef COLLIE_SIG_BKGD_DIST_ADAPTER_HH
#define COLLIE_SIG_BKGD_DIST_ADAPTER_HH

#include "BinnedModel.hh"

#include <string>

class SigBkgdDist;

namespace collie {

BinnedModel makeBinnedModel(const SigBkgdDist& distribution,
                            const std::string& fallbackChannelName = "combined");

} // namespace collie

#endif
