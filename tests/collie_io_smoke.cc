#include "CollieDistribution.hh"

#include <TClass.h>

#include <iostream>
#include <string>

int main() {
  CollieDistribution distribution;
  if (std::string(distribution.GetName()) != "Default") {
    std::cerr << "Unexpected default distribution name\n";
    return 1;
  }

  const TClass* dictionary = TClass::GetClass("CollieDistribution");
  if (dictionary == nullptr) {
    std::cerr << "ROOT dictionary for CollieDistribution was not loaded\n";
    return 2;
  }
  if (dictionary->GetClassVersion() != 1) {
    std::cerr << "Unexpected CollieDistribution dictionary version\n";
    return 3;
  }

  return 0;
}
