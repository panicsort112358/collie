#include <TFile.h>
#include <TKey.h>
#include <TNamed.h>
#include <TTree.h>

#include <iostream>
#include <set>
#include <string>

int main(int argc, char* argv[]) {
  if (argc != 3)
    return 2;
  TFile input(argv[1], "READ");
  if (input.IsZombie()) {
    std::cerr << "Cannot open result file\n";
    return 3;
  }

  std::set<std::string> keys;
  TIter next(input.GetListOfKeys());
  while (TKey* key = static_cast<TKey*>(next()))
    keys.insert(key->GetName());
  if (keys != std::set<std::string>{"SCAN", "collie_run_manifest"}) {
    std::cerr << "Unexpected top-level result objects\n";
    return 4;
  }
  TTree* tree = dynamic_cast<TTree*>(input.Get("SCAN"));
  TNamed* manifest = dynamic_cast<TNamed*>(input.Get("collie_run_manifest"));
  if (tree == nullptr || tree->GetEntries() <= 0 || manifest == nullptr) {
    std::cerr << "Incomplete result file\n";
    return 5;
  }
  const std::string expected = "\"command\":\"" + std::string(argv[2]) + "\"";
  if (std::string(manifest->GetTitle()).find(expected) == std::string::npos) {
    std::cerr << "Unexpected manifest command\n";
    return 6;
  }
  return 0;
}
