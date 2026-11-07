#include "CollieCli.hh"

#include "BinnedModel.hh"
#include "CLcompute.hh"
#include "CLfast.hh"
#include "CLfit.hh"
#include "CLfit2.hh"
#include "CLpoint.hh"
#include "CLsyst.hh"
#include "CollieChannel.hh"
#include "CollieLoader.hh"
#include "CrossSectionCalc.hh"
#include "CrossSectionLimit.hh"
#include "RunContext.hh"
#include "SigBkgdDist.hh"
#include "SigBkgdDistAdapter.hh"

#include <TBranch.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TList.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>
#include <vector>

namespace collie {
namespace cli {
namespace {

class CliError : public std::runtime_error {
public:
  explicit CliError(const std::string& message)
      : std::runtime_error(message) {}
};

class HistogramAttachmentGuard {
public:
  explicit HistogramAttachmentGuard(bool enabled)
      : previous_(TH1::AddDirectoryStatus()) {
    TH1::AddDirectory(enabled);
  }

  ~HistogramAttachmentGuard() {
    TH1::AddDirectory(previous_);
  }

private:
  bool previous_;
};

std::string pointString(const Point& point) {
  std::ostringstream output;
  output << point.var1 << '/' << point.var2 << '/' << point.var3;
  return output.str();
}

std::string trim(const std::string& value) {
  const std::string whitespace = " \t\r\n";
  const std::size_t first = value.find_first_not_of(whitespace);
  if (first == std::string::npos)
    return "";
  const std::size_t last = value.find_last_not_of(whitespace);
  return value.substr(first, last - first + 1);
}

bool fileExists(const std::string& path) {
  struct stat status;
  return ::stat(path.c_str(), &status) == 0;
}

std::string requireValue(const std::vector<std::string>& arguments, std::size_t& index,
                         const std::string& option) {
  if (++index >= arguments.size()) {
    throw CliError(option + " requires a value");
  }
  return arguments[index];
}

long long parseSigned(const std::string& value, const std::string& option) {
  std::size_t parsed = 0;
  try {
    const long long result = std::stoll(value, &parsed, 10);
    if (parsed != value.size())
      throw std::invalid_argument("trailing");
    return result;
  } catch (const std::exception&) {
    throw CliError(option + " expects an integer, got '" + value + "'");
  }
}

std::uint64_t parseUnsigned(const std::string& value, const std::string& option) {
  if (value.empty() || value[0] == '-') {
    throw CliError(option + " expects a non-negative integer");
  }
  std::size_t parsed = 0;
  try {
    const unsigned long long result = std::stoull(value, &parsed, 10);
    if (parsed != value.size())
      throw std::invalid_argument("trailing");
    return static_cast<std::uint64_t>(result);
  } catch (const std::exception&) {
    throw CliError(option + " expects a non-negative integer, got '" + value + "'");
  }
}

double parseDouble(const std::string& value, const std::string& option) {
  std::size_t parsed = 0;
  try {
    const double result = std::stod(value, &parsed);
    if (parsed != value.size() || !std::isfinite(result)) {
      throw std::invalid_argument("invalid");
    }
    return result;
  } catch (const std::exception&) {
    throw CliError(option + " expects a finite number, got '" + value + "'");
  }
}

int checkedInt(long long value, const std::string& option) {
  if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max()) {
    throw CliError(option + " is outside the supported integer range");
  }
  return static_cast<int>(value);
}

Point parsePoint(std::string value, int& dimensions) {
  std::replace(value.begin(), value.end(), ':', ',');
  std::vector<std::string> fields;
  std::size_t begin = 0;
  while (begin <= value.size()) {
    const std::size_t end = value.find(',', begin);
    fields.push_back(value.substr(begin, end - begin));
    if (end == std::string::npos)
      break;
    begin = end + 1;
  }
  if (fields.empty() || fields.size() > 3) {
    throw CliError("--point expects V1, V1,V2, or V1,V2,V3");
  }
  dimensions = static_cast<int>(fields.size());
  Point point;
  point.var1 = checkedInt(parseSigned(fields[0], "--point"), "--point");
  if (fields.size() > 1) {
    point.var2 = checkedInt(parseSigned(fields[1], "--point"), "--point");
  }
  if (fields.size() > 2) {
    point.var3 = checkedInt(parseSigned(fields[2], "--point"), "--point");
  }
  return point;
}

void appendInputList(const std::string& path, std::vector<std::string>& output) {
  std::ifstream input(path.c_str());
  if (!input)
    throw CliError("cannot open input list '" + path + "'");

  std::string line;
  while (std::getline(input, line)) {
    const std::size_t comment = line.find('#');
    if (comment != std::string::npos)
      line.erase(comment);
    line = trim(line);
    if (!line.empty())
      output.push_back(line);
  }
}

std::vector<std::string> normalizedArguments(int argc, char* argv[]) {
  std::vector<std::string> output;
  for (int i = 0; i < argc; ++i) {
    const std::string value(argv[i]);
    const std::size_t equals = value.find('=');
    if (i > 1 && value.compare(0, 2, "--") == 0 && equals != std::string::npos) {
      output.push_back(value.substr(0, equals));
      output.push_back(value.substr(equals + 1));
    } else {
      output.push_back(value);
    }
  }
  return output;
}

bool parseRunOption(const std::vector<std::string>& arguments, std::size_t& index,
                    RunConfig& config) {
  const std::string& option = arguments[index];
  if (option == "--seed") {
    config.seed = parseUnsigned(requireValue(arguments, index, option), option);
    config.seedExplicit = true;
    return true;
  }
  if (option == "--workers" || option == "-j") {
    const std::uint64_t workers = parseUnsigned(requireValue(arguments, index, option), option);
    if (workers > std::numeric_limits<std::size_t>::max()) {
      throw CliError(option + " is outside the supported range");
    }
    if (workers > RunContext::maximumWorkers()) {
      throw CliError(option + " must not exceed " + std::to_string(RunContext::maximumWorkers()));
    }
    config.workers = static_cast<std::size_t>(workers);
    config.workersExplicit = true;
    return true;
  }
  if (option == "--verbose" || option == "-v") {
    config.verbose = true;
    return true;
  }
  if (option == "--force" || option == "-f") {
    config.force = true;
    return true;
  }
  return false;
}

bool parseInputOption(const std::vector<std::string>& arguments, std::size_t& index,
                      InputConfig& config) {
  const std::string& option = arguments[index];
  if (option == "--point" || option == "-p") {
    config.selection.all = false;
    config.selection.point =
        parsePoint(requireValue(arguments, index, option), config.selection.dimensions);
    return true;
  }
  if (option == "--input-list" || option == "-L") {
    appendInputList(requireValue(arguments, index, option), config.files);
    return true;
  }
  return false;
}

void rejectUnknownOption(const std::string& option) {
  throw CliError("unknown option '" + option + "'");
}

void requireInputs(const std::vector<std::string>& files, const std::string& command) {
  if (files.empty())
    throw CliError(command + " requires at least one input");
}

void requireOutput(const std::string& output, const std::string& command) {
  if (output.empty())
    throw CliError(command + " requires --output FILE");
}

LimitMethod parseLimitMethod(const std::string& value) {
  if (value == "fast")
    return LimitMethod::Fast;
  if (value == "syst" || value == "systematics") {
    return LimitMethod::Systematics;
  }
  if (value == "profile" || value == "fit")
    return LimitMethod::Profile;
  if (value == "profile2" || value == "fit2")
    return LimitMethod::Profile2;
  throw CliError("--method expects fast, systematics, profile, or profile2");
}

std::string limitMethodName(LimitMethod method) {
  switch (method) {
    case LimitMethod::Fast:
      return "fast";
    case LimitMethod::Systematics:
      return "systematics";
    case LimitMethod::Profile:
      return "profile";
    case LimitMethod::Profile2:
      return "profile2";
  }
  return "unknown";
}

InspectConfig parseInspect(const std::vector<std::string>& arguments) {
  InspectConfig config;
  for (std::size_t i = 2; i < arguments.size(); ++i) {
    const std::string& option = arguments[i];
    if (option == "--details" || option == "-d") {
      config.details = true;
    } else if (option == "--input-list" || option == "-L") {
      appendInputList(requireValue(arguments, i, option), config.files);
    } else if (!option.empty() && option[0] == '-') {
      rejectUnknownOption(option);
    } else {
      config.files.push_back(option);
    }
  }
  requireInputs(config.files, "inspect");
  return config;
}

ValidateConfig parseValidate(const std::vector<std::string>& arguments) {
  ValidateConfig config;
  for (std::size_t i = 2; i < arguments.size(); ++i) {
    const std::string& option = arguments[i];
    if (option == "--strict") {
      config.strict = true;
    } else if (option == "--input-list" || option == "-L") {
      appendInputList(requireValue(arguments, i, option), config.files);
    } else if (!option.empty() && option[0] == '-') {
      rejectUnknownOption(option);
    } else {
      config.files.push_back(option);
    }
  }
  requireInputs(config.files, "validate");
  return config;
}

LimitConfig parseLimit(const std::vector<std::string>& arguments) {
  LimitConfig config;
  for (std::size_t i = 2; i < arguments.size(); ++i) {
    const std::string& option = arguments[i];
    if (parseRunOption(arguments, i, config.run) || parseInputOption(arguments, i, config.input)) {
      continue;
    }
    if (option == "--output" || option == "-o") {
      config.output = requireValue(arguments, i, option);
    } else if (option == "--method" || option == "-m") {
      config.method = parseLimitMethod(requireValue(arguments, i, option));
    } else if (option == "--precision") {
      config.precision =
          checkedInt(parseSigned(requireValue(arguments, i, option), option), option);
    } else if (option == "--iterations") {
      config.iterations =
          checkedInt(parseSigned(requireValue(arguments, i, option), option), option);
    } else if (option == "--cl") {
      config.confidenceLevel = parseDouble(requireValue(arguments, i, option), option);
    } else if (option == "--accuracy") {
      config.accuracy = parseDouble(requireValue(arguments, i, option), option);
    } else if (option == "--search-start") {
      config.searchStart = parseDouble(requireValue(arguments, i, option), option);
    } else if (option == "--sigma") {
      config.sigma = checkedInt(parseSigned(requireValue(arguments, i, option), option), option);
    } else if (option == "--signal-exclusion") {
      config.signalExclusion = parseDouble(requireValue(arguments, i, option), option);
    } else if (option == "--expected-only") {
      config.calculateExpected = true;
      config.calculateObserved = false;
    } else if (option == "--observed-only") {
      config.calculateExpected = false;
      config.calculateObserved = true;
    } else if (option == "--histogram-statistics") {
      config.useHistogramStatistics = true;
    } else if (option == "--fit-signal") {
      config.fitSignal = true;
    } else if (!option.empty() && option[0] == '-') {
      rejectUnknownOption(option);
    } else {
      config.input.files.push_back(option);
    }
  }

  requireInputs(config.input.files, "limit");
  requireOutput(config.output, "limit");
  if (config.precision < 0 || config.precision > 4) {
    throw CliError("--precision must be between 0 and 4");
  }
  if (config.iterations < 0) {
    throw CliError("--iterations must be non-negative");
  }
  if (config.iterations > CLcompute::MAX_MONTE_CARLO_ITERATIONS) {
    throw CliError("--iterations exceeds the safe trial-count limit");
  }
  if (!(config.confidenceLevel > 0.5 && config.confidenceLevel < 1.0)) {
    throw CliError("--cl must be greater than 0.5 and less than 1");
  }
  if (!(config.accuracy > 0.0 && config.accuracy < 1.0 - config.confidenceLevel)) {
    throw CliError("--accuracy must be positive and smaller than 1 - CL");
  }
  if (config.searchStart <= 0.0) {
    throw CliError("--search-start must be positive");
  }
  if (config.sigma < -2 || config.sigma > 2) {
    throw CliError("--sigma must be between -2 and 2");
  }
  if (config.signalExclusion < 0.0) {
    throw CliError("--signal-exclusion must be non-negative");
  }
  return config;
}

DataHypothesis parseHypothesis(const std::string& value) {
  if (value == "observed" || value == "data") {
    return DataHypothesis::Observed;
  }
  if (value == "background" || value == "b-only") {
    return DataHypothesis::BackgroundOnly;
  }
  if (value == "signal+background" || value == "s+b") {
    return DataHypothesis::SignalPlusBackground;
  }
  throw CliError("--hypothesis expects observed, background, or signal+background");
}

std::string hypothesisName(DataHypothesis value) {
  switch (value) {
    case DataHypothesis::Observed:
      return "observed";
    case DataHypothesis::BackgroundOnly:
      return "background";
    case DataHypothesis::SignalPlusBackground:
      return "signal+background";
  }
  return "unknown";
}

FitConfig parseFit(const std::vector<std::string>& arguments) {
  FitConfig config;
  for (std::size_t i = 2; i < arguments.size(); ++i) {
    const std::string& option = arguments[i];
    if (parseRunOption(arguments, i, config.run) || parseInputOption(arguments, i, config.input)) {
      continue;
    }
    if (option == "--output" || option == "-o") {
      config.output = requireValue(arguments, i, option);
    } else if (option == "--precision") {
      config.precision =
          checkedInt(parseSigned(requireValue(arguments, i, option), option), option);
    } else if (option == "--lower") {
      config.lowerBound = parseDouble(requireValue(arguments, i, option), option);
    } else if (option == "--upper") {
      config.upperBound = parseDouble(requireValue(arguments, i, option), option);
    } else if (option == "--step") {
      config.granularity = parseDouble(requireValue(arguments, i, option), option);
    } else if (option == "--scan") {
      config.scan = true;
    } else if (option == "--plots") {
      config.makePlots = true;
    } else if (option == "--pseudo-experiments") {
      config.pseudoExperiments =
          checkedInt(parseSigned(requireValue(arguments, i, option), option), option);
    } else if (option == "--pseudo-signal") {
      config.pseudoSignal = parseDouble(requireValue(arguments, i, option), option);
    } else if (option == "--hypothesis") {
      config.hypothesis = parseHypothesis(requireValue(arguments, i, option));
    } else if (!option.empty() && option[0] == '-') {
      rejectUnknownOption(option);
    } else {
      config.input.files.push_back(option);
    }
  }

  requireInputs(config.input.files, "fit");
  requireOutput(config.output, "fit");
  if (config.precision < 0 || config.precision > 4) {
    throw CliError("--precision must be between 0 and 4");
  }
  if (config.granularity <= 0.0) {
    throw CliError("--step must be positive");
  }
  if (config.upperBound <= config.lowerBound) {
    throw CliError("--upper must be greater than --lower");
  }
  if (config.pseudoExperiments < 0) {
    throw CliError("--pseudo-experiments must be non-negative");
  }
  return config;
}

CombineConfig parseCombine(const std::vector<std::string>& arguments) {
  CombineConfig config;
  for (std::size_t i = 2; i < arguments.size(); ++i) {
    const std::string& option = arguments[i];
    if (parseRunOption(arguments, i, config.run))
      continue;
    if (option == "--output" || option == "-o") {
      config.output = requireValue(arguments, i, option);
    } else if (option == "--input-list" || option == "-L") {
      appendInputList(requireValue(arguments, i, option), config.inputs);
    } else if (option == "--allow-duplicate-points") {
      config.allowDuplicatePoints = true;
    } else if (!option.empty() && option[0] == '-') {
      rejectUnknownOption(option);
    } else {
      config.inputs.push_back(option);
    }
  }
  requireInputs(config.inputs, "combine");
  requireOutput(config.output, "combine");
  return config;
}

void printGeneralHelp(std::ostream& output) {
  output << "Collie statistical inference\n\n"
         << "Usage: collie COMMAND [OPTIONS] INPUT...\n\n"
         << "Commands:\n"
         << "  inspect    Describe Collie inputs or result files\n"
         << "  limit      Calculate expected and observed confidence limits\n"
         << "  fit        Fit a signal-strength model\n"
         << "  combine    Concatenate compatible SCAN result trees\n"
         << "  validate   Check file structure and numerical model invariants\n\n"
         << "Run 'collie COMMAND --help' for command-specific options.\n";
}

void printCommandHelp(const std::string& command, std::ostream& output) {
  if (command == "inspect") {
    output << "Usage: collie inspect [--details] INPUT...\n\n"
           << "  -d, --details          Print every mass point and nominal yields\n"
           << "  -L, --input-list FILE  Read one input path per line\n";
  } else if (command == "validate") {
    output << "Usage: collie validate [--strict] INPUT...\n\n"
           << "      --strict           Treat numerical warnings as errors\n"
           << "  -L, --input-list FILE  Read one input path per line\n";
  } else if (command == "limit") {
    output << "Usage: collie limit -o RESULT.root [OPTIONS] INPUT...\n\n"
           << "  -o, --output FILE         Atomic ROOT output (required)\n"
           << "  -p, --point V1[,V2,V3]    Calculate one exact mass point\n"
           << "  -m, --method METHOD       fast, systematics, profile, profile2\n"
           << "      --precision N         Effort from 0 (fast) to 4 (fine)\n"
           << "      --iterations N        Override Monte Carlo trial count\n"
           << "      --cl VALUE            Confidence level (default 0.95)\n"
           << "      --accuracy VALUE      Limit-search tolerance\n"
           << "      --search-start VALUE  Initial signal multiplier\n"
           << "      --sigma N             Expected band from -2 to +2\n"
           << "      --expected-only       Skip the observed limit\n"
           << "      --observed-only       Skip the expected limit\n"
           << "      --seed N              Reproducible master RNG seed\n"
           << "  -j, --workers N           Bounded worker count (0 = automatic)\n"
           << "  -L, --input-list FILE     Read one input path per line\n"
           << "  -f, --force               Replace an existing output\n"
           << "  -v, --verbose             Print inference progress\n\n"
           << "The 'fast' method ignores nuisance uncertainties and is intended\n"
           << "for input checks, not publication results.\n";
  } else if (command == "fit") {
    output << "Usage: collie fit -o RESULT.root [OPTIONS] INPUT...\n\n"
           << "  -o, --output FILE            Atomic ROOT output (required)\n"
           << "  -p, --point V1[,V2,V3]       Fit one exact mass point\n"
           << "      --hypothesis TYPE        observed, background, signal+background\n"
           << "      --scan                   Also scan confidence levels\n"
           << "      --lower VALUE            Scan lower bound\n"
           << "      --upper VALUE            Scan upper bound\n"
           << "      --step VALUE             Scan granularity\n"
           << "      --precision N            Effort from 0 to 4\n"
           << "      --pseudo-experiments N   Run N fit pseudo-experiments\n"
           << "      --pseudo-signal VALUE    Injected pseudo-experiment strength\n"
           << "      --plots                  Persist diagnostic fit plots\n"
           << "      --seed N                 Reproducible master RNG seed\n"
           << "  -j, --workers N              Bounded worker count (0 = automatic)\n"
           << "  -L, --input-list FILE        Read one input path per line\n"
           << "  -f, --force                  Replace an existing output\n"
           << "  -v, --verbose                Print fit details\n";
  } else if (command == "combine") {
    output << "Usage: collie combine -o RESULT.root [OPTIONS] INPUT...\n\n"
           << "  -o, --output FILE              Atomic ROOT output (required)\n"
           << "  -L, --input-list FILE          Read one input path per line\n"
           << "      --allow-duplicate-points   Permit repeated parameter points\n"
           << "  -f, --force                    Replace an existing output\n";
  } else {
    throw CliError("unknown command '" + command + "'");
  }
}

class AtomicRootOutput {
public:
  AtomicRootOutput(const std::string& target, bool force)
      : target_(target),
        force_(force),
        committed_(false) {
    if (!force_ && fileExists(target_)) {
      throw CliError("output '" + target_ + "' already exists (use --force to replace it)");
    }

    std::string pattern = target_ + ".tmp.XXXXXX";
    std::vector<char> buffer(pattern.begin(), pattern.end());
    buffer.push_back('\0');
    const int descriptor = ::mkstemp(buffer.data());
    if (descriptor < 0) {
      throw CliError("cannot create temporary output beside '" + target_ +
                     "': " + std::strerror(errno));
    }
    const mode_t processMask = ::umask(0);
    ::umask(processMask);
    if (::fchmod(descriptor, static_cast<mode_t>(0666 & ~processMask)) != 0) {
      const int failure = errno;
      ::close(descriptor);
      std::remove(buffer.data());
      throw CliError("cannot set permissions on temporary output beside '" + target_ +
                     "': " + std::strerror(failure));
    }
    if (::close(descriptor) != 0) {
      const int failure = errno;
      std::remove(buffer.data());
      throw CliError("cannot close temporary output beside '" + target_ +
                     "': " + std::strerror(failure));
    }
    temporary_ = buffer.data();

    file_.reset(TFile::Open(temporary_.c_str(), "RECREATE"));
    if (!file_ || file_->IsZombie()) {
      std::remove(temporary_.c_str());
      throw CliError("ROOT could not create output '" + temporary_ + "'");
    }
  }

  ~AtomicRootOutput() {
    if (file_) {
      file_->Close();
      file_.reset();
    }
    if (!committed_ && !temporary_.empty()) {
      std::remove(temporary_.c_str());
    }
  }

  TFile& file() {
    return *file_;
  }

  void commit() {
    file_->Write("", TObject::kOverwrite);
    file_->Close();
    file_.reset();
    if (!force_ && fileExists(target_)) {
      throw CliError("output '" + target_ + "' appeared while the calculation was running");
    }
    if (std::rename(temporary_.c_str(), target_.c_str()) != 0) {
      throw CliError("cannot publish output '" + target_ + "': " + std::strerror(errno));
    }
    committed_ = true;
  }

private:
  std::string target_;
  std::string temporary_;
  bool force_;
  bool committed_;
  std::unique_ptr<TFile> file_;
};

std::string jsonEscape(const std::string& value) {
  std::ostringstream output;
  for (std::string::const_iterator character = value.begin(); character != value.end();
       ++character) {
    switch (*character) {
      case '"':
        output << "\\\"";
        break;
      case '\\':
        output << "\\\\";
        break;
      case '\b':
        output << "\\b";
        break;
      case '\f':
        output << "\\f";
        break;
      case '\n':
        output << "\\n";
        break;
      case '\r':
        output << "\\r";
        break;
      case '\t':
        output << "\\t";
        break;
      default:
        if (static_cast<unsigned char>(*character) < 0x20) {
          output << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                 << static_cast<unsigned int>(static_cast<unsigned char>(*character)) << std::dec;
        } else {
          output << *character;
        }
    }
  }
  return output.str();
}

std::string fileChecksum(const std::string& path) {
  std::ifstream input(path.c_str(), std::ios::binary);
  if (!input)
    throw CliError("cannot read input '" + path + "'");

  std::uint64_t hash = UINT64_C(14695981039346656037);
  char buffer[64 * 1024];
  while (input) {
    input.read(buffer, sizeof(buffer));
    const std::streamsize count = input.gcount();
    for (std::streamsize i = 0; i < count; ++i) {
      hash ^= static_cast<unsigned char>(buffer[i]);
      hash *= UINT64_C(1099511628211);
    }
  }
  if (!input.eof())
    throw CliError("failed while checksumming '" + path + "'");

  std::ostringstream output;
  output << std::hex << std::setw(16) << std::setfill('0') << hash;
  return output.str();
}

std::string makeManifest(const std::string& command, const std::vector<std::string>& inputs,
                         const std::string& settings) {
  std::ostringstream output;
  output << "{\"schema\":1,\"command\":\"" << jsonEscape(command) << "\""
         << ",\"run\":" << RunContext::global().manifestJson("collie-" + command)
         << ",\"root_version\":\"" << jsonEscape(gROOT->GetVersion()) << "\""
         << ",\"checksum\":\"fnv1a64\",\"inputs\":[";
  for (std::size_t i = 0; i < inputs.size(); ++i) {
    if (i != 0)
      output << ',';
    output << "{\"path\":\"" << jsonEscape(inputs[i]) << "\",\"digest\":\""
           << fileChecksum(inputs[i]) << "\"}";
  }
  output << "],\"settings\":" << settings << '}';
  return output.str();
}

void writeManifest(TFile& file, const std::string& manifest) {
  file.cd();
  TNamed object("collie_run_manifest", manifest.c_str());
  if (object.Write("", TObject::kOverwrite) <= 0) {
    throw CliError("failed to write the run manifest");
  }
}

std::size_t jsonMemberValue(const std::string& document, const std::string& member,
                            const std::string& path) {
  const std::string key = "\"" + member + "\"";
  const std::size_t keyPosition = document.find(key);
  if (keyPosition == std::string::npos) {
    throw CliError("manifest in '" + path + "' has no '" + member + "' field");
  }
  std::size_t position = keyPosition + key.size();
  while (position < document.size() &&
         std::isspace(static_cast<unsigned char>(document[position]))) {
    ++position;
  }
  if (position >= document.size() || document[position] != ':') {
    throw CliError("manifest in '" + path + "' has malformed '" + member + "' field");
  }
  do {
    ++position;
  } while (position < document.size() &&
           std::isspace(static_cast<unsigned char>(document[position])));
  return position;
}

std::string jsonStringMember(const std::string& document, const std::string& member,
                             const std::string& path) {
  std::size_t position = jsonMemberValue(document, member, path);
  if (position >= document.size() || document[position] != '"') {
    throw CliError("manifest in '" + path + "' has a non-string '" + member + "' field");
  }
  const std::size_t begin = ++position;
  bool escaped = false;
  for (; position < document.size(); ++position) {
    const char value = document[position];
    if (!escaped && value == '"') {
      return document.substr(begin, position - begin);
    }
    if (!escaped && value == '\\') {
      escaped = true;
    } else {
      escaped = false;
    }
  }
  throw CliError("manifest in '" + path + "' has an unterminated '" + member + "' field");
}

std::string jsonObjectMember(const std::string& document, const std::string& member,
                             const std::string& path) {
  const std::size_t begin = jsonMemberValue(document, member, path);
  if (begin >= document.size() || document[begin] != '{') {
    throw CliError("manifest in '" + path + "' has a non-object '" + member + "' field");
  }

  std::size_t depth = 0;
  bool inString = false;
  bool escaped = false;
  for (std::size_t position = begin; position < document.size(); ++position) {
    const char value = document[position];
    if (inString) {
      if (!escaped && value == '"')
        inString = false;
      if (!escaped && value == '\\') {
        escaped = true;
      } else {
        escaped = false;
      }
      continue;
    }
    if (value == '"') {
      inString = true;
    } else if (value == '{') {
      ++depth;
    } else if (value == '}' && --depth == 0) {
      return document.substr(begin, position - begin + 1);
    }
  }
  throw CliError("manifest in '" + path + "' has an unterminated '" + member + "' object");
}

struct ResultProvenance {
  std::string command;
  std::string settings;
};

ResultProvenance resultProvenance(TFile& file, const std::string& path) {
  TNamed* manifest = dynamic_cast<TNamed*>(file.Get("collie_run_manifest"));
  if (manifest == nullptr) {
    throw CliError("'" + path + "' has no collie_run_manifest");
  }
  const std::string document = manifest->GetTitle();
  ResultProvenance provenance;
  provenance.command = jsonStringMember(document, "command", path);
  provenance.settings = jsonObjectMember(document, "settings", path);
  return provenance;
}

std::string discoverChannelName(TFile& file) {
  std::set<std::string> candidates;
  TIter next(file.GetListOfKeys());
  while (TKey* key = static_cast<TKey*>(next())) {
    const std::string name = key->GetName();
    if (name != "SCAN" && name != "collie_run_manifest") {
      candidates.insert(name);
    }
  }
  if (candidates.size() != 1) {
    std::ostringstream message;
    message << "expected exactly one Collie channel in '" << file.GetName() << "', found "
            << candidates.size();
    throw CliError(message.str());
  }
  return *candidates.begin();
}

std::string discoverChannelName(const std::string& path) {
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if (!file || file->IsZombie()) {
    throw CliError("ROOT could not open input '" + path + "'");
  }
  return discoverChannelName(*file);
}

std::vector<Point> loaderPoints(CollieLoader& loader, const std::string& path) {
  const int count = loader.getNMasspoints();
  if (count <= 0) {
    throw CliError("'" + path + "' contains no mass points");
  }
  std::vector<int> first(count);
  std::vector<int> second(count);
  std::vector<int> third(count);
  if (!loader.getMasspointList(count, first.data(), second.data(), third.data())) {
    throw CliError("could not read mass-point indices from '" + path + "'");
  }

  std::vector<Point> points;
  points.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    Point point;
    point.var1 = first[i];
    point.var2 = second[i];
    point.var3 = third[i];
    points.push_back(point);
  }
  std::sort(points.begin(), points.end());
  if (std::adjacent_find(points.begin(), points.end()) != points.end()) {
    throw CliError("'" + path + "' contains duplicate mass-point indices");
  }
  return points;
}

struct LoadedInput {
  std::string path;
  std::string channel;
  std::unique_ptr<CollieLoader> loader;
  std::vector<Point> points;
};

class LoadedInputs {
public:
  ~LoadedInputs() {
    for (std::vector<LoadedInput>::iterator input = inputs.begin(); input != inputs.end();
         ++input) {
      if (input->loader)
        input->loader->close();
    }
  }

  std::vector<LoadedInput> inputs;
  std::vector<Point> points;
};

std::unique_ptr<LoadedInputs> loadInputs(const std::vector<std::string>& paths) {
  std::unique_ptr<LoadedInputs> loaded(new LoadedInputs);
  loaded->inputs.reserve(paths.size());

  std::set<std::string> channels;
  for (std::vector<std::string>::const_iterator path = paths.begin(); path != paths.end(); ++path) {
    LoadedInput input;
    input.path = *path;
    input.channel = discoverChannelName(*path);
    if (!channels.insert(input.channel).second) {
      throw CliError("channel '" + input.channel + "' appears more than once in the input set");
    }

    input.loader.reset(new CollieLoader);
    const std::string options = "name='" + input.channel + "'";
    if (!input.loader->open(input.path.c_str(), options.c_str())) {
      throw CliError("Collie could not load channel '" + input.channel + "' from '" + input.path +
                     "'");
    }
    input.points = loaderPoints(*input.loader, input.path);
    if (loaded->inputs.empty()) {
      loaded->points = input.points;
    } else if (input.points != loaded->points) {
      throw CliError("mass-point grid in '" + input.path +
                     "' does not exactly match the first channel");
    }
    loaded->inputs.push_back(std::move(input));
  }
  return loaded;
}

std::vector<Point> selectedPoints(const LoadedInputs& inputs, const PointSelection& selection) {
  if (selection.all)
    return inputs.points;
  std::vector<Point> matches;
  for (std::vector<Point>::const_iterator point = inputs.points.begin();
       point != inputs.points.end(); ++point) {
    const bool firstMatches = point->var1 == selection.point.var1;
    const bool secondMatches = selection.dimensions < 2 || point->var2 == selection.point.var2;
    const bool thirdMatches = selection.dimensions < 3 || point->var3 == selection.point.var3;
    if (firstMatches && secondMatches && thirdMatches)
      matches.push_back(*point);
  }
  if (matches.empty()) {
    throw CliError("requested point " + pointString(selection.point) +
                   " is not present in the input grid");
  }
  if (matches.size() > 1) {
    throw CliError("requested point is ambiguous; provide all independent "
                   "variable coordinates");
  }
  return matches;
}

std::unique_ptr<SigBkgdDist> combinedDistribution(LoadedInputs& inputs, const Point& point) {
  std::unique_ptr<SigBkgdDist> combined(
      inputs.inputs[0].loader->get(point.var1, point.var2, point.var3));
  if (!combined) {
    throw CliError("failed to load point " + pointString(point) + " from '" +
                   inputs.inputs[0].path + "'");
  }
  for (std::size_t i = 1; i < inputs.inputs.size(); ++i) {
    std::unique_ptr<SigBkgdDist> channel(
        inputs.inputs[i].loader->get(point.var1, point.var2, point.var3));
    if (!channel) {
      throw CliError("failed to load point " + pointString(point) + " from '" +
                     inputs.inputs[i].path + "'");
    }
    if (combined->append(*channel) <= 0) {
      throw CliError("failed to append channel '" + inputs.inputs[i].channel + "' at point " +
                     pointString(point));
    }
  }
  return combined;
}

std::vector<std::string> distributionProblems(const SigBkgdDist& distribution,
                                              bool includeWarnings) {
  std::vector<std::string> problems;
  if (distribution.nbins() <= 0) {
    problems.push_back("distribution has no bins");
    return problems;
  }

  const BinnedModel model = makeBinnedModel(distribution);
  const std::vector<std::string> structural = model.validate();
  problems.insert(problems.end(), structural.begin(), structural.end());

  for (int bin = 0; bin < distribution.nbins(); ++bin) {
    const double values[] = {distribution.signal(bin), distribution.bkgd(bin),
                             distribution.data(bin), distribution.getSignalStatErr(bin),
                             distribution.getBkgdStatErr(bin)};
    for (std::size_t value = 0; value < sizeof(values) / sizeof(values[0]); ++value) {
      if (!std::isfinite(values[value])) {
        std::ostringstream message;
        message << "bin " << bin << " contains a non-finite value";
        problems.push_back(message.str());
        break;
      }
    }
    if (includeWarnings && (distribution.signal(bin) < 0.0 || distribution.bkgd(bin) < 0.0 ||
                            distribution.data(bin) < 0.0)) {
      std::ostringstream message;
      message << "bin " << bin << " contains a negative nominal yield";
      problems.push_back(message.str());
    }
    if (includeWarnings &&
        (distribution.getSignalStatErr(bin) < 0.0 || distribution.getBkgdStatErr(bin) < 0.0)) {
      std::ostringstream message;
      message << "bin " << bin << " contains a negative statistical error";
      problems.push_back(message.str());
    }
  }
  if (includeWarnings && distribution.totSignal() <= 0.0) {
    problems.push_back("total signal yield is not positive");
  }
  return problems;
}

int effortForPrecision(int precision) {
  switch (precision) {
    case 0:
      return CLcompute::LEVEL_VERYFAST;
    case 1:
      return CLcompute::LEVEL_FAST;
    case 2:
      return CLcompute::LEVEL_STANDARD;
    case 3:
      return CLcompute::LEVEL_FINE;
    case 4:
      return CLcompute::LEVEL_VERYFINE;
    default:
      throw CliError("internal precision error");
  }
}

std::unique_ptr<CLcompute> makeCalculator(const LimitConfig& config) {
  std::unique_ptr<CLcompute> calculator;
  switch (config.method) {
    case LimitMethod::Fast:
      calculator.reset(new CLfast);
      break;
    case LimitMethod::Systematics:
      calculator.reset(new CLsyst);
      break;
    case LimitMethod::Profile:
      calculator.reset(new CLfit);
      break;
    case LimitMethod::Profile2:
      calculator.reset(new CLfit2);
      break;
  }
  calculator->useHistoStats(config.useHistogramStatistics);
  if (config.iterations > 0) {
    std::ostringstream options;
    options << "its=" << config.iterations;
    calculator->configure(options.str().c_str());
  }
  if (config.method == LimitMethod::Profile) {
    CLfit* profile = static_cast<CLfit*>(calculator.get());
    profile->fitSignal(config.fitSignal);
    profile->logSigExclusion(config.signalExclusion);
  }
  return calculator;
}

void configureRun(const RunConfig& config) {
  RunContext& context = RunContext::global();
  if (config.seedExplicit)
    context.setSeed(config.seed);
  if (config.workersExplicit)
    context.setWorkers(config.workers);
  gROOT->SetBatch(kTRUE);
}

std::string limitSettings(const LimitConfig& config) {
  std::ostringstream output;
  output << "{\"method\":\"" << limitMethodName(config.method) << "\""
         << ",\"precision\":" << config.precision << ",\"iterations\":" << config.iterations
         << ",\"confidence_level\":" << std::setprecision(17) << config.confidenceLevel
         << ",\"accuracy\":" << config.accuracy << ",\"search_start\":" << config.searchStart
         << ",\"sigma\":" << config.sigma
         << ",\"expected\":" << (config.calculateExpected ? "true" : "false")
         << ",\"observed\":" << (config.calculateObserved ? "true" : "false")
         << ",\"histogram_statistics\":" << (config.useHistogramStatistics ? "true" : "false")
         << ",\"fit_signal\":" << (config.fitSignal ? "true" : "false")
         << ",\"signal_exclusion\":" << config.signalExclusion << '}';
  return output.str();
}

std::string fitSettings(const FitConfig& config) {
  std::ostringstream output;
  output << "{\"hypothesis\":\"" << hypothesisName(config.hypothesis) << "\""
         << ",\"precision\":" << config.precision
         << ",\"scan\":" << (config.scan ? "true" : "false")
         << ",\"plots\":" << (config.makePlots ? "true" : "false")
         << ",\"lower\":" << std::setprecision(17) << config.lowerBound
         << ",\"upper\":" << config.upperBound << ",\"step\":" << config.granularity
         << ",\"pseudo_experiments\":" << config.pseudoExperiments
         << ",\"pseudo_signal\":" << config.pseudoSignal << '}';
  return output.str();
}

int executeLimit(const LimitConfig& config) {
  configureRun(config.run);
  std::unique_ptr<LoadedInputs> inputs = loadInputs(config.input.files);
  const std::vector<Point> points = selectedPoints(*inputs, config.input.selection);
  std::unique_ptr<CLcompute> calculator = makeCalculator(config);

  CrossSectionLimit limit;
  limit.setup(calculator.get());
  limit.setCLlevel(config.confidenceLevel);
  limit.setAccuracy(config.accuracy);
  limit.setPrecision(config.precision);
  limit.setNSigma(config.sigma);
  limit.setSearchSeed(config.searchStart);
  limit.calculateExpected(config.calculateExpected);
  limit.calculateObserved(config.calculateObserved);
  limit.setVerbose(config.run.verbose);

  AtomicRootOutput output(config.output, config.run.force);
  output.file().cd();
  TTree tree("SCAN", "Collie confidence-limit scan");
  CLpoint result;
  result.branch(&tree);

  for (std::vector<Point>::const_iterator point = points.begin(); point != points.end(); ++point) {
    std::unique_ptr<SigBkgdDist> distribution = combinedDistribution(*inputs, *point);
    const std::vector<std::string> problems = distributionProblems(*distribution, false);
    if (!problems.empty()) {
      throw CliError("invalid model at point " + pointString(*point) + ": " + problems.front());
    }

    result.reset(point->var1, point->var2, point->var3);
    if (!calculator->calculateCLs(*distribution, result, effortForPrecision(config.precision))) {
      throw CliError("confidence calculation failed at point " + pointString(*point));
    }
    if (!limit.calculate(*distribution, result)) {
      throw CliError("limit search failed at point " + pointString(*point));
    }
    tree.Fill();
    std::cout << pointString(*point) << "  expected=" << result.xsec_medfactor
              << "  observed=" << result.xsec_obsfactor << '\n';
  }

  output.file().cd();
  if (tree.Write("", TObject::kOverwrite) <= 0) {
    throw CliError("failed to write the SCAN tree");
  }
  tree.SetDirectory(NULL);
  writeManifest(output.file(), makeManifest("limit", config.input.files, limitSettings(config)));
  output.commit();
  std::cout << "Wrote " << points.size() << " limit result" << (points.size() == 1 ? "" : "s")
            << " to " << config.output << '\n';
  return 0;
}

int executeFit(const FitConfig& config) {
  configureRun(config.run);
  std::unique_ptr<LoadedInputs> inputs = loadInputs(config.input.files);
  const std::vector<Point> points = selectedPoints(*inputs, config.input.selection);
  std::unique_ptr<CLcompute> calculator(new CLfast);

  CrossSectionCalc fit;
  fit.setup(calculator.get());
  fit.setPrecision(config.precision);
  fit.setLowerBound(config.lowerBound);
  fit.setUpperBound(config.upperBound);
  fit.setGranularity(config.granularity);
  fit.runPseudoExp(config.scan);
  fit.setVerbose(config.run.verbose);
  if (config.hypothesis == DataHypothesis::BackgroundOnly) {
    fit.setBOnlyHypothesis(true);
  } else if (config.hypothesis == DataHypothesis::SignalPlusBackground) {
    fit.setSBHypothesis(true);
  }

  AtomicRootOutput output(config.output, config.run.force);
  output.file().cd();
  TTree tree("SCAN", "Collie signal-strength fits");
  CLpoint result;
  result.branch(&tree);

  for (std::vector<Point>::const_iterator point = points.begin(); point != points.end(); ++point) {
    std::unique_ptr<SigBkgdDist> distribution = combinedDistribution(*inputs, *point);
    const std::vector<std::string> problems = distributionProblems(*distribution, false);
    if (!problems.empty()) {
      throw CliError("invalid model at point " + pointString(*point) + ": " + problems.front());
    }

    result.reset(point->var1, point->var2, point->var3);
    const auto calculate = [&]() {
      if (!fit.calculate(*distribution, result, config.makePlots))
        return 1;
      if (config.pseudoExperiments > 0 &&
          !fit.testFitPE(*distribution, config.pseudoSignal, config.pseudoExperiments)) {
        return 2;
      }
      return 0;
    };

    int fitStatus = 0;
    const bool persistDiagnostics = config.makePlots || config.scan || config.pseudoExperiments > 0;
    if (persistDiagnostics) {
      const std::string directoryName = "point_" + std::to_string(point->var1) + "_" +
                                        std::to_string(point->var2) + "_" +
                                        std::to_string(point->var3);
      TDirectory* directory = output.file().mkdir(directoryName.c_str());
      if (directory == nullptr) {
        throw CliError("failed to create diagnostic directory for point " + pointString(*point));
      }
      TDirectory::TContext directoryContext(directory);
      fitStatus = calculate();
    } else {
      // Legacy fit code creates transient fixed-name histograms even when
      // plotting is disabled. Keep those objects out of the result file.
      TDirectory::TContext directoryContext;
      gROOT->cd();
      HistogramAttachmentGuard attachment(false);
      fitStatus = calculate();
    }

    if (fitStatus == 1) {
      throw CliError("signal-strength fit failed at point " + pointString(*point));
    }
    if (fitStatus == 2) {
      throw CliError("fit pseudo-experiments failed at point " + pointString(*point));
    }
    tree.Fill();
    std::cout << pointString(*point) << "  signal_strength=" << result.fit_sigScale
              << "  error=" << result.fit_sigScale_Err << '\n';
    if (config.run.verbose) {
      fit.print();
      fit.printEMAT();
    }
  }

  output.file().cd();
  if (tree.Write("", TObject::kOverwrite) <= 0) {
    throw CliError("failed to write the SCAN tree");
  }
  tree.SetDirectory(NULL);
  writeManifest(output.file(), makeManifest("fit", config.input.files, fitSettings(config)));
  output.commit();
  std::cout << "Wrote " << points.size() << " fit result" << (points.size() == 1 ? "" : "s")
            << " to " << config.output << '\n';
  return 0;
}

std::set<std::string> branchSignature(TTree& tree) {
  std::set<std::string> signature;
  TIter next(tree.GetListOfBranches());
  while (TBranch* branch = static_cast<TBranch*>(next())) {
    std::ostringstream entry;
    entry << branch->GetName() << ':';
    if (branch->GetClassName() && branch->GetClassName()[0] != '\0') {
      entry << branch->GetClassName();
    } else {
      TLeaf* leaf = branch->GetLeaf(branch->GetName());
      entry << (leaf ? leaf->GetTypeName() : "unknown");
    }
    signature.insert(entry.str());
  }
  return signature;
}

std::set<Point> resultPoints(TTree& tree, const std::string& path) {
  const char* required[] = {"var1", "var2", "var3"};
  for (std::size_t i = 0; i < sizeof(required) / sizeof(required[0]); ++i) {
    if (!tree.GetBranch(required[i])) {
      throw CliError("SCAN tree in '" + path + "' has no '" + required[i] + "' branch");
    }
  }

  int first = 0;
  int second = 0;
  int third = 0;
  tree.SetBranchAddress("var1", &first);
  tree.SetBranchAddress("var2", &second);
  tree.SetBranchAddress("var3", &third);
  std::set<Point> points;
  const Long64_t entries = tree.GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    if (tree.GetEntry(entry) <= 0) {
      throw CliError("failed to read SCAN entry in '" + path + "'");
    }
    Point point;
    point.var1 = first;
    point.var2 = second;
    point.var3 = third;
    if (!points.insert(point).second) {
      throw CliError("'" + path + "' repeats point " + pointString(point));
    }
  }
  tree.ResetBranchAddresses();
  return points;
}

int executeCombine(const CombineConfig& config) {
  configureRun(config.run);
  std::set<std::string> expectedSignature;
  std::set<Point> seenPoints;
  ResultProvenance expectedProvenance;
  Long64_t totalEntries = 0;

  for (std::size_t i = 0; i < config.inputs.size(); ++i) {
    std::unique_ptr<TFile> input(TFile::Open(config.inputs[i].c_str(), "READ"));
    if (!input || input->IsZombie()) {
      throw CliError("ROOT could not open input '" + config.inputs[i] + "'");
    }
    TTree* tree = dynamic_cast<TTree*>(input->Get("SCAN"));
    if (!tree) {
      throw CliError("'" + config.inputs[i] + "' contains no SCAN tree");
    }
    const ResultProvenance provenance = resultProvenance(*input, config.inputs[i]);
    const std::set<std::string> signature = branchSignature(*tree);
    if (i == 0) {
      expectedSignature = signature;
      expectedProvenance = provenance;
    } else if (signature != expectedSignature) {
      throw CliError("SCAN schema in '" + config.inputs[i] + "' does not match the first input");
    } else if (provenance.command != expectedProvenance.command) {
      throw CliError("result command in '" + config.inputs[i] + "' does not match the first input");
    } else if (provenance.settings != expectedProvenance.settings) {
      throw CliError("analysis settings in '" + config.inputs[i] +
                     "' do not match the first input");
    }
    const std::set<Point> points = resultPoints(*tree, config.inputs[i]);
    if (!config.allowDuplicatePoints) {
      for (std::set<Point>::const_iterator point = points.begin(); point != points.end(); ++point) {
        if (!seenPoints.insert(*point).second) {
          throw CliError("point " + pointString(*point) + " appears in more than one input");
        }
      }
    }
    totalEntries += tree->GetEntries();
  }

  AtomicRootOutput output(config.output, config.run.force);
  output.file().cd();
  TChain chain("SCAN");
  for (std::vector<std::string>::const_iterator input = config.inputs.begin();
       input != config.inputs.end(); ++input) {
    if (chain.Add(input->c_str()) == 0) {
      throw CliError("could not add '" + *input + "' to the result chain");
    }
  }
  TTree* merged = chain.CloneTree(-1, "fast");
  if (!merged || merged->GetEntries() != totalEntries) {
    throw CliError("failed to concatenate all SCAN entries");
  }
  merged->SetName("SCAN");
  merged->SetTitle("Combined Collie results");
  if (merged->Write("", TObject::kOverwrite) <= 0) {
    throw CliError("failed to write the combined SCAN tree");
  }

  std::ostringstream settings;
  settings << "{\"allow_duplicate_points\":" << (config.allowDuplicatePoints ? "true" : "false")
           << ",\"source_command\":\"" << jsonEscape(expectedProvenance.command) << "\""
           << ",\"source_settings\":" << expectedProvenance.settings << '}';
  writeManifest(output.file(), makeManifest("combine", config.inputs, settings.str()));
  output.commit();
  std::cout << "Combined " << totalEntries << " entries from " << config.inputs.size()
            << " files into " << config.output << '\n';
  return 0;
}

void inspectResult(TFile& file, const std::string& path, TTree& tree) {
  std::cout << path << "\n"
            << "  type: Collie result\n"
            << "  entries: " << tree.GetEntries() << "\n"
            << "  branches: " << tree.GetListOfBranches()->GetEntries() << '\n';
  TNamed* manifest = dynamic_cast<TNamed*>(file.Get("collie_run_manifest"));
  if (manifest) {
    std::cout << "  manifest: " << manifest->GetTitle() << '\n';
  }
}

void inspectCollie(const std::string& path, bool details) {
  std::unique_ptr<LoadedInputs> loaded = loadInputs(std::vector<std::string>(1, path));
  LoadedInput& input = loaded->inputs[0];
  CollieChannel* channel = input.loader->accessChannel();
  if (!channel)
    throw CliError("loaded channel is unexpectedly null");

  std::cout << path << "\n"
            << "  type: Collie input\n"
            << "  channel: " << channel->getChannelName() << "\n"
            << "  luminosity: " << channel->getLuminosity() << "\n"
            << "  independent variables: " << channel->getNIndepVariables();
  for (int i = 0; i < channel->getNIndepVariables(); ++i) {
    std::cout << (i == 0 ? " (" : ", ") << channel->getIndepVariableName(i);
  }
  if (channel->getNIndepVariables() > 0)
    std::cout << ')';
  std::cout << "\n  signals: " << channel->getNSignals();
  for (int i = 0; i < channel->getNSignals(); ++i) {
    std::cout << (i == 0 ? " (" : ", ") << channel->getSignalName(i);
  }
  if (channel->getNSignals() > 0)
    std::cout << ')';
  std::cout << "\n  backgrounds: " << channel->getNBackgrounds();
  for (int i = 0; i < channel->getNBackgrounds(); ++i) {
    std::cout << (i == 0 ? " (" : ", ") << channel->getBackgroundName(i);
  }
  if (channel->getNBackgrounds() > 0)
    std::cout << ')';
  std::cout << "\n  models: " << channel->getNModels()
            << "\n  mass points: " << input.points.size();
  if (!input.points.empty()) {
    std::cout << " (" << pointString(input.points.front());
    if (input.points.size() > 1) {
      std::cout << " .. " << pointString(input.points.back());
    }
    std::cout << ')';
  }
  std::cout << '\n';

  if (details) {
    for (std::vector<Point>::const_iterator point = input.points.begin();
         point != input.points.end(); ++point) {
      std::unique_ptr<SigBkgdDist> distribution(
          input.loader->get(point->var1, point->var2, point->var3));
      if (!distribution) {
        throw CliError("failed to inspect point " + pointString(*point));
      }
      std::cout << "    " << pointString(*point) << "  bins=" << distribution->nbins()
                << "  signal=" << distribution->totSignal()
                << "  background=" << distribution->totBkgd()
                << "  data=" << distribution->totData() << '\n';
    }
  }
}

int executeInspect(const InspectConfig& config) {
  for (std::size_t i = 0; i < config.files.size(); ++i) {
    std::unique_ptr<TFile> file(TFile::Open(config.files[i].c_str(), "READ"));
    if (!file || file->IsZombie()) {
      throw CliError("ROOT could not open input '" + config.files[i] + "'");
    }
    TTree* tree = dynamic_cast<TTree*>(file->Get("SCAN"));
    if (tree) {
      inspectResult(*file, config.files[i], *tree);
    } else {
      file->Close();
      file.reset();
      inspectCollie(config.files[i], config.details);
    }
    if (i + 1 != config.files.size())
      std::cout << '\n';
  }
  return 0;
}

struct Diagnostics {
  std::vector<std::string> warnings;
  std::vector<std::string> errors;
};

void validateResult(TFile& file, const std::string& path, Diagnostics& diagnostics) {
  TTree* tree = dynamic_cast<TTree*>(file.Get("SCAN"));
  if (!tree) {
    diagnostics.errors.push_back(path + ": missing SCAN tree");
    return;
  }
  if (tree->GetEntries() <= 0) {
    diagnostics.warnings.push_back(path + ": SCAN tree is empty");
  }
  bool schemaValid = true;
  const char* required[] = {"var1", "var2", "var3"};
  for (std::size_t i = 0; i < sizeof(required) / sizeof(required[0]); ++i) {
    if (!tree->GetBranch(required[i])) {
      diagnostics.errors.push_back(path + ": SCAN has no " + required[i] + " branch");
      schemaValid = false;
    }
  }
  if (schemaValid) {
    try {
      resultPoints(*tree, path);
    } catch (const std::exception& error) {
      diagnostics.errors.push_back(error.what());
    }
  }
}

int executeValidate(const ValidateConfig& config) {
  Diagnostics diagnostics;
  std::vector<std::pair<std::string, std::vector<Point>>> grids;

  for (std::vector<std::string>::const_iterator path = config.files.begin();
       path != config.files.end(); ++path) {
    std::unique_ptr<TFile> file(TFile::Open(path->c_str(), "READ"));
    if (!file || file->IsZombie()) {
      diagnostics.errors.push_back(*path + ": ROOT cannot open file");
      continue;
    }
    if (file->Get("SCAN")) {
      validateResult(*file, *path, diagnostics);
      continue;
    }
    file->Close();
    file.reset();

    try {
      std::unique_ptr<LoadedInputs> loaded = loadInputs(std::vector<std::string>(1, *path));
      grids.push_back(std::make_pair(*path, loaded->points));
      LoadedInput& input = loaded->inputs[0];
      for (std::vector<Point>::const_iterator point = input.points.begin();
           point != input.points.end(); ++point) {
        std::unique_ptr<SigBkgdDist> distribution(
            input.loader->get(point->var1, point->var2, point->var3));
        if (!distribution) {
          diagnostics.errors.push_back(*path + ": cannot load point " + pointString(*point));
          continue;
        }
        const std::vector<std::string> hardProblems = distributionProblems(*distribution, false);
        for (std::vector<std::string>::const_iterator problem = hardProblems.begin();
             problem != hardProblems.end(); ++problem) {
          diagnostics.errors.push_back(*path + ": point " + pointString(*point) + ": " + *problem);
        }
        const std::vector<std::string> allProblems = distributionProblems(*distribution, true);
        for (std::vector<std::string>::const_iterator problem = allProblems.begin();
             problem != allProblems.end(); ++problem) {
          if (std::find(hardProblems.begin(), hardProblems.end(), *problem) == hardProblems.end()) {
            const std::string message = *path + ": point " + pointString(*point) + ": " + *problem;
            if (config.strict) {
              diagnostics.errors.push_back(message);
            } else {
              diagnostics.warnings.push_back(message);
            }
          }
        }
      }
    } catch (const std::exception& error) {
      diagnostics.errors.push_back(error.what());
    }
  }

  if (grids.size() > 1) {
    for (std::size_t i = 1; i < grids.size(); ++i) {
      if (grids[i].second != grids[0].second) {
        diagnostics.errors.push_back(grids[i].first +
                                     ": mass-point grid does not match the first Collie input");
      }
    }
  }

  for (std::vector<std::string>::const_iterator warning = diagnostics.warnings.begin();
       warning != diagnostics.warnings.end(); ++warning) {
    std::cerr << "warning: " << *warning << '\n';
  }
  for (std::vector<std::string>::const_iterator error = diagnostics.errors.begin();
       error != diagnostics.errors.end(); ++error) {
    std::cerr << "error: " << *error << '\n';
  }
  if (!diagnostics.errors.empty()) {
    std::cerr << "Validation failed: " << diagnostics.errors.size() << " error(s), "
              << diagnostics.warnings.size() << " warning(s)\n";
    return 1;
  }
  std::cout << "Validated " << config.files.size() << " file"
            << (config.files.size() == 1 ? "" : "s") << ": " << diagnostics.warnings.size()
            << " warning(s)\n";
  return 0;
}

bool wantsHelp(const std::vector<std::string>& arguments) {
  return std::find(arguments.begin() + std::min<std::size_t>(2, arguments.size()), arguments.end(),
                   "--help") != arguments.end() ||
         std::find(arguments.begin() + std::min<std::size_t>(2, arguments.size()), arguments.end(),
                   "-h") != arguments.end();
}

} // namespace

int run(int argc, char* argv[]) {
  try {
    const std::vector<std::string> arguments = normalizedArguments(argc, argv);
    if (arguments.size() < 2 || arguments[1] == "--help" || arguments[1] == "-h") {
      printGeneralHelp(std::cout);
      return arguments.size() < 2 ? 2 : 0;
    }

    const std::string command = arguments[1];
    if (wantsHelp(arguments)) {
      printCommandHelp(command, std::cout);
      return 0;
    }
    if (command == "inspect")
      return executeInspect(parseInspect(arguments));
    if (command == "validate") {
      return executeValidate(parseValidate(arguments));
    }
    if (command == "limit")
      return executeLimit(parseLimit(arguments));
    if (command == "fit")
      return executeFit(parseFit(arguments));
    if (command == "combine")
      return executeCombine(parseCombine(arguments));
    throw CliError("unknown command '" + command + "'");
  } catch (const CliError& error) {
    std::cerr << "collie: " << error.what() << '\n' << "Try 'collie --help' for usage.\n";
    return 2;
  } catch (const std::exception& error) {
    std::cerr << "collie: " << error.what() << '\n';
    return 1;
  }
}

} // namespace cli
} // namespace collie
