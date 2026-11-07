#include "CollieIOFile.hh"
#include "TRandom.h"
#include "TClass.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <cstdlib>

//#include <limits.h>
//#include <stdlib.h>
//

using namespace std;

namespace Param {
const char* iso_names[150] = {0};

struct G {
  map<int, double> gmap;
  void Read(ifstream& is) {
    iso_names[76] = "Ge76";
    iso_names[82] = "Se82";
    iso_names[100] = "Mo100";
    iso_names[130] = "Te130";
    iso_names[136] = "Xe136";
    int i;
    double g;
    while (true) {
      is >> i >> g;
      if (is.fail())
        break;
      gmap[i] = g;
    }
  }
  double operator[](int i) const {
    map<int, double>::const_iterator g = gmap.find(i);
    if (g != gmap.end())
      return g->second;
    throw range_error("Isotope not defined");
  }
};

struct E {
  string name;
  int isotope;
  string file;
  struct Error {
    string name;
    double err;
    bool unconstrained;
  };
  struct Hist {
    enum { data, signal, background } type;
    vector<int> systs;
  };
  double lifetime;
  map<int, Error> errs;
  map<string, Hist> hists;
  map<int, pair<string, string>> hist_errs_abs;
  map<int, pair<string, string>> hist_errs;
  void Read(ifstream& is, const string& conffile) {
    string t;
    while (true) {
      is >> t;
      if (is.fail())
        break;
      if (t == "N") {
        string n;
        is >> n;
        if (is.fail())
          break;
        name = n;
      } else if (t == "I") {
        int i;
        is >> i;
        if (is.fail())
          break;
        isotope = i;
      } else if (t == "F") {
        string f;
        is >> f;
        if (is.fail())
          break;
        char full_path[4096];
        const string input_path = conffile + "/" + f;
        if (!realpath(input_path.c_str(), full_path))
          throw runtime_error("Cannot resolve input file: " + input_path);
        cout << (conffile + "/" + f) << " --> " << full_path << endl;
        file = full_path;
        //free(full_path);
      } else if (t == "S" || t == "D" || t == "B") {
        string n;
        double l;
        int e;
        vector<int> ee;
        string eee;
        is >> n;
        cout << "Histo name: " << n << endl;
        if (is.fail())
          break;
        if (t == "S")
          is >> l;
        if (is.fail())
          break;
        if (t != "D") {
          getline(is, eee);
          if (is.fail())
            break;
          istringstream eeee(eee);
          while (true) {
            eeee >> e;
            if (eeee.fail())
              break;
            ee.push_back(e);
            cout << "-- added syst " << e << endl;
          }
        }
        hists[n].type = (t == "D" ? Hist::data : t == "S" ? Hist::signal : Hist::background);
        hists[n].systs = ee;
        if (t == "S")
          lifetime = l;
      } else if (t == "E") {
        int i;
        string n;
        double e;
        is >> i >> n >> e;
        if (is.fail())
          break;
        errs[i].name = n;
        errs[i].err = e;
        errs[i].unconstrained = false;
      } else if (t == "U") {
        int i;
        string n;
        double e;
        is >> i >> n >> e;
        if (is.fail())
          break;
        errs[i].name = n;
        errs[i].err = e;
        errs[i].unconstrained = true;
      } else if (t == "HS" || t == "HF") {
        int i;
        string n, h_p, h_m;
        is >> i >> n >> h_p >> h_m;
        if (is.fail())
          break;
        errs[i].name = n;
        errs[i].err = -1.;
        if (t == "HS") {
          hist_errs_abs[i].first = h_p;
          hist_errs_abs[i].second = h_m;
        } else {
          hist_errs[i].first = h_p;
          hist_errs[i].second = h_m;
        }
      }
    }
  }
};

struct M {
  string name;
  map<int, double> mval;
  void Read(ifstream& is) {
    is >> name;
    int i;
    double m;
    while (true) {
      is >> i >> m;
      if (is.fail())
        break;
      mval[i] = m;
    }
  }
  double operator[](int i) const {
    map<int, double>::const_iterator m = mval.find(i);
    if (m != mval.end())
      return m->second;
    throw range_error("Isotope not defined");
  }
};

class Reader {
public:
  G gs;
  vector<M> ms;
  vector<E> es;
  void Read(const string& file) {
    ifstream is(file.c_str());
    string t;
    is >> t;
    if (is.fail())
      throw runtime_error("end of file");
    if (t == "G") {
      gs.Read(is);
    } else if (t == "M") {
      M m;
      m.Read(is);
      ms.push_back(m);
    } else if (t == "E") {
      E e;
      e.Read(is,
             file.substr(0, file.find_last_of('/') != string::npos ? file.find_last_of('/') : 0));
      es.push_back(e);
    }
  }
};
} // namespace Param

namespace {
bool IsTH1D(TH1* h) {
  //return h->Class()->InheritsFrom("TH1D");
  return true;
}
}; // namespace

int run(int argc, char* argv[]) {

  vector<string> paramfs;
  string outfn;
  char c;
  double scale = 1.;
  double target_hl = -1.;
  double target_mass = -1.;
  int NME_TE = 0;
  int NME_MO = 0;
  int NME_XE = 0;
  while ((c = getopt(argc, argv, "o:s:h:m:T:M:X:")) != -1) {
    switch (c) {
      case 'o':
        outfn = optarg;
        break;
      case 'h':
        istringstream(optarg) >> target_hl;
        break;
      case 'm':
        istringstream(optarg) >> target_mass;
        break;
      case 's':
        istringstream(optarg) >> scale;
        break;
      case 'M':
        istringstream(optarg) >> NME_MO;
        break;
      case 'X':
        istringstream(optarg) >> NME_XE;
        break;
      case 'T':
        istringstream(optarg) >> NME_TE;
        break;
      default:
        break;
    }
  }
  for (int i = optind; i < argc; ++i)
    paramfs.push_back(argv[i]);

  Param::Reader rr;

  for (vector<string>::iterator p = paramfs.begin(); p != paramfs.end(); ++p) {
    cout << "reading " << *p << endl;
    rr.Read(*p);
  }

  // load generic NME
  rr.ms.clear();
  for (int nme_ge = 0; nme_ge <= 100; nme_ge += 2) {
    for (int nme_te = NME_TE; nme_te <= NME_TE; ++nme_te) {
      for (int nme_mo = NME_MO; nme_mo <= NME_MO; ++nme_mo) {
        int num = 90000000 + nme_mo * 100000 + nme_te * 1000 + nme_ge;
        Param::M m;
        m.name = Form("%d", num);
        m.mval[136] = 0.1 * NME_XE;
        m.mval[76] = 0.1 * nme_ge;
        m.mval[100] = nme_mo;
        m.mval[130] = nme_te;
        rr.ms.push_back(m);
      }
    }
  }

  const unsigned int n_nme = rr.ms.size();
  cout << "Number of NMEs: " << n_nme << endl;
  //for(int i = 0; i < n_nme; ++i) cout << i << ": " << rr.ms[i].name << endl;

  const unsigned int n_exp = rr.es.size();
  for (int i = 0; i < n_exp; ++i) {
    Param::E& exp = rr.es[i];
    string chan_name = exp.name;
    string filename = "collie_input_";
    filename += chan_name + ".root";

    if (exp.hists.empty())
      continue;

    TFile infile(exp.file.c_str());
    TH1* h = (TH1*)infile.Get(exp.hists.begin()->first.c_str());
    if (!h) {
      cerr << "Cannot find histogram named " << exp.hists.begin()->first << " in " << exp.file
           << endl;
      return 1;
    }
    if (!IsTH1D(h)) {
      cerr << "Error, all input histograms need to be TH1D" << endl;
      return 2;
    }

    CollieIOFile* cfile = new CollieIOFile();
    cfile->initFile(filename.c_str(), chan_name.c_str());
    TFile* ffile = TFile::CurrentFile();

    int Nbins = h->GetNbinsX();
    double Xmin = h->GetBinLowEdge(1);
    double Xmax = h->GetBinLowEdge(Nbins + 1);

    cfile->setInputHist(Xmin, Xmax, Nbins);

    TH1D* h_sig = NULL;
    string sig_name;
    vector<int> sig_sys;
    TH1D* h_data = NULL;
    string data_name;
    vector<TH1D*> h_bgs;
    vector<string> bg_names;
    vector<vector<int>> bg_sys;
    vector<double> valpha;

    for (map<string, Param::E::Hist>::const_iterator hh = exp.hists.begin(); hh != exp.hists.end();
         ++hh) {
      if (hh->second.type == Param::E::Hist::signal) {
        h_sig = (TH1D*)infile.Get(hh->first.c_str());
        if (!h_sig) {
          cerr << "Cannot find histogram " << hh->first << " in file!!!" << endl;
          throw invalid_argument("No histogram");
        }
        if (!IsTH1D(h_sig)) {
          cerr << "Error, all histograms need to be TH1D" << endl;
          return 2;
        }
        sig_name = hh->first;
        sig_sys = hh->second.systs;
      } else if (hh->second.type == Param::E::Hist::data) {
        h_data = (TH1D*)infile.Get(hh->first.c_str());
        if (!h_data) {
          cerr << "Cannot find histogram " << hh->first << " in file!!!" << endl;
          throw invalid_argument("No histogram");
        }
        if (!IsTH1D(h_data)) {
          cerr << "Error, all histograms need to be TH1D" << endl;
          return 2;
        }
        data_name = hh->first;
      } else if (hh->second.type == Param::E::Hist::background) {
        TH1D* h_bg = (TH1D*)infile.Get(hh->first.c_str());
        if (!h_bg) {
          cerr << "Cannot find histogram " << hh->first << " in file!!!" << endl;
          throw invalid_argument("No histogram");
        }
        if (!IsTH1D(h_bg)) {
          cerr << "Error, all histograms need to be TH1D" << endl;
          return 2;
        }
        h_bg->SetTitle(hh->first.c_str());
        h_bgs.push_back(h_bg);
        bg_names.push_back(hh->first);
        bg_sys.push_back(hh->second.systs);
        valpha.push_back(-1.);
      }
    }
    if (!h_sig || !h_data) {
      cerr << "Each experiment must define both signal and data histograms" << endl;
      return 1;
    }
    cfile->createChannel(bg_names);

    for (vector<string>::iterator i = bg_names.begin(); i != bg_names.end(); ++i) {
      cout << "BG name : " << *i << endl;
    }

    int iso = exp.isotope;

    const double neu_mass = target_mass > 0. ? target_mass : 500e-3;
    const double ele_mass = 510.9989e3;
    //mass points
    for (int j = 0; j < n_nme; ++j) {
      Param::M& nme = rr.ms[j];
      const double scale_factor_500meV =
          (rr.gs[iso] * nme[iso] * nme[iso] * neu_mass * neu_mass / ele_mass / ele_mass) *
          exp.lifetime;

      vector<TH1*> cloned_h;

      TH1D* h_sig_scaled = (TH1D*)h_sig->Clone("h_sig_scaled");
      cloned_h.push_back(h_sig_scaled);
      if (nme[iso] >= 0.)
        h_sig_scaled->Scale(scale_factor_500meV);
      else
        h_sig_scaled->Scale(0.);
      if (scale != 1.)
        h_sig_scaled->Scale(scale);
      TH1D* h_data_cloned = (TH1D*)h_data->Clone("h_data_cloned");
      cloned_h.push_back(h_data_cloned);
      if (nme[iso] < 0.)
        h_data_cloned->Scale(0.);
      vector<TH1D*> h_bgs_cloned;
      for (vector<TH1D*>::iterator h = h_bgs.begin(); h != h_bgs.end(); ++h) {
        h_bgs_cloned.push_back((TH1D*)(*h)->Clone((string((*h)->GetName()) + "_cloned").c_str()));
        cloned_h.push_back(h_bgs_cloned.back());
        if (nme[iso] < 0.)
          h_bgs_cloned.back()->Scale(0.);
      }
      cfile->createMassPoint(j + 1, h_data_cloned, h_sig_scaled, -1., h_bgs_cloned, valpha);

      for (vector<int>::iterator e = sig_sys.begin(); e != sig_sys.end(); ++e) {
        if (exp.errs[*e].err >= 0.) {
          cfile->createFlatSigSystematic((exp.errs[*e].name + "_" + nme.name).c_str(),
                                         exp.errs[*e].err, exp.errs[*e].err, j + 1);
          if (exp.errs[*e].err > 0.25)
            cfile->setLogNormalFlag((exp.errs[*e].name + "_" + nme.name).c_str(), true, j + 1);
          if (exp.errs[*e].unconstrained)
            cfile->setSigFloatFlag((exp.errs[*e].name + "_" + nme.name).c_str(), true, j + 1);
        } else {
          const string& n = exp.errs[*e].name;
          map<int, pair<string, string>>::const_iterator err = exp.hist_errs_abs.find(*e);
          if (err != exp.hist_errs_abs.end()) {
            const string& p = err->second.first;
            const string& m = err->second.second;
            TH1D* h_err_pos = (TH1D*)infile.Get(p.c_str())->Clone((p + "_cloned").c_str());
            cloned_h.push_back(h_err_pos);
            h_err_pos->Scale(scale_factor_500meV);
            TH1D* h_err_neg = (TH1D*)infile.Get(m.c_str())->Clone((m + "_cloned").c_str());
            cloned_h.push_back(h_err_neg);
            h_err_neg->Scale(scale_factor_500meV);
            cfile->createShapeSigSystematic((n + "_" + nme.name).c_str(), h_err_pos, h_err_neg,
                                            j + 1);
          } else {
            err = exp.hist_errs.find(*e);
            if (err == exp.hist_errs.end())
              throw range_error("No systematic");
            const string& p = err->second.first;
            const string& m = err->second.second;
            TH1D* h_err_pos = (TH1D*)infile.Get(p.c_str())->Clone((p + "_cloned").c_str());
            cloned_h.push_back(h_err_pos);
            h_err_pos->Scale(scale_factor_500meV);
            TH1D* h_err_neg = (TH1D*)infile.Get(m.c_str())->Clone((m + "_cloned").c_str());
            cloned_h.push_back(h_err_neg);
            h_err_neg->Scale(scale_factor_500meV);
            cfile->createSigSystematic((n + "_" + nme.name).c_str(), h_err_pos, h_err_neg, j + 1);
          }
        }
      }

      int n_bg = 0;
      for (vector<vector<int>>::iterator b = bg_sys.begin(); b != bg_sys.end(); ++b) {
        for (vector<int>::iterator e = b->begin(); e != b->end(); ++e) {
          if (exp.errs[*e].err >= 0.) {
            cfile->createFlatBkgdSystematic(n_bg, (exp.errs[*e].name + "_" + nme.name).c_str(),
                                            exp.errs[*e].err, exp.errs[*e].err, j + 1);
            if (exp.errs[*e].err > 0.25)
              cfile->setLogNormalFlag((exp.errs[*e].name + "_" + nme.name).c_str(), true, j + 1);
            if (exp.errs[*e].unconstrained)
              cfile->setBkgdFloatFlag(n_bg, (exp.errs[*e].name + "_" + nme.name).c_str(), true,
                                      j + 1);
          } else {
            const string& n = exp.errs[*e].name;
            map<int, pair<string, string>>::const_iterator err = exp.hist_errs_abs.find(*e);
            if (err != exp.hist_errs_abs.end()) {
              const string& p = err->second.first;
              const string& m = err->second.second;
              TH1D* h_err_pos = (TH1D*)infile.Get(p.c_str())->Clone((p + "_cloned").c_str());
              cloned_h.push_back(h_err_pos);
              TH1D* h_err_neg = (TH1D*)infile.Get(m.c_str())->Clone((m + "_cloned").c_str());
              cloned_h.push_back(h_err_neg);
              cfile->createShapeBkgdSystematic(n_bg, (n + "_" + nme.name).c_str(), h_err_pos,
                                               h_err_neg, j + 1);
            } else {
              err = exp.hist_errs.find(*e);
              if (err == exp.hist_errs.end())
                throw range_error("No systematic");
              const string& p = err->second.first;
              const string& m = err->second.second;
              TH1D* h_err_pos = (TH1D*)infile.Get(p.c_str())->Clone((p + "_cloned").c_str());
              cloned_h.push_back(h_err_pos);
              TH1D* h_err_neg = (TH1D*)infile.Get(m.c_str())->Clone((m + "_cloned").c_str());
              cloned_h.push_back(h_err_neg);
              cfile->createBkgdSystematic(n_bg, (n + "_" + nme.name).c_str(), h_err_pos, h_err_neg,
                                          j + 1);
            }
          }
        }
        n_bg++;
      }

      for (vector<TH1*>::iterator hhhh = cloned_h.begin(); hhhh != cloned_h.end(); ++hhhh) {
        delete *hhhh;
      }
    }

    // half life only scaling
    TObject* o;
    TH1D* h_sig_scaled = new TH1D;
    o = h_sig_scaled;
    h_sig->Copy(*o);
    h_sig_scaled->SetName("h_sig_scaled");
    TH1D* h_data_cloned = new TH1D;
    o = h_data_cloned;
    h_data->Copy(*o);
    h_data_cloned->SetName("h_data_cloned");
    vector<TH1D*> h_bgs_cloned;
    for (vector<TH1D*>::iterator h = h_bgs.begin(); h != h_bgs.end(); ++h) {
      TH1D* h_bg_cloned = new TH1D;
      o = h_bg_cloned;
      (*h)->Copy(*o);
      h_bg_cloned->SetName((string((*h)->GetName()) + "_cloned").c_str());
      h_bgs_cloned.push_back(h_bg_cloned);
      //h_bgs_cloned.push_back((TH1D*)(*h)->Clone((string((*h)->GetName()) + "_cloned").c_str()));
    }
    if (target_hl > 0.)
      scale = exp.lifetime / target_hl;
    if (scale != 1.)
      h_sig_scaled->Scale(scale);
    cfile->createMassPoint(0, h_data_cloned, h_sig_scaled, -1., h_bgs_cloned, valpha);

    for (vector<int>::iterator e = sig_sys.begin(); e != sig_sys.end(); ++e) {
      if (exp.errs[*e].err >= 0.) {
        cfile->createFlatSigSystematic(exp.errs[*e].name.c_str(), exp.errs[*e].err,
                                       exp.errs[*e].err, 0);
        if (exp.errs[*e].err > 0.25)
          cfile->setLogNormalFlag(exp.errs[*e].name.c_str(), true, 0);
        if (exp.errs[*e].unconstrained)
          cfile->setSigFloatFlag(exp.errs[*e].name.c_str(), true, 0);
      } else {
        const string& n = exp.errs[*e].name;
        map<int, pair<string, string>>::const_iterator err = exp.hist_errs_abs.find(*e);
        if (err != exp.hist_errs_abs.end()) {
          const string& p = err->second.first;
          const string& m = err->second.second;
          TH1D* h_err_pos = (TH1D*)infile.Get(p.c_str())->Clone((p + "_cloned").c_str());
          if (scale != 1.)
            h_err_pos->Scale(scale);
          TH1D* h_err_neg = (TH1D*)infile.Get(m.c_str())->Clone((m + "_cloned").c_str());
          if (scale != 1.)
            h_err_neg->Scale(scale);
          cfile->createShapeSigSystematic(n.c_str(), h_err_pos, h_err_neg, 0);
        } else {
          err = exp.hist_errs.find(*e);
          if (err == exp.hist_errs.end())
            throw range_error("No systematic");
          const string& p = err->second.first;
          const string& m = err->second.second;
          TH1D* h_err_pos = (TH1D*)infile.Get(p.c_str())->Clone((p + "_cloned").c_str());
          if (scale != 1.)
            h_err_pos->Scale(scale);
          TH1D* h_err_neg = (TH1D*)infile.Get(m.c_str())->Clone((m + "_cloned").c_str());
          if (scale != 1.)
            h_err_neg->Scale(scale);
          cfile->createSigSystematic(n.c_str(), h_err_pos, h_err_neg, 0);
        }
      }
    }

    int n_bg = 0;
    for (vector<vector<int>>::iterator b = bg_sys.begin(); b != bg_sys.end(); ++b) {
      for (vector<int>::iterator e = b->begin(); e != b->end(); ++e) {
        if (exp.errs[*e].err >= 0.) {
          cfile->createFlatBkgdSystematic(n_bg, exp.errs[*e].name.c_str(), exp.errs[*e].err,
                                          exp.errs[*e].err, 0);
          if (exp.errs[*e].err > 0.25)
            cfile->setLogNormalFlag(exp.errs[*e].name.c_str(), true, 0);
          if (exp.errs[*e].unconstrained)
            cfile->setBkgdFloatFlag(n_bg, exp.errs[*e].name.c_str(), true, 0);
        } else {
          const string& n = exp.errs[*e].name;
          map<int, pair<string, string>>::const_iterator err = exp.hist_errs_abs.find(*e);
          if (err != exp.hist_errs_abs.end()) {
            const string& p = err->second.first;
            const string& m = err->second.second;
            TH1D* h_err_pos = (TH1D*)infile.Get(p.c_str())->Clone((p + "_cloned").c_str());
            TH1D* h_err_neg = (TH1D*)infile.Get(m.c_str())->Clone((m + "_cloned").c_str());
            cfile->createShapeBkgdSystematic(n_bg, n.c_str(), h_err_pos, h_err_neg, 0);
          } else {
            err = exp.hist_errs.find(*e);
            if (err == exp.hist_errs.end())
              throw range_error("No systematic");
            const string& p = err->second.first;
            const string& m = err->second.second;
            TH1D* h_err_pos = (TH1D*)infile.Get(p.c_str())->Clone((p + "_cloned").c_str());
            TH1D* h_err_neg = (TH1D*)infile.Get(m.c_str())->Clone((m + "_cloned").c_str());
            cfile->createBkgdSystematic(n_bg, n.c_str(), h_err_pos, h_err_neg, 0);
          }
        }
      }
      n_bg++;
    }

    TDirectory* td = ffile->mkdir("params");
    td->cd();
    (new TNamed("HL_nom", Form("%g", target_hl > 0. ? target_hl : (exp.lifetime / scale))))
        ->Write();
    (new TNamed("mass_nom", Form("%g", neu_mass)))->Write();
    for (int j = 0; j < n_nme; ++j)
      (new TNamed(Form("NME_%d", j + 1), rr.ms[j].name.c_str()))->Write();

    cfile->storeFile();
  }

  return 0;
}

int main(int argc, char* argv[]) {
  try {
    return run(argc, argv);
  } catch (const std::exception& error) {
    cerr << "ERROR: " << error.what() << endl;
    return 1;
  }
}
