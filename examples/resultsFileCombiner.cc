#include <TFile.h>
#include <TTree.h>
#include <CollieLoader.hh>
#include <CLfit.hh>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>

using namespace std;

/*
This program is a utility to combine ROOT output files from Collie that were run on a subset of
parameter space.  For example, for a channel with 10 mass points, a user could run 20 jobs in
parallel (exp., obs., x10 points).  All of these output files can be merged to one global file
using this utility.

*/

void usage() {

  printf("********************************\n");
  printf("\n");
  printf("Usage A: ./resultsFileCombiner.exe newFile.root expList.txt obsList.txt");
  printf("\n");
  printf("Usage B: ./resultsFileCombiner.exe newFile.root expList.txt");
  printf("\n");
  printf(
      "Usage C: ./resultsFileCombiner.exe newFile.root expList.txt obsList.txt exp_m2_List.txt exp_m1_List.txt exp_p1_List.txt exp_p2_List.txt");
  printf("\n");
  printf("********************************\n");
}

vector<TFile*> load_file_list(const string& list_name) {
  ifstream fname_stream(list_name.c_str());
  if (!fname_stream) {
    cerr << "Error, could not open " << list_name << endl;
    abort();
  }

  vector<TFile*> fvec;

  string fname;
  while (fname_stream >> fname) {
    cout << "Reading " << fname << endl;
    TFile* f = new TFile(fname.c_str());

    if (!f->IsOpen()) {
      cerr << "Error, could not open " << fname << endl;
      abort();
    }
    fvec.push_back(f);
  } //end read file list

  return fvec;
} //end load_file_list

void combine(const char* outName, string expList, string obsList, string exp_m2_List,
             string exp_m1_List, string exp_p1_List, string exp_p2_List) {

  bool expObs = (obsList.size() > 0);

  vector<TFile*> expFiles = load_file_list(expList);

  vector<TFile*> obsFiles;
  if (expObs)
    obsFiles = load_file_list(obsList);

  if (expObs) {
    if (expFiles.size() != obsFiles.size()) {
      cout << "Error: Expected and Observed file lists are not the same length!!" << endl;
      return;
    }
  }

  bool expBand = (exp_m2_List.size() > 0);
  vector<TFile*> expFiles_m2, expFiles_m1;
  vector<TFile*> expFiles_p2, expFiles_p1;
  if (expBand) {

    expFiles_m2 = load_file_list(exp_m2_List);
    if (expFiles_m2.size() != expFiles.size()) {
      cerr << "Error: expected and expected_m2 flist lists are not the same length!!" << endl;
      return;
    } //end error

    expFiles_m1 = load_file_list(exp_m1_List);
    if (expFiles_m1.size() != expFiles.size()) {
      cerr << "Error: expected and expected_m1 flist lists are not the same length!!" << endl;
      return;
    } //end error

    expFiles_p1 = load_file_list(exp_p1_List);
    if (expFiles_p1.size() != expFiles.size()) {
      cerr << "Error: expected and expected_p1 flist lists are not the same length!!" << endl;
      return;
    } //end error

    expFiles_p2 = load_file_list(exp_p2_List);
    if (expFiles_p2.size() != expFiles.size()) {
      cerr << "Error: expected and expected_p2 flist lists are not the same length!!" << endl;
      return;
    } //end error

  } //end if expBand

  int m_var1, m_var2, m_var3;
  double m_cls_obs, m_clb_obs, m_clsb_obs, m_clb_sb;
  double m_cls_med, m_clb_med, m_clsb_med;
  double m_clb_sb_p2s, m_clb_sb_p1s, m_clb_sb_m2s, m_clb_sb_m1s;
  double m_cls_med_p2s, m_cls_med_p1s, m_cls_med_m2s, m_cls_med_m1s;

  double m_llrobs, m_llrb, m_llrsb;
  double m_llrb_p2s, m_llrb_p1s, m_llrb_m2s, m_llrb_m1s;
  double m_llrsb_p2s, m_llrsb_p1s, m_llrsb_m2s, m_llrsb_m1s;
  double m_xsec_cl, m_xsec_obs, m_xsec_med, m_signal_scale;
  double m_xsec_med_m2, m_xsec_med_m1;
  double m_xsec_med_p2, m_xsec_med_p1;

  TFile outFile(outName, "RECREATE");
  TTree outTree("SCAN", "SCAN");

  CLpoint clresults;
  clresults.branch(&outTree);

  for (uint i = 0; i < expFiles.size(); i++) {
    TFile* infile = expFiles[i];
    TFile* infile2 = NULL;
    if (expObs)
      infile2 = obsFiles[i];

    TFile* infile_exp_m1(NULL);
    TFile* infile_exp_m2(NULL);
    TFile* infile_exp_p1(NULL);
    TFile* infile_exp_p2(NULL);
    if (expBand) {
      infile_exp_m1 = expFiles_m1[i];
      infile_exp_m2 = expFiles_m2[i];
      infile_exp_p2 = expFiles_p2[i];
      infile_exp_p1 = expFiles_p1[i];
    }

    TTree* inTree = (TTree*)infile->Get("SCAN");
    if (inTree == NULL) {
      printf("Cannot open the 'SCAN' Tree in  '%s'\n", infile->GetName());
      return;
    }

    TTree* inTree2 = NULL;
    if (expObs) {
      inTree2 = (TTree*)infile2->Get("SCAN");
      if (inTree2 == NULL) {
        printf("Cannot open the 'SCAN' Tree in '%s'\n", infile2->GetName());
        return;
      }
    }

    TTree* intree_exp_m1(NULL);
    TTree* intree_exp_m2(NULL);
    TTree* intree_exp_p1(NULL);
    TTree* intree_exp_p2(NULL);
    if (expBand) {

      intree_exp_m1 = (TTree*)infile_exp_m1->Get("SCAN");
      if (intree_exp_m1 == NULL) {
        cerr << "Error, cannot get 'SCAN' Tree in " << infile_exp_m1->GetName() << endl;
      } //end error

      intree_exp_m2 = (TTree*)infile_exp_m2->Get("SCAN");
      if (intree_exp_m2 == NULL) {
        cerr << "Error, cannot get 'SCAN' Tree in " << infile_exp_m2->GetName() << endl;
      } //end error

      intree_exp_p1 = (TTree*)infile_exp_p1->Get("SCAN");
      if (intree_exp_p1 == NULL) {
        cerr << "Error, cannot get 'SCAN' Tree in " << infile_exp_p1->GetName() << endl;
      } //end error

      intree_exp_p2 = (TTree*)infile_exp_p2->Get("SCAN");
      if (intree_exp_p2 == NULL) {
        cerr << "Error, cannot get 'SCAN' Tree in " << infile_exp_p2->GetName() << endl;
      } //end error

    } //end if expBand

    inTree->SetBranchAddress("var1", &m_var1);
    inTree->SetBranchAddress("var2", &m_var2);
    inTree->SetBranchAddress("var3", &m_var3);

    inTree->SetBranchAddress("cls_obs", &m_cls_obs);
    inTree->SetBranchAddress("clb_obs", &m_clb_obs);
    inTree->SetBranchAddress("clsb_obs", &m_clsb_obs);
    inTree->SetBranchAddress("clb_med", &m_clb_med);
    inTree->SetBranchAddress("clsb_med", &m_clsb_med);

    inTree->SetBranchAddress("clb_sb", &m_clb_sb);
    inTree->SetBranchAddress("clb_sb_m2s", &m_clb_sb_m2s);
    inTree->SetBranchAddress("clb_sb_m1s", &m_clb_sb_m1s);
    inTree->SetBranchAddress("clb_sb_p1s", &m_clb_sb_p1s);
    inTree->SetBranchAddress("clb_sb_p2s", &m_clb_sb_p2s);

    inTree->SetBranchAddress("cls_med", &m_cls_med);
    inTree->SetBranchAddress("cls_med_p2s", &m_cls_med_p2s);
    inTree->SetBranchAddress("cls_med_p1s", &m_cls_med_p1s);
    inTree->SetBranchAddress("cls_med_m2s", &m_cls_med_m2s);
    inTree->SetBranchAddress("cls_med_m1s", &m_cls_med_m1s);

    inTree->SetBranchAddress("llrobs", &m_llrobs);
    inTree->SetBranchAddress("llrsb", &m_llrsb);
    inTree->SetBranchAddress("llrsb_p2s", &m_llrsb_p2s);
    inTree->SetBranchAddress("llrsb_p1s", &m_llrsb_p1s);
    inTree->SetBranchAddress("llrsb_m1s", &m_llrsb_m1s);
    inTree->SetBranchAddress("llrsb_m2s", &m_llrsb_m2s);

    inTree->SetBranchAddress("llrb", &m_llrb);
    inTree->SetBranchAddress("llrb_p2s", &m_llrb_p2s);
    inTree->SetBranchAddress("llrb_p1s", &m_llrb_p1s);
    inTree->SetBranchAddress("llrb_m1s", &m_llrb_m1s);
    inTree->SetBranchAddress("llrb_m2s", &m_llrb_m2s);

    inTree->SetBranchAddress("xsec_cl", &m_xsec_cl);

    if (expObs)
      inTree2->SetBranchAddress("xsec_obsfactor", &m_xsec_obs);
    else
      inTree->SetBranchAddress("xsec_obsfactor", &m_xsec_obs);

    if (expBand) {
      //if variations are stored in separate trees they are put in the xsec_medfactor variable,
      //  not the xsec_medfactor_[pm][12]s variables
      intree_exp_m1->SetBranchAddress("xsec_medfactor", &m_xsec_med_m1);
      intree_exp_m2->SetBranchAddress("xsec_medfactor", &m_xsec_med_m2);
      intree_exp_p1->SetBranchAddress("xsec_medfactor", &m_xsec_med_p1);
      intree_exp_p2->SetBranchAddress("xsec_medfactor", &m_xsec_med_p2);
    } //end if separate trees for +- 1/2 sigma limits

    else { //otherwise, always get from same place as expected limit Tree
      inTree->SetBranchAddress("xsec_medfactor_m1s", &m_xsec_med_m1);
      inTree->SetBranchAddress("xsec_medfactor_m2s", &m_xsec_med_m2);
      inTree->SetBranchAddress("xsec_medfactor_m1s", &m_xsec_med_p1);
      inTree->SetBranchAddress("xsec_medfactor_m2s", &m_xsec_med_p2);
    }

    inTree->SetBranchAddress("xsec_medfactor", &m_xsec_med);
    inTree->SetBranchAddress("signal_scale", &m_signal_scale);

    int len = inTree->GetEntries();
    for (int i = 0; i < len; i++) {
      inTree->GetEntry(i);
      if (expObs)
        inTree2->GetEntry(i);
      if (expBand) {
        intree_exp_m1->GetEntry(i);
        intree_exp_m2->GetEntry(i);
        intree_exp_p1->GetEntry(i);
        intree_exp_p2->GetEntry(i);
      }

      clresults.reset(m_var1, m_var2, m_var3);

      clresults.cls_obs = m_cls_obs;
      clresults.clb_obs = m_clb_obs;
      clresults.clsb_obs = m_clsb_obs;

      clresults.cls_med = m_cls_med;
      clresults.clb_med = m_clb_med;
      clresults.clsb_med = m_clsb_med;

      clresults.cls_med_p2s = m_cls_med_p2s;
      clresults.cls_med_p1s = m_cls_med_p1s;
      clresults.cls_med_m2s = m_cls_med_m2s;
      clresults.cls_med_m1s = m_cls_med_m1s;

      clresults.llrobs = m_llrobs;
      clresults.llrsb = m_llrsb;
      clresults.llrsb_p2s = m_llrsb_p2s;
      clresults.llrsb_p1s = m_llrsb_p1s;
      clresults.llrsb_m1s = m_llrsb_m1s;
      clresults.llrsb_m2s = m_llrsb_m2s;
      clresults.llrb = m_llrb;
      clresults.llrb_p2s = m_llrb_p2s;
      clresults.llrb_p1s = m_llrb_p1s;
      clresults.llrb_m1s = m_llrb_m1s;
      clresults.llrb_m2s = m_llrb_m2s;

      clresults.xsec_cl = m_xsec_cl;
      clresults.xsec_obsfactor = m_xsec_obs;
      clresults.xsec_medfactor = m_xsec_med;

      clresults.xsec_medfactor_m1s = m_xsec_med_m1;
      clresults.xsec_medfactor_m2s = m_xsec_med_m2;
      clresults.xsec_medfactor_p1s = m_xsec_med_p1;
      clresults.xsec_medfactor_p2s = m_xsec_med_p2;

      outTree.Fill();
    }
    infile->Close();
    if (expObs)
      infile2->Close();
  }
  outFile.Write();
}

int main(int argc, char* argv[]) {
  timeval a, b;
  gettimeofday(&a, NULL);

  if (argc == 8)
    combine(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
  else if (argc == 4)
    combine(argv[1], argv[2], argv[3], string(""), string(""), string(""), string(""));
  else if (argc == 3)
    combine(argv[1], argv[2], "", "", "", "", "");
  else
    usage();

  gettimeofday(&b, NULL);
  double deltat = a.tv_sec * 1000.0 + a.tv_usec / 1000.0;
  deltat = b.tv_sec * 1000.0 + b.tv_usec / 1000.0 - deltat;
  printf(" %f sec run time\n", deltat / 1000);
  printf("\n");
  return 0;
}
