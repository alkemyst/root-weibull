#include <TCanvas.h>
#include <TGaxis.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMath.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TRandom3.h>

/* HOW TO TEST

   root -l 
   .L weibull.cpp+
   macro();

   OR
   
   root -l 
   .L weibull.cpp+
   double beta=2.;
   double eta=1977.;
   int n=100;
   int maxFailures=5;
   generateFailuresCountTruncated(beta, eta, n, maxFailures);
   macro();

   OR
   
   root -l
   generateFailuresCountTruncated(1, 1000, 300, 2); macro(false);

*/

TRandom3 myDice;
const double inverseBetaPrecision = 1e-4;
const int inverseBetaMaxCycles = 40;

const double alpha_high_cut = 0.999;
const double alpha_low_cut = 1-alpha_high_cut;

const double ranks_high_cut = 0.99;
const double ranks_low_cut = 1e-3;

TCanvas* myCanvas = NULL;
TGraph* rawData = NULL;
TGraph* adjustedData = NULL;

std::vector<double> failure_list;

double inverseBeta(double targetp, double alpha, double beta) {
  double x=targetp;
  double prevx_lower=0;
  double prevx_higher=1;
  double p=TMath::BetaDistI(x, alpha, beta);

  int counter = 0;
  while (fabs(p-targetp)>inverseBetaPrecision) {
    // std::cout << "prevx_lower=" << prevx_lower
    //           << ", prevx_higher=" << prevx_higher
    //           << ", x=" << x
    //           << ", p=" << p << std::endl;
    if (p>targetp) {
      prevx_higher=x;
      x = (x+prevx_lower)/2.;
    } else {
      prevx_lower=x;
      x = (x+prevx_higher)/2.;
    }
    p=TMath::BetaDistI(x, alpha, beta);
    if (counter++>inverseBetaMaxCycles) {
      std::cerr << "WARNING: it took too much to compute inverseBeta(" << targetp
                << ", " << alpha
                << "," << beta <<"). I quit after " << counter << " cycles with p=" << p
                << " instead of p=" << targetp << std::endl;
      return x;
    }
  }
  return x;
}

double ranks(int j_failures, int n_samples) {
  if (j_failures>n_samples) return -1;
  if (n_samples<=0) return -1;
  if (j_failures<=0) return -1; // TODO: add error here
  return inverseBeta(0.5,j_failures,n_samples+1-j_failures);
}

int loadFailures() {
  failure_list.clear();
  failure_list.push_back(137);
  failure_list.push_back(290);
  failure_list.push_back(511);
  failure_list.push_back(803);
  failure_list.push_back(1205);
  failure_list.push_back(-2000);
  failure_list.push_back(-2000);
  failure_list.push_back(-2000);
  failure_list.push_back(-2000);
  failure_list.push_back(-2000);
  return failure_list.size();
}

void generateFailuresCountTruncated(double beta, double eta, int n, int maxFailures) {
  double cycles=0;
  double cyclestep=eta/10000.;
  int iFailures=0;
  failure_list.clear();
  double hazard;
  while ((iFailures<maxFailures)&&(n!=0)) {
    cycles+=cyclestep;
    hazard = beta/eta*pow(cycles/eta, beta-1);
    for (int i=0; i<n; ++i) {
      if (myDice.Rndm()<hazard*cyclestep) {
        iFailures++;
        failure_list.push_back(cycles);
        n--;
      }
    }
  }

  for (int i=0; i<n; ++i) {
    failure_list.push_back(-cycles-cyclestep);
  }
}

double horizontalScale(double cycles) {
  return log10(cycles);
}

double verticalScale(double rank) {
  return log10(log(1/(1-rank)));
}

void computeRanks(double& minCycles, double& maxCycles) {
  std::sort (failure_list.begin(), failure_list.end());
  int n=failure_list.size();
  int j=0;
  double x, y;
  double rank;
  if (rawData) delete rawData;
  if (adjustedData) delete adjustedData;
  rawData = new TGraph();
  adjustedData = new TGraph();
  adjustedData->SetMarkerStyle(8);
  for (auto & it : failure_list) {
    if (it>=0) {
      j++;
      if (j==1) minCycles=it;
      maxCycles=it;
      rank=ranks(j,n);
      std::cout << it << ", j=" << j << ", n=" << n <<", rank=" << rank << std::endl;
      rawData->SetPoint(rawData->GetN(), it, rank);
      x = horizontalScale(it);
      y = verticalScale(rank);
      adjustedData->SetPoint(adjustedData->GetN(), x, y);
    }
  }
}

bool rankRegression(TGraph* aGraph, double& a, double& b) {
  TGraph invertedData;
  double x, y;
  for (int i=0; i<aGraph->GetN(); ++i) {
    aGraph->GetPoint(i, x, y);
    invertedData.SetPoint(i, y, x);
  }
  TFitResultPtr result = invertedData.Fit("pol1", "SQN");
  TFitResult* res_p = result.Get();
  if (!res_p) return false;
  a = res_p->GetParams()[1];
  b = res_p->GetParams()[0];
  b = -b/a; // Move to the inverted plane coordinates
  a = 1/a;  //
  return true;
}

Double_t myHorizontalScale(Double_t *x, Double_t *par) {
  return horizontalScale(x[0]);
}

Double_t myVerticalScale(Double_t *x, Double_t *par) {
  return verticalScale(x[0]);
}

void createCanvasAxes(double coord_x0, double coord_x1, double coord_y0, double coord_y1) {
  myCanvas = new TCanvas("c1","Examples of Gaxis",700,500);
  myCanvas->Draw();
  
  //TF1* nonlinear_x = new TF1("nonlinear_x", "log10(x)", coord_x0, coord_x1);
  TF1* nonlinear_x = new TF1("nonlinear_x", myHorizontalScale, coord_x0, coord_x1, 0);
  TF1* nonlinear_y = new TF1("nonlinear_y", myVerticalScale, coord_y0, coord_y1, 0);

  double margin = 0.1;
  double pad_x0 = nonlinear_x->Eval(coord_x0);
  double pad_y0 = nonlinear_y->Eval(coord_y0);
  double pad_x1 = nonlinear_x->Eval(coord_x1);
  double pad_y1 = nonlinear_y->Eval(coord_y1);
  double pad_dx = pad_x1-pad_x0;
  double pad_dy = pad_y1-pad_y0;

  myCanvas->Range(pad_x0-pad_dx*margin,  pad_y0-pad_dy*margin,  pad_x1+pad_dx*margin,  pad_y1+pad_dy*margin);
  TGaxis* A1 = new TGaxis(pad_x0, pad_y0, pad_x1, pad_y0, "nonlinear_x", 50510, "G");
  A1->SetTitle("Life [cycles]");
  A1->Draw();
  TGaxis* A2 = new TGaxis(pad_x0, pad_y0, pad_x0, pad_y1, "nonlinear_y", 50510, "");
  A2->SetTitle("Occurrence CFD");
  A2->Draw();

}

// line goes y = a*x + b
void drawCL(double alpha, double a, double b,
            double min_cycles_graph=-1, double max_cycles_graph=-1,
            double min_rank_graph=-1, double max_rank_graph=-1) {
  if ((alpha<alpha_low_cut)||(alpha>alpha_high_cut)) return;

  if (min_rank_graph<0) min_rank_graph=ranks_low_cut;
  if (max_rank_graph<0) max_rank_graph=ranks_high_cut;

  min_rank_graph=verticalScale(min_rank_graph);
  max_rank_graph=verticalScale(max_rank_graph);
  
  bool doHorizontalClipLow = (min_cycles_graph>0);
  bool doHorizontalClipHigh = (max_cycles_graph>0);

  min_cycles_graph=horizontalScale(min_cycles_graph);
  max_cycles_graph=horizontalScale(max_cycles_graph);

  double xx, yy, dummy;
  int n=failure_list.size();
  TGraph* bound = new TGraph();

  for (double jj=0.1; jj<n+1; jj+=0.01) {
    yy=inverseBeta(alpha, jj, n+1-jj);
    yy=verticalScale(yy);
    xx = (yy - b)/a;

    if ((doHorizontalClipLow)&&(xx<min_cycles_graph)) continue;
    if ((doHorizontalClipHigh)&&(xx>max_cycles_graph)) break;
    
    yy=inverseBeta(0.5, jj, n+1-jj);
    yy=verticalScale(yy);

    if (yy>max_rank_graph) break;

    if (yy>min_rank_graph)
      bound->SetPoint(bound->GetN(), xx, yy);
  }

  bound->SetLineColor(kBlue);
  bound->Draw("l");
}



/* 

 double estimateProbabilityPoint(double eta, double beta, double cycles, double arank, int j, int n) {
   std::cout << "Point cycles=" << cycles;
   double predictedRank=1-exp(-pow(cycles/eta,beta));
   std::cout << ", predictedRank=" << predictedRank;
   std::cout << ", arank=" << arank;
   double prob=TMath::BetaDistI(predictedRank, j, n+1-j);
   std::cout << "  P="<<prob<<std::endl;

   double cyclesToFail=pow(-log(1-predictedRank), 1/beta)*eta;

   return prob;
 }

*/

// double estimateProbabilityPoints(double eta, double beta) {
//   int j=adjustedData->GetN();
//   int n=failure_list.size();
//   double cycles, arank;
//   std::cout << "Estimating probability for points with eta=" << eta
//             << ", beta=" << beta << ", j=" << j << ", n=" << n << std::endl;
//   for (int i=0; i<rawData->GetN(); ++i) {
//     rawData->GetPoint(i, cycles, arank);
//     estimateProbabilityPoint(eta, beta, cycles, arank, j, n);
//   }
//   return 0;
// }

double estimateConfidence(double alpha_confidence, double beta_characteristic, bool timeTruncated=true) {
  // int j=adjustedData->GetN();
  int r=failure_list.size();
  if (timeTruncated) r++;
  double cycles, rank;
  double totalCycles=0;
  for (int i=0; i<rawData->GetN(); ++i) {
    rawData->GetPoint(i, cycles, rank);
    totalCycles+=pow(cycles, beta_characteristic);
  }
  double chi2= TMath::ChisquareQuantile(alpha_confidence, 2*r);
  //std::cout << "totalCycles=" << totalCycles << std::endl;
  //std::cout << "chi^2(alpha, " << ((timeTruncated) ? "2r+1" : "2r") <<" )=" << chi2 << std::endl;
  std::cout << "Eta > " << pow(2*totalCycles/chi2, 1/beta_characteristic)  << " ("<< alpha_confidence*100<<"% CL)" << std::endl;
  return pow(2*totalCycles/chi2, 1/beta_characteristic);
}

// bool draw = true; double clLimit=0.1;
void macro(bool draw = true, double clLimit=0.1) {
  double minCycles, maxCycles;
  if (!failure_list.size()) loadFailures();
  computeRanks(minCycles, maxCycles);
  // Variable horizontal axis
  double lmc = log10(minCycles);
  double lMc = log10(maxCycles);
  minCycles=pow(10, floor(lmc)-0.00001);
  maxCycles=pow(10, ceil(lMc)+0.00001);
  // Fixed vertical axis for the moment [ forever? :-) ]
  double minRanks = ranks_low_cut;
  double maxRanks = ranks_high_cut;
  double a, b;
  if (rankRegression(adjustedData, a, b)) {
    std::cout << "y = a*x + b" << std::endl;
    std::cout << "beta = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "eta = " << pow(10,-b/a) << std::endl;
    estimateConfidence(0.9, a, false);

    if (draw) {
      createCanvasAxes(minCycles, maxCycles, minRanks, maxRanks);
      double line_min_x = horizontalScale(minCycles);
      double line_max_x = horizontalScale(maxCycles);
      double line_min_y = verticalScale(minRanks);
      double line_max_y = verticalScale(maxRanks);
      { // Let's adjust the line's range not to exceed the graph
        double y_at_minx = a*line_min_x+b;
        double y_at_maxx = a*line_max_x+b;
        if (y_at_minx < line_min_y ) line_min_x = (line_min_y - b)/a;
        if (y_at_maxx > line_max_y ) line_max_x = (line_max_y - b)/a;
      }
      TF1 *myStraight = new TF1("myStraight", "pol1", line_min_x, line_max_x);
      myStraight->SetParameter(0, b);
      myStraight->SetParameter(1, a);

      myStraight->Draw("same");
      adjustedData->Draw("P");  
      drawCL(clLimit,a,b, minCycles, maxCycles, minRanks, maxRanks);
      drawCL(1-clLimit,a,b, minCycles, maxCycles, minRanks, maxRanks);
    }
  } else {
    std::cout << "Sorry, for some kind of weird reason, I could not fit" << std::endl;
  }
}
