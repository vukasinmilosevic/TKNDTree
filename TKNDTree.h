/// @(#)root/mathcore:$Id: TKNDTree.h,v 1.0
/// Author: Vukasin Milosevic, Lorenzo Moneta

/*************************************************************************
* Copyright (C) 2015 ROOT Math Team                                     *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/



/*  Kernel estimation is one of the non-parametric methods used for estimation of probability density function. Its first ROOT implementation, as part of RooFit package, has one major issue, its evaluation time is extremely slow making in almost unusable. The goal of this class TKNDTree, which follows the original idea of kernel estimation, is to greatly improve the evaluation time (using the TKTree class for storing the data and creating different user-controlled modes of evaluation) and add the interpolation option, for 2D case, with the help of the new Delaunnay2D class.   */

/// Header file for the class TKNDTree

#ifndef TKNDTree_H
#define TKNDTree_H
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include <map>
#include <vector>
#include <string>
#include "TKDTree.h"
#include "TStopwatch.h" 
#include "Math/Delaunay2D.h"

class TKNDTree
{
     ///Input data array (x1,y1,z1,...)
	std::vector<Double_t> fData;
    ///Corresponding weights (w1,w2,w3,...)
	std::vector<Double_t> fWgt;
    TString fOptions;
    ///Number of dimensions
	Int_t fNDim;
    ///Number of events (data points)
	Int_t fNEvents;
	Double_t fSqrt2pi;
	Double_t fD;
    Double_t fN;
    ///Sum of all the data point weights
	Double_t fNEventsW;
    /// User inputed mode of calculation of bandwiht and evaluation of the pdf
	Int_t fMode;
    
 /*fMode == 1: Gauss is used to calculate the adptive bandwidth and GaussAll is used to evaluate the pdf
   fMode == 2: Gauss is used for both the calculation of the adaptive bandwith and for the evaluation
   fMode == 3: GaussAll is used for both the calculation of the adaptive bandwith and for the evaluation*/
    
    /// Maximum number of points inside each bin
	Int_t fNBins;

	// std::vector<Double_t> fX;
 	std::vector<Double_t> fX0, fX1, fX2;
	std::vector<Double_t> fMean, fSigma;
	Double_t fWidthFactor;
    /// User inputed sigma range
	Double_t fNSigma;
	Bool_t fRotate;
	 std::vector<Double_t> fRho;

	 TMatrixDSym* fCovMat; ///Covariance matrix
	 TMatrixDSym* fCorrMat; /// Correlation matrix
  	 TMatrixD* fRotMat; /// Rotation matrix
 	 TVectorD* fSigmaR; 
  	 TVectorD* fDx;
 	 Double_t fSigmaAvgR;

	 std::vector<std::vector<Double_t> > fDataPts;
	 std::vector<TVectorD> fDataPtsR;
	 std::vector<std::vector<Double_t> > fWeights0;/// Fixed bandwidth
	 std::vector<std::vector<Double_t> > fWeights1;/// Adaptive bandwidth
	 std::vector<std::vector<Double_t> >* fWeights;
     TKDTreeID *fTree;/// TKDTree class member used to store the data
	 std::vector<Double_t> fAvrP;
	 std::vector<Double_t> fAvrW;
	 Double_t **fTData; /// Data structure needed for the creation of TKDTree
	 std::vector<std::vector<Double_t> > fAvrPoints;/// Average data point per bin
	 std::vector<std::vector<Double_t> > fAvrWeights;// Average bandwith per bin
	 std::vector<Double_t> fAvrWgt;// Average weight of data per bin
	 int fTotalNodes;
	 int fNNodes;
	 int fSigmaRange;
		double *fPDF;
     ROOT::Math::Delaunay2D *fInterpolate;
	
 	
public:
  	TKNDTree(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, Int_t mode, Int_t nbins );
	TKNDTree(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, std::vector<Double_t>& wgt, Int_t mode, Int_t nbins );
  	TKNDTree(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, Int_t mode, Int_t nbins, int sigmarange);
	TKNDTree(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, std::vector<Double_t>& wgt, Int_t mode, Int_t nbins, int sigmarange );
	Double_t evaluate(double *x);
    Double_t Interpol(double *x);
protected:

    void     CreatePdf() ;
    void     SetOptions();
	void     Initialize();
	void     LoadDataSet();/// Loads the data points, creates the TKDTree and calculates all the variables needed for further calculation (mean and sigma values,rotation matrix,...)
	void     CalculateBandWidth();
	void     AverageW(); /// Calculates the average bin point, bandwith and corresponding weight per bin
    Double_t Gauss(double *x, std::vector<std::vector<Double_t> >& weights) ;/// Evaluation function that uses TKDTree::FindInRange method in calculation
	Double_t GaussAll(double *x);/// Evaluation function that uses the average bin point for calculation
    void ComputePDF() ;/// Precomputes the PDF for each data point
    double * GetPoint(Int_t i);
};

#endif
