
#include <iostream>
#include <algorithm>
#include <string>
#include "TMath.h"
#include "TMatrixDSymEigen.h"
#include "TKNDTree.h"
#include "TError.h"
#include "TKDTree.h"
#include "TStopwatch.h"
#include "Math/Delaunay2D.h"

using namespace std;

ClassImp(TKNDTree)

TKNDTree::TKNDTree(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, Int_t mode, Int_t nbins ) :
fData(data),
fOptions(options),
fWidthFactor(rho),
fNSigma(nSigma),
fRotate(rotate),
fNEvents(nEvents),
fNDim(nDim),
fMode(mode),
fNBins(nbins)

{
    fWgt.resize(data.size(),0);
    
    for (UInt_t i=0; i< data.size(); i++)
    fWgt[i]=1;
    fSigmaRange=5;
    CreatePdf();
}

TKNDTree::TKNDTree(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, std::vector<Double_t>& wgt ,Int_t mode, Int_t nbins) :
fData(data),
fOptions(options),
fWidthFactor(rho),
fNSigma(nSigma),
fRotate(rotate),
fNEvents(nEvents),
fNDim(nDim),
fWgt(wgt),
fMode(mode),
fNBins(nbins)

{
    fSigmaRange=5;
    CreatePdf();
    
    
}


TKNDTree::TKNDTree(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, Int_t mode, Int_t nbins, int sigmarange ) :
fData(data),
fOptions(options),
fWidthFactor(rho),
fNSigma(nSigma),
fRotate(rotate),
fNEvents(nEvents),
fNDim(nDim),
fMode(mode),
fNBins(nbins),
fSigmaRange(sigmarange)
{
    fWgt.resize(data.size(),0);
    
    for (UInt_t i=0; i< data.size(); i++)
    fWgt[i]=1;
    CreatePdf();
}

TKNDTree::TKNDTree(std::vector<Double_t>& data, TString options, Double_t rho, Double_t nSigma, Bool_t rotate,Int_t nEvents, Int_t nDim, std::vector<Double_t>& wgt ,Int_t mode, Int_t nbins, int sigmarange) :
fData(data),
fOptions(options),
fWidthFactor(rho),
fNSigma(nSigma),
fRotate(rotate),
fNEvents(nEvents),
fNDim(nDim),
fWgt(wgt),
fMode(mode),
fNBins(nbins),
fSigmaRange(sigmarange)

{
    CreatePdf();
    
    
}








void TKNDTree::CreatePdf()
{
    // evaluation order of constructor.
    
    
    // set options
    SetOptions();
    // initialization
    Initialize();
    
    
    // copy dataset, calculate sigma_i's, determine min and max event weight
    LoadDataSet();
    
    // determine static and/or adaptive bandwidth
    CalculateBandWidth();
    // Calculates average bandwidth in 5 sigma range
    if (fMode==1)  AverageW();//Calculates average point and h in bin
    if ((fMode==3)&&(fOptions.Contains("b"))) AverageW();
    ComputePDF();  

}

void TKNDTree::SetOptions()
{
    // set the configuration
    
    fOptions.ToLower();
    
    if( fOptions.Contains("a") ) { fWeights = &fWeights1; }
    else  { fWeights = &fWeights0; }
    
}
void TKNDTree::Initialize()
{
    // initialization
   fSqrt2pi  = sqrt(2.0*TMath::Pi()) ;
    
    
    
    fD = static_cast<Double_t>(fNDim);
    
    vector<Double_t> dummy(fNDim,0.);
    //fDataPts.resize(fNEvents,dummy);
    fWeights0.resize(fNEvents,dummy);
    
    
    if (fWidthFactor>0) { fRho.resize(fNDim,fWidthFactor); }
    
    fX0.resize(fNDim,0.);
    fX1.resize(fNDim,0.);
    fX2.resize(fNDim,0.);
    
    fMean.resize(fNDim,0.);
    fSigma.resize(fNDim,0.);
    
    fCovMat = 0;
    fCorrMat= 0;
    fRotMat = 0;
    fSigmaR = 0;
    fDx = new TVectorD(fNDim); fDx->Zero();
    //fDataPtsR.resize(fNEvents,*fDx);
}


double * TKNDTree::GetPoint(Int_t i)
{
   return fData.data() + i*fNDim;

}

void TKNDTree::LoadDataSet()
{
    // copy the dataset and calculate some useful variables
    
    // first some initialization
    
    fNEventsW=0;
    fTData= new Double_t *[fNDim];
    for (int i=0; i<fNDim; i++)
    {
        fTData[i]=new Double_t[fNEvents];
    }
    
    TMatrixD mat(fNDim,fNDim);
    if (!fCovMat) fCovMat = new TMatrixDSym(fNDim);
    if (!fCorrMat) fCorrMat= new TMatrixDSym(fNDim);
    if (!fRotMat) fRotMat = new TMatrixD(fNDim,fNDim);
    if (!fSigmaR) fSigmaR = new TVectorD(fNDim);
    
    mat.Zero();
    fCovMat->Zero();
    fCorrMat->Zero();
    fRotMat->Zero();
    fSigmaR->Zero();
    
    
    for (Int_t j=0; j<fNDim; j++)
    {
        fX0[j]=fX1[j]=fX2[j]=0.;
    }
    
    
    
    for (Int_t i=0; i<fNEvents; i++)
    {
        
        //std::vector<Double_t> point(fData.begin() + i*fNDim, fData.begin() + (i+1)*fNDim);
        double *  point = GetPoint(i);
       // FindPoint(i,point);
		//fDataPts[i]=point;//Creates matrix [nPoints,nDim] with all the data points
        
        for (int j=0;j<fNDim; j++)
        {
            fTData[j][i]=point[j];
        // Creates data matrix needed for TKDTree
            
        }
      // cout<<"OK"<<endl;
        fNEventsW +=fWgt[i] ;
        
        for (Int_t j=0; j<fNDim; j++) {
            for (Int_t k=0; k<fNDim; k++)
            {
                mat(j,k) += point[j] * point[k] * fWgt[i];
            }
            
            fX0[j] += 1. * fWgt[i];
            fX1[j] += point[j] * fWgt[i];
            fX2[j] += point[j] * point[j] * fWgt[i] ;
            if (fX2[j]!=fX2[j]) exit(3);
            
        }
    }
    
    fN = TMath::Power(4./(fNEventsW*(fD+2.)), 1./(fD+4.)) ;
    
    for (Int_t j=0; j<fNDim; j++)
        {

            fMean[j] = fX1[j]/fX0[j];
            fSigma[j] = sqrt(fX2[j]/fX0[j]-fMean[j]*fMean[j]);
            
        }
    
    TMatrixDSym covMatRho(fNDim); // covariance matrix times rho parameters
    for (Int_t j=0; j<fNDim; j++) {
        for (Int_t k=0; k<fNDim; k++) {
            (*fCovMat)(j,k) = mat(j,k)/fX0[j] - fMean[j]*fMean[k];
            covMatRho(j,k) = (mat(j,k)/fX0[j] - fMean[j]*fMean[k]) * fRho[j] * fRho[k];
            
        }
    }
    fCovMat->Print();
    // find decorrelation matrix and eigenvalues (R)
    TMatrixDSymEigen evCalculatorRho(covMatRho);
    *fRotMat = evCalculatorRho.GetEigenVectors();
    *fRotMat = fRotMat->T(); // transpose
    *fSigmaR = evCalculatorRho.GetEigenValues();
    
    
    // set rho = 1 because sigmaR now contains rho
    for (Int_t j=0; j<fNDim; j++) { fRho[j] = 1.; }
    
    for (Int_t j=0; j<fNDim; j++) { (*fSigmaR)[j] = sqrt((*fSigmaR)[j]); }
    
    for (Int_t j=0; j<fNDim; j++) {
        for (Int_t k=0; k<fNDim; k++)
        (*fCorrMat)(j,k) = (*fCovMat)(j,k)/(fSigma[j]*fSigma[k]) ;
    }
    
    
    
    if (!fRotate) {
        fRotMat->Print();
        fSigmaR->Print();
        TMatrixD haar(fNDim,fNDim);
        TMatrixD unit(TMatrixD::kUnit,haar);
        *fRotMat = unit;
        for (Int_t j=0; j<fNDim; j++) { (*fSigmaR)[j] = fSigma[j]; }
        fRotMat->Print();
        fSigmaR->Print();
    }
    
    // use raw sigmas (without rho) for sigmaAvgR
    TMatrixDSymEigen evCalculator(*fCovMat);
    TVectorD sigmaRraw = evCalculator.GetEigenValues();
    for (Int_t j=0; j<fNDim; j++) { sigmaRraw[j] = sqrt(sigmaRraw[j]); }
    
    fSigmaAvgR=1.;
    for (Int_t j=0; j<fNDim; j++) { fSigmaAvgR *= sigmaRraw[j]; }
    fSigmaAvgR = TMath::Power(fSigmaAvgR, 1./fD) ;
    
    
    
    
    fTree=new TKDTreeID(fNEvents,fNDim,fNBins, fTData);
    fTree->Build();
    
    
}

void TKNDTree::AverageW()

{
    
    fNNodes=fTree->GetNNodes();
    fTotalNodes=fTree->GetTotalNodes();
    vector<Double_t> dummy(fNDim,0.);
    fAvrPoints.resize(fTotalNodes-fNNodes,dummy);
    fAvrWeights.resize(fTotalNodes-fNNodes,dummy);
    
    
    std::vector<Double_t> avr(fNDim);
    std::vector<Double_t> avrW(fNDim);
    
    
    
    int count=0;
    for (int i=fNNodes; i<fTotalNodes; i++)
     {
        int* Index=fTree->GetPointsIndexes(i);
        Double_t w=0;
        
        for (int k=0; k<fNDim ; k++)
        
        {
            avr[k]=0;
            avrW[k]=0;
            
        }
        for (int j=0; j <fTree->GetNPointsNode(i) ;j++)
        
        { w+=fWgt[Index[j]];
            double *point=GetPoint(Index[j]);
            for (int k=0; k<fNDim ; k++)
            
            {
                avr[k]+=point[k]*fWgt[Index[j]];
                
                if(fOptions.Contains("a"))
                
                {
                    if (fMode==1)
                    avrW[k]+=fWeights1[Index[j]][k];
                    if (fMode==3) avrW[k]+=fWeights0[Index[j]][k];
                }
                else {
                    avrW[k]+=fWeights0[Index[j]][k];
                     }
                
                
                
            }
            
        }
        for (int k=0; k<fNDim ; k++)
        
        {
            
            avr[k]=(1.0*avr[k])/w;
            avrW[k]=(1.0*avrW[k]/fTree->GetNPointsNode(i));
             //cout<<avr[k]<<endl;
        }
       
        fAvrPoints[count]=avr;
        
        fAvrWeights[count]=avrW;
        count++;
         
    }
    
    
}






void TKNDTree::CalculateBandWidth()
{
    
    cout<< "Entering Calculate bandwidth" <<endl;
    // non-adaptive bandwidth
    // (default, and needed to calculate adaptive bandwidth)
    
    
    for (Int_t i=0; i<fNEvents; i++) {
        vector<Double_t>& weight = fWeights0[i];
        for (Int_t j=0; j<fNDim; j++) {
            weight[j] = fRho[j] * fN * (*fSigmaR)[j] ;
            
        }
    }
    
    
    fWeights = &fWeights0;
    
    // adaptive width
    if(fOptions.Contains("a")) {
        
        double sqrt12=sqrt(12.);
        double sqrtSigmaAvgR=sqrt(fSigmaAvgR);
        cout<< fMode<< endl;
        
        if (fMode!=3){
            
        vector<Double_t> dummy(fNDim,0.);
        fWeights1.resize(fNEvents,dummy);
        
        
        for(Int_t i=0; i<fNEvents; ++i) {
            
            double *point=GetPoint(i);
            
          // vector<Double_t>& x = point;
            Double_t f = TMath::Power( Gauss(point,fWeights0)/fNEventsW , -1./(2.*fD) ) ;
            for (Int_t j=0; j<fNDim; j++) {
                Double_t norm = (fRho[j]*fN*(*fSigmaR)[j]) / sqrtSigmaAvgR ;
                
                fWeights1[i][j] = norm * f / sqrt12 ; // note additional factor of sqrt(12) compared with HEP-EX/0011057
                
                            }
          
        }
             fWeights = &fWeights1;
        }
        else {
            
            
            AverageW();
        
           for(Int_t i=0; i<(int)(fTotalNodes-fNNodes); ++i) {
               
               double *point = GetPoint(i);
               
               vector<Double_t>& x = fAvrPoints[i];
               // vector<Double_t>& x = fAvrPoints[i];
                Double_t f = TMath::Power( GaussAll(x.data())/fNEventsW , -1./(2.*fD) ) ;
               
               for (Int_t j=0; j<fNDim; j++) {
                   
                    Double_t norm = (fRho[j]*fN*(*fSigmaR)[j]) / sqrtSigmaAvgR ;
                    fAvrWeights[i][j]=norm * f / sqrt12 ;
                    
                   
                    
           
            
        }
	}
           
        }
       
    }
    
    
   
}


Double_t TKNDTree::Gauss(double *x, vector<vector<Double_t> >& weights)
{
    // loop over all closest point to x, as determined by loopRange()
    
    if(fNEvents==0) return 0.;
    
    Double_t z=0.;
    
    
    Double_t Sigma=0;
    for (Int_t j=0; j<fNDim; j++) {
        Sigma+=fSigma[j]*fSigma[j];
    }
    
    Sigma=sqrt(Sigma);
    
    
    std::vector<Int_t> res;
    fTree->FindInRange(x,fSigmaRange*Sigma,res);
    
    
    
   
    
    for (UInt_t i=0; i<res.size(); i++) {
        
        Double_t g(1.);
        
        double *point = GetPoint(res[i]);
        
       // const vector<Double_t>& point = fDataPts[res[i]];
        const vector<Double_t>& weight = weights[res[i]];
        
        for (Int_t j=0; j<fNDim; j++) {
            (*fDx)[j] = x[j]-point[j];
        }
        
        if (fNDim>1) {
            *fDx *= *fRotMat; // rotate to decorrelated frame!
        }
        
        for (Int_t j=0; j<fNDim; j++) {
            Double_t r = (*fDx)[j]; //x[j] - point[j];
            Double_t c = 1./(2.*weight[j]*weight[j]);
            
            g *= exp( -c*r*r );
            g *= 1./(fSqrt2pi*weight[j]);
            //cout << "c= "<<c<< " weight[j] " <<weight[j]<<" r "<< r<< " g= " <<g <<endl;
            
        }
        z += (g*fWgt[i]);
    }
    //cout << "z= " << z<< endl;
    return z;
    
    
    
}

Double_t TKNDTree::evaluate(double *x)
{
     Double_t val=0;
    
    if (fMode==1) val = GaussAll(x);
    if (fMode==2) val = Gauss(x,*fWeights);
    if (fMode==3) val = GaussAll(x);
    //cout<<"returning "<<val<<endl;
    
    if (val>=1E-20)
    return val ;
    else
    return (1E-20) ;
}
Double_t TKNDTree::Interpol(double *x)
{
    Double_t val=0;
    val=fInterpolate->Interpolate(x[0],x[1]);
    
    if (val>=1E-20)
        return val ;
    else
        return (1E-20) ;
}


Double_t TKNDTree::GaussAll(double *x)
{
    // loop over all closest point to x, as determined by loopRange()
    
    if(fNEvents==0) return 0.;
    
    Double_t z=0.;
    
    
    int Node=fTree->FindNode(x);
 
    for (int i=0; i< (int)(fTotalNodes-fNNodes); i++)
    {
        Double_t g(1.);
        for (Int_t j=0; j<fNDim; j++) {
            (*fDx)[j] = x[j]-fAvrPoints[i][j];
            //cout << "dX["<< j<< "]= " << (*fDx)[j]<< cout<< endl;
        }
        
        if (fNDim>1) {
            *fDx *= *fRotMat; // rotate to decorrelated frame!
        }
        
        for (Int_t j=0; j<fNDim; j++) {
            Double_t r = (*fDx)[j]; //x[j] - point[j];
            Double_t c = 1./(2.*fAvrWeights[i][j]*fAvrWeights[i][j]);
            //cout<< "g=" << g<<endl;
            g *= exp( -c*r*r);
            g *= 1./(fSqrt2pi*fAvrWeights[i][j]);
          //  cout<<"TEZINE  =   "<<fAvrWeights[i][j]<<endl;
        }
    
        
        //cout<< "Number of points in the node = " << fTree->GetNPointsNode(i+1+fNNodes)<< endl;
        z += (g*1)*fTree->GetNPointsNode(i+fNNodes);
        
        
    }
    //cout << "z= " << z<< endl;
    return z;
    
    
    
}

void TKNDTree::ComputePDF()

{
    double *x,*y;
    Int_t n=0;
if ((fMode==1)||(fMode==3)) 

	{
        n=fTotalNodes-fNNodes;
        
        
        fPDF= new Double_t [n];
        x= new Double_t [n];
        y= new Double_t [n];
        
	for(int i=0;i<fTotalNodes-fNNodes;i++)
		{
            
            x[i]=fAvrPoints[i][0];
            y[i]=fAvrPoints[i][1];
		fPDF[i]=evaluate(fAvrPoints[i].data());
		//cout<< "X= "<<fAvrPoints[i][0]<< "; Y= "<<fAvrPoints[i][1] <<"; Eval= "  <<fPDF[i]<<endl;
		//cout<<"--------------------------------"<<endl;
		}
        
        
        if (fNDim==2)
           fInterpolate=new ROOT::Math::Delaunay2D(n,x,y,fPDF);

}

	else
    {   n=fNEvents;
        
        fPDF= new Double_t [n];
        x= new Double_t [n];
        y= new Double_t [n];
        
       
        for(int i=0;i<fNEvents;i++)
        {
            double *point = GetPoint(i);
            
            x[i]=point[0];
            y[i]=point[1];
            fPDF[i]=evaluate(point);
            //cout<< "X= "<<fDataPts[i][0]<< "; Y= "<<fDataPts[i][1] <<"; Eval= "  <<fPDF[i]<<endl;
            //cout<<"--------------------------------"<<endl;
            if (fNDim==2)
               fInterpolate=new ROOT::Math::Delaunay2D(n,x,y,fPDF);
        }
        
        
    }

}



