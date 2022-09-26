#pragma once

#include<Eigen/Dense>
#include<vector>
#include<algorithm>
#include"InputData.h"

using namespace Eigen;

class HyperParams {
public:
	// Default constructor
	HyperParams() {};

	// Default virtual destructor
	~HyperParams() {};

	// Set sizes
	void setSizes(const unsigned int& nCovariates) {

		_mu0.setZero(nCovariates);
		_Tau0.setZero(nCovariates, nCovariates);
		_workSqrtTau0.setZero(nCovariates, nCovariates);
		_R0.setZero(nCovariates, nCovariates);
		_workInverseR0.setZero(nCovariates, nCovariates);
	}

	// Set defaults
	void setDefaults(const inputData& dataset) {

		unsigned int nSubjects = dataset.nSubjects();
		unsigned int nCov = dataset.nCovariates();

		_UpperNComp = 100;
	
		// For alpha
		_f_df1 = 6.0;
		_f_df2 = 3.0;

		_Comp_ck = 1.0;
		_Comp_ak = 4.0;
		_Comp_bk = 3.0;


		// Values for mu, Tau0, R0, kappa0, kappa1
		// In the following it is useful to have the rows of X as
		// Eigen vectors
		vector<VectorXd> xi(nSubjects);
		for (unsigned int i = 0; i<nSubjects; i++) {
			xi[i].setZero(nCov);
			for (unsigned int j = 0; j<nCov; j++) {
				xi[i](j) = dataset.continuousX(i, j);
			}
		}

		// First compute the hyper prior parameters
		// First compute the hyper parameter for mu_c
		VectorXd mu0 = VectorXd::Zero(nCov);

		MatrixXd Sigma0 = MatrixXd::Zero(nCov, nCov);
		for (unsigned int j = 0; j<nCov; j++) {
			double meanX = 0.0;
			double minX = 0;
			double maxX = 0;
			for (unsigned int i = 0; i<nSubjects; i++) {
				double tmpX = xi[i](j);
				meanX += tmpX;
				if (tmpX<minX || i == 0) {
					minX = tmpX;
				}
				if (tmpX>maxX || i == 0) {
					maxX = tmpX;
				}
			}
			meanX /= (double)(nSubjects);
			double rangeX = maxX - minX;
			mu0(j) = meanX;
			Sigma0(j, j) = rangeX*rangeX;
		}

		_mu0 = mu0;

		MatrixXd Tau0 = Sigma0.inverse();
		_Tau0 = Tau0;
		LLT<MatrixXd> llt;
		_workSqrtTau0 = (llt.compute(Tau0)).matrixU();
		_workLogDetTau0 = log(Tau0.determinant());

	
		//Hyperparameter for mu_c
	
		/*
		if (options.useSparsePrior()) {
			MatrixXd Tau00 = Sigma0.inverse();
			_Tau00 = Tau00;
			LLT<MatrixXd> llt;
			_workSqrtTau00 = (llt.compute(Tau00)).matrixU();
			_workLogDetTau00 = log(Tau00.determinant());
		}

		if (options.useHyperpriorR1()) {
			MatrixXd Tau00 = Sigma0.inverse();
			_Tau00 = Tau00;
			LLT<MatrixXd> llt;
			_workSqrtTau00 = (llt.compute(Tau00)).matrixU();
			_workLogDetTau00 = log(Tau00.determinant());
		} */

		// Now we compute the hyper parameters for Tau_c
		// First we compute the sample covariance
		MatrixXd R = MatrixXd::Zero(nCov, nCov);

		for (unsigned int i = 0; i<nSubjects; i++) {
			R += (xi[i] - mu0)*((xi[i] - mu0).transpose());
		}
		R /= (double)(nSubjects-1);
		//maybe use another coefficient
		R *= (double)(nCov);
		_R0 = R;
		_workInverseR0 = R.inverse();
		_workLogDetR0 = log(R.determinant());


	
		/*if (useSparsePrior) {

			for (unsigned int j = 0; j < nContCovs; j++) {
				R0_sparse(j, j) = 10.0;

			}

			for (unsigned int j = 0; j < (nContCovs - 1); j++) {
				for (unsigned int k = (j + 1); k < nContCovs; k++) {
					R0_sparse(j, k) = 30.0;
					R0_sparse(k, j) = 30.0;
				}
			}
			_R0_sparse = R0_sparse;
		}*/



		//For kappa1	
		_kappa0 = nCov + 2;		

	}

	double shapeAlpha() const {
		return _shapeAlpha;
	}

	void shapeAlpha(const double& s) {
		_shapeAlpha = s;
	}

	double rateAlpha() const {
		return _rateAlpha;
	}

	void rateAlpha(const double& r) {
		_rateAlpha = r;
	}

	/*
	double shapeKappa1() const {
		return _shapeKappa1;
	}

	void shapeKappa1(const double& s) {
		_shapeKappa1 = s;
	}

	double scaleKappa1() const {
		return _scaleKappa1;
	}

	void scaleKappa1(const double& r) {
		_scaleKappa1 = r;
	} */


	/// \brief Return the hyper parameter Tau0Mu0
	const VectorXd& mu0() const {
		return _mu0;
	}

	/// \brief Set the hyper parameter Tau0Mu0
	void mu0(const VectorXd& m0) {
		_mu0 = m0;
	}

	/// \brief Return the hyper parameter Tau0
	const MatrixXd& Tau0() const {
		return _Tau0;

	}

	/// \brief Set the hyper parameter Tau0
	void Tau0(const MatrixXd& T0) {
		_Tau0 = T0;
		_workLogDetTau0 = log(T0.determinant());
		LLT<MatrixXd> llt;
		_workSqrtTau0 = (llt.compute(T0)).matrixU();
	}


	/*
	// \brief Return the hyper parameter Tau00
	const MatrixXd& Tau00() const {
		return _Tau00;

	}

	/// \brief Set the hyper parameter Tau00
	void Tau00(const MatrixXd& T0) {
		_Tau00 = T0;
		_workLogDetTau00 = log(T0.determinant());
		LLT<MatrixXd> llt;
		_workSqrtTau00 = (llt.compute(T0)).matrixU();
	}  */


	/// \brief Return the hyper parameter R0
	const MatrixXd& R0() const {
		return _R0;
	}

	/// \brief Set the hyper parameter R0
	void R0(const MatrixXd& R) {
		_R0 = R;
		_workLogDetR0 = log(R.determinant());
		_workInverseR0 = R.inverse();
	}


	/*
	//used when useSparsePrior=TRUE
	const MatrixXd& R0_sparse() const {
		return _R0_sparse;
	}

	void R0_sparse(const MatrixXd& R) {
		_R0_sparse = R;
	}  */

	/// \brief Return the hyper parameter kappa0
	const double& kappa0() const {
		return _kappa0;
	}

	/// \brief Set the hyper parameter kappa0
	void kappa0(const double& k0) {
		_kappa0 = k0;
	}

	/*
	/// \brief Return the hyper parameter kappa1
	const double& kappa1() const {
		return _kappa1;
	}

	/// \brief Set the hyper parameter kappa1
	void kappa1(const double& k0) {
		_kappa1 = k0;
	} */

	const MatrixXd& workSqrtTau0() const {
		return _workSqrtTau0;
	}

	double workLogDetTau0() const {
		return _workLogDetTau0;
	}

	const MatrixXd& workInverseR0() const {
		return _workInverseR0;
	}

	double workLogDetR0() const {
		return _workLogDetR0;
	}

	vector<double> initAlloc() const {
		return _initAlloc;
	}

	void initAlloc(const vector<double>& c) {
		_initAlloc = c;
	}

	double initAlloc(const unsigned int& j) const {
		return _initAlloc[j];
	}

	unsigned int UpperNComp () const{
		return _UpperNComp;
	}

	void UpperNComp(const unsigned int& NcompBound) {
		_UpperNComp = NcompBound;
	}

	double Comp_ck() const {
		return _Comp_ck;
	}

	double Comp_ak() const {
		return _Comp_ak;
	}

	double Comp_bk() const {
		return _Comp_bk;
	}

	void Comp_ck(const unsigned int& ck) {
		_Comp_ck = ck;
	}

	void Comp_ak(const unsigned int& ak) {
		_Comp_ak = ak;
	}

	void Comp_bk(const unsigned int& bk) {
		_Comp_bk = bk;
	}

	double f_df1() const {
		return _f_df1;
	}

	void f_df1(const double& df1) {
		_f_df1 = df1;
	}

	double f_df2() const {
		return _f_df2;
	}

	void f_df2(const double& df2) {
		_f_df2 = df2;
	}


	// Copy operator
	HyperParams& operator=(const HyperParams& hyperParams) {
		_shapeAlpha = hyperParams.shapeAlpha();
		_rateAlpha = hyperParams.rateAlpha();
		//_shapeKappa1 = hyperParams.shapeKappa1();
		//_scaleKappa1 = hyperParams.scaleKappa1();
		_mu0 = hyperParams.mu0();
		_Tau0 = hyperParams.Tau0();
		//_Tau00 = hyperParams.Tau00();
		
		_R0 = hyperParams.R0();
		//_R00 = hyperParams.R00();
		_kappa0 = hyperParams.kappa0();
		//_kappa00 = hyperParams.kappa00();
		//_kappa1 = hyperParams.kappa1();

		_workSqrtTau0 = hyperParams.workSqrtTau0();
		_workLogDetTau0 = hyperParams.workLogDetTau0();
		//_workSqrtTau00 = hyperParams.workSqrtTau00();
		//_workLogDetTau00 = hyperParams.workLogDetTau00();
		_workInverseR0 = hyperParams.workInverseR0();
		_workLogDetR0 = hyperParams.workLogDetR0();
		_initAlloc = hyperParams.initAlloc();
		//_R0_sparse = hyperParams.R0_sparse();
		_UpperNComp = hyperParams.UpperNComp();
		_Comp_ck = hyperParams.Comp_ck();
		_Comp_ak = hyperParams.Comp_ak();
		_Comp_bk = hyperParams.Comp_bk();
		_f_df1 = hyperParams.f_df1();
		_f_df2= hyperParams.f_df2();
		return *this;
	}


private:
	// Hyper parameters for prior for alpha
	// prior is alpha ~ Gamma(shape,rate)
	double _shapeAlpha;
	double _rateAlpha;

	// Hyper parameters for prior for kappa1
	// prior is kappa1 ~ InvGamma(shapeKappa1,scaleKappa1)
	//double _shapeKappa1;
	//double _scaleKappa1;


	// Hyper parameters for prior for mu (for Normal covariates)
	// Prior is mu ~ N(mu0,inv(Tau0)) or Prior is mu ~ N(mu0,inv(Tau)/nu0) if the Normal inverse Wishart prior is chosen
	VectorXd _mu0;
	MatrixXd _Tau0;
	//MatrixXd _Tau00;

	// When useHyperpriorR1=FALSE, 
	//hyper parameters for prior of Tau (for Normal covariates)
	// Prior is Tau ~ Wishart(R0.inverse(),kappa0) (which has mean R0*kappa0).
	// When useHyperpriorR1=TRUE, 
	//hyper parameters for prior for R1 (for Normal covariates)
	// Prior is Tau ~ Wishart(R0,kappa0)
	// Prior is R1 ~ Wishart(R0,kappa0)
	MatrixXd _R0;
	//MatrixXd _R00;
	double _kappa0;
	//double _kappa00;
	//double _kappa1;



	//Some working variables for speeding up linear algebra
	double _workLogDetTau0;
	MatrixXd _workSqrtTau0;

	//double _workLogDetTau00;
	//MatrixXd _workSqrtTau00;

	double _workLogDetR0;
	MatrixXd _workInverseR0;

	//double _workLogDetR00;
	//MatrixXd _workInverseR00;

	//only used in the sparse prior case
	//MatrixXd _R0_sparse;
	//MatrixXd _R1_init;

	// Initial allocations 
	vector<double> _initAlloc;

	unsigned int _UpperNComp;

	double _Comp_ck;
	double _Comp_ak;
	double _Comp_bk;

	double _f_df1;
	double _f_df2;

};



