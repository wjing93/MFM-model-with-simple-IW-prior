#pragma once

#include<Eigen/Dense>
#include<vector>
#include<algorithm>

#include"HyperParams.h"
#include"pdfDis.h"


using namespace Eigen;


class Params {

public:
	/// \brief Default constructor
	Params() {};

	/// \brief Destructor
	~Params() {};

	/// \brief Function to set the sizes of the various class members
	void setSizes(const unsigned int& nSubjects,
		const unsigned int& nCovariates,
		const unsigned int& nCompInit) {

		// Initially make the maximum number of clusters  the bigger or
		// the initial number of clusters and 150.
		// This is only used for initial memory allocation
		// And will ensure that at this initial time sufficient
		// space is allocated to prevent future allocations
		unsigned int NComp = nCompInit;

		_NComp = NComp;

		// Resize all the objects and set to 0
		_logPsi.resize(NComp);
		_Psi.resize(NComp);
		for (unsigned int c = 0; c<NComp; c++) {
			_logPsi[c] = 0.0;
			_Psi[c] = 0.0;
		}
		_mu.resize(NComp);
		_workNXInCluster.resize(NComp, 0);
		_Tau.resize(NComp);
		_workSqrtTau.resize(NComp);
		_workLogDetTau.resize(NComp);
		_Sigma.resize(NComp);


		for (unsigned int c = 0; c<NComp; c++) {

			_mu[c].setZero(nCovariates);

			_Tau[c].setZero(nCovariates, nCovariates);
			_workSqrtTau[c].setZero(nCovariates, nCovariates);
			_Sigma[c].setZero(nCovariates, nCovariates);

		}


		_z.resize(nSubjects);
		_workContinuousX.resize(nSubjects);
		_workLogPXiGivenZi.resize(nSubjects);

	}

	/// \brief Return the number of clusters
	unsigned int NComp() const {
		return _NComp;
	}

	/// \brief Set the number of components
	void NComp(const unsigned int& nComp) {

		_NComp = nComp;

		// Check if we need to do a resize of the
		// number of various vectors
		unsigned int prevNComp = _logPsi.size();
		unsigned int nCov = nCovariates();

		_mu.resize(nComp);

		_logPsi.resize(nComp);
		_Psi.resize(nComp);
		_workNXInCluster.resize(nComp);

		_Tau.resize(nComp);
		_workSqrtTau.resize(nComp);
		_workLogDetTau.resize(nComp);
		_Sigma.resize(nComp);

		if (nComp>prevNComp) {
			//unsigned int nCov = nCovariates();
			for (unsigned int c = prevNComp; c < nComp; c++) {
				_workNXInCluster[c] = 0;
				_mu[c].setZero(nCov);

				_Tau[c].setZero(nCov, nCov);
				_workSqrtTau[c].setZero(nCov, nCov);
				_Sigma[c].setZero(nCov, nCov);

			}
		}
	}


	/// \brief Return the probabilities of allocation
	vector<double> logPsi() const {
		return _logPsi;
	}

	/// \brief Return the probability of allocation to cluster c
	double logPsi(const unsigned int& c) const {
		return _logPsi[c];
	}

	vector<double> Psi() const {
		return _Psi;
	}

	/// \brie Return the probability of allocation to cluster c
	double Psi(const unsigned int& c) const {
		return _Psi[c];
	}


	/// \brief Set the probability of allocation to cluster c
	void Psi(const unsigned int& c, const double& PsiVal) {
		_Psi[c] = PsiVal;
		_logPsi[c] = log(PsiVal);
	}

	/// \brief Set the vector of probability of allocations to clusters
	void Psi(const vector<double>& PsiVec) {
		_Psi = PsiVec;
		for (unsigned int c=0; c < NComp(); c++) {
			_logPsi[c] = log(PsiVec[c]);
		}
		
	}


	/// \brief Return the vector of normal means mu
	const vector<VectorXd>& mu() const {
		return _mu;
	}

	/// \brief Return the normal mean mu for cluster c
	const VectorXd& mu(const unsigned int& c) const {
		return _mu[c];
	}

	/// \brief Return the normal mean mu for cluster c covariate j
	double mu(const unsigned int& c, const unsigned int& j) const {
		return _mu[c](j);
	}

	unsigned int nSubjects() {
		return _z.size();
	}

	unsigned int nCovariates() {
		return _mu[0].size();
	}


	/// \brief Set the normal mean for cluster c
	void mu(const unsigned int& c, const VectorXd& muVec, bool update) {

		_mu[c] = muVec;

		if (update) {
			unsigned int nSbj = nSubjects();
			unsigned int nCov = nCovariates();

			VectorXd xi = VectorXd::Zero(nCov);

			for (unsigned int i = 0; i<nSbj; i++) {
				if (z(i) == (int)c) {
					for (unsigned int j = 0; j<nCov; j++) {
						xi(j) = workContinuousX(i, j);
					}

					_workLogPXiGivenZi[i] = logPdfMultivarNormal(nCov, xi, muVec, workSqrtTau(c), workLogDetTau(c));


				}
			}

		}
		
	}

	

	/// \brief Return the vector of precision matrices Tau
	const vector<MatrixXd>& Tau() const {
		return _Tau;
	}

	/// \brief Return the precision matrix Tau for cluster c
	const MatrixXd& Tau(const unsigned int& c) const {
		return _Tau[c];
	}

	/// \brief Set the precision matrix Tau for cluster c
	void Tau(const unsigned int& c, const MatrixXd& TauMat) {

		_Tau[c] = TauMat;
		Sigma(c, TauMat.inverse());
		workLogDetTau(c, log(TauMat.determinant()));
		LLT<MatrixXd> llt;
		MatrixXd sqrtTau = (llt.compute(TauMat)).matrixU();
		workSqrtTau(c, sqrtTau);

		unsigned int nSbj = nSubjects();
		unsigned int nCov = nCovariates();

		for (unsigned int i = 0; i<nSbj; i++) {
			VectorXd xi = VectorXd::Zero(nCov);
			if (z(i) == (int)c) {
				for (unsigned int j = 0; j<nCov; j++) {
					xi(j) = workContinuousX(i, j);
				}
				_workLogPXiGivenZi[i] = logPdfMultivarNormal(nCov, xi, mu(c), workSqrtTau(c), workLogDetTau(c));
			}
		}
	}

	/// \brief Return the covariance matrix Sigma element j1,j2 for cluster c
	double Tau(const unsigned int& c, const unsigned int& j1, const unsigned int& j2) const {
		return _Tau[c](j1, j2);
	}


	/// \brief Return the scale matrix R1 for the Wishart distribution of Tau_c 
	/*
	const MatrixXd& R1() const {
		return _R1;
	}

	/// \brief Set the scale matrix R1 for the Wishart distribution of Tau_c
	void R1(const MatrixXd& R1Mat) {
		_R1 = R1Mat;
		workLogDetR1(log(R1Mat.determinant()));
		workInverseR1(R1Mat.inverse());
	}

	/// \brief Return the element j1,j2 for R1
	double R1(const unsigned int& j1, const unsigned int& j2) const {
		return _R1(j1, j2);
	}

	const VectorXd& mu00() const {
		return _mu00;
	}

	void mu00(const VectorXd& muVec) {
		_mu00 = muVec;
	}

	double mu00(const unsigned int& j) const {
		return _mu00(j);
	}
	*/


	/// \brief Return the vector of covariance matrices Sigma
	const vector<MatrixXd>& Sigma() const {
		return _Sigma;
	}

	/// \brief Return the covariance matrix Sigma for cluster c
	const MatrixXd& Sigma(const unsigned int& c) const {
		return _Sigma[c];
	}

	/// \brief Set the covariance matrix Sigma for cluster c
	void Sigma(const unsigned int& c, const MatrixXd& SigmaMat) {
		_Sigma[c] = SigmaMat;
	}

	/// \brief Return the covariance matrix Sigma element j1,j2 for cluster c
	double Sigma(const unsigned int& c, const unsigned int& j1, const unsigned int& j2) const {
		return _Sigma[c](j1, j2);
	}


	/// \brief Return the hyper parameter alpha
	double alpha() const {
		return _alpha;
	}

	/// \brief Set the hyper parameter alpha
	void alpha(const double& alphaVal) {
		_alpha = alphaVal;
	}

	/// \brief Return the hyper parameter epsilon

	/// \brief Return the allocation variables
	const vector<int>& z() const {
		return _z;
	}

	/// \brief Return the allocation variable of the ith subject
	int z(const unsigned int& i) const {
		return _z[i];
	}

	/// \brief Set the ith allocation variable to cluster c
	void z(const unsigned int& i, const int& c, bool update) {

		if (update) {
			unsigned int nCov = nCovariates();

			VectorXd xi = VectorXd::Zero(nCov);
			for (unsigned int j = 0; j < nCov; j++) {
				xi(j) = workContinuousX(i, j);
			}

			_workLogPXiGivenZi[i] = logPdfMultivarNormal(nCov, xi, mu(c), workSqrtTau(c), workLogDetTau(c));
		}

		_z[i] = c;

	}


	/*
	const HyperParams& hyperParams() const {
		return _hyperParams;
	}

	HyperParams& hyperParams() {
		return _hyperParams;
	}

	void hyperParams(const HyperParams& hyPar) {
		_hyperParams = hyPar;
	} */

	const vector<unsigned int>& workNXInCluster() const {
		return _workNXInCluster;
	}

	const unsigned int& workNXInCluster(unsigned int c) const {
		return _workNXInCluster[c];
	}

	void workNXInCluster(const vector<unsigned int>& nXInCluster) {
		for (unsigned int c = 0; c<_workNXInCluster.size(); c++) {
			if (c<nXInCluster.size()) {
				_workNXInCluster[c] = nXInCluster[c];
			}
			else {
				_workNXInCluster[c] = 0;
			}
		}
	}

	void workNXInCluster(const unsigned int& c, const unsigned int& n) {
		_workNXInCluster[c] = n;
	}


	const vector<vector<double> >& workContinuousX() const {
		return _workContinuousX;
	}

	double workContinuousX(const unsigned int& i, const unsigned int& j) const {
		return _workContinuousX[i][j];
	}


	void workContinuousX(const vector<vector<double> >& X) {
		_workContinuousX = X;
	}

	void workContinuousX(const unsigned int& i, const unsigned int& j, const double& x) {
		_workContinuousX[i][j] = x;
	}

	/// \brief Return the conditional probabilities
	const vector<double>& workLogPXiGivenZi() const {
		return _workLogPXiGivenZi;
	}

	double workLogPXiGivenZi(const unsigned int& i) const {
		return _workLogPXiGivenZi[i];
	}

	void workLogPXiGivenZi(const unsigned int& i, const double& newVal) {
		_workLogPXiGivenZi[i] = newVal;
	}


	/*
	unsigned int workNClusInit() const {
		return _workNClusInit;
	}

	void workNClusInit(const unsigned int& nClusInit) {
		_workNClusInit = nClusInit;
	} */

	const vector<MatrixXd>& workSqrtTau() const {
		return _workSqrtTau;
	}


	const MatrixXd& workSqrtTau(const unsigned int& c) const {
		return _workSqrtTau[c];
	}

	void workSqrtTau(const unsigned int& c, const MatrixXd& sqrtTau) {
		_workSqrtTau[c] = sqrtTau;
	}

	const vector<double>& workLogDetTau() const {
		return _workLogDetTau;
	}

	double workLogDetTau(const unsigned int& c) const {
		return _workLogDetTau[c];
	}

	void workLogDetTau(const unsigned int& c, const double& logDetTau) {
		_workLogDetTau[c] = logDetTau;
	}

	/*
	const double& workLogDetR1() const {
		return _workLogDetR1;
	}

	void workLogDetR1(const double& logDetR1) {
		_workLogDetR1 = logDetR1;
	}

	const MatrixXd& workInverseR1() const {
		return _workInverseR1;
	}

	void workInverseR1(const MatrixXd& R1) {
		_workInverseR1 = R1;
	} */

	//swap all the cluster specific parameters
	void switchLabels(int c1, int c2) {
		//Covariate parameters including working parameters

		VectorXd muTmp = _mu[c1];
		_mu[c1] = _mu[c2];
		_mu[c2] = muTmp;

		MatrixXd SigmaTmp = _Sigma[c1];
		_Sigma[c1] = _Sigma[c2];
		_Sigma[c2] = SigmaTmp;

		MatrixXd TauTmp = _Tau[c1];
		_Tau[c1] = _Tau[c2];
		_Tau[c2] = TauTmp;
		MatrixXd sqrtTauTmp = _workSqrtTau[c1];
		_workSqrtTau[c1] = _workSqrtTau[c2];
		_workSqrtTau[c2] = sqrtTauTmp;
		double logDetTauTmp = _workLogDetTau[c1];
		_workLogDetTau[c1] = _workLogDetTau[c2];
		_workLogDetTau[c2] = logDetTauTmp;


		//component weight switch
		double PsiTmp = _Psi[c2];
		_Psi[c2] = _Psi[c1];
		_Psi[c1] = PsiTmp;

		double logPsiTmp = _logPsi[c2];
		_logPsi[c2] = _logPsi[c1];
		_logPsi[c1] = logPsiTmp;


		for (unsigned int i = 0; i<nSubjects(); i++) {
			if (_z[i] == (int)c1) {
				_z[i] = c2;
			}
			else if (_z[i] == (int)c2) {
				_z[i] = c1;
			}
		}
		double workNXInClusterTmp = _workNXInCluster[c1];
		_workNXInCluster[c1] = _workNXInCluster[c2];
		_workNXInCluster[c2] = workNXInClusterTmp;	

	}

	unsigned int NClust() const {
		return _NClust;
	}

	void NClust(unsigned int k) {
		_NClust = k;
	}

	void RelabelClust() {
		int Ncomp = NComp();
		int Nempty = 0;
		int Nfill = 0;
		for (int c = 0; c < Ncomp; c++) {
			if (workNXInCluster(c) == 0) {
				Nempty++;
			}
			else {
				Nfill++;
			}
		}

		//
		NClust(Nfill);


		if (Nempty > 0) {
			std::vector<int> FillInd(Nfill, 0);
			std::vector<int> EmptyInd(Nempty, 0);

			int fillcount = 0;
			int emptycount = 0;
			for (int c = 0; c < Ncomp; c++) {
				if (workNXInCluster(c) == 0) {
					EmptyInd[emptycount] = c;
					emptycount++;
				}
				else {
					FillInd[fillcount] = c;
					fillcount++;
				}
			}

			for (int i = 0; i < Nfill; i++) {
				int currFill = FillInd[i];
				for (int j = 0; j < Nempty; j++) {
					int currEmpty = EmptyInd[j];
					if (currEmpty < currFill) {
						switchLabels(currEmpty, currFill);
						FillInd[i] = currEmpty;
						EmptyInd[j] = currFill;
						std::sort(EmptyInd.begin(), EmptyInd.end());
						break;
					}
				}
			}

		}
		//Obtain the indices for filled and unfilled clusters
	
	}
	

		
	/// \brief Copy operator
	Params& operator=(const Params& params) {
		//number of components
		_NComp = params.NComp();
		//number of clusters
		_NClust = params.NClust();
		_logPsi = params.logPsi();
		_mu = params.mu();
		_Tau = params.Tau();
		//_R1 = params.R1();
		//_mu00 = params.mu00();
		_Sigma = params.Sigma();
		_alpha = params.alpha();
		_z = params.z();
		//_hyperParams = params.hyperParams();
		_workNXInCluster = params.workNXInCluster();
		_workContinuousX = params.workContinuousX();
		_workLogPXiGivenZi = params.workLogPXiGivenZi();
		//_workNClusInit = params.workNClusInit();
		_workLogDetTau = params.workLogDetTau();
		//_workLogDetTau00 = params.workLogDetTau00();
		//_workLogDetR1 = params.workLogDetR1();
		//_workInverseR1 = params.workInverseR1();
		_workSqrtTau = params.workSqrtTau();
	
		return *this;
	}



private:
	/// \brief The number of components
	unsigned int _NComp;

	/// \brief The number of clusters
	unsigned int _NClust;

	/// \brief Vector of probabilities of assignment to clusters
	vector<double> _logPsi;
	vector<double>_Psi;

	/// \brief A vector of Eigen dynamic vectors containing covariate means
	/// for the case of Normal covariates
	vector<VectorXd> _mu;

	/// \brief A vector of Eigen dynamic vectors containing covariate covariance
	/// matrices for the case of Normal covariates
	vector<MatrixXd> _Sigma;

	/// \brief A vector of Eigen dynamic vectors containing covariate precision
	/// matrices for the case of Normal covariates
	vector<MatrixXd> _Tau;

	/// \brief A matrix of parameters for R1 where Tau ~ Wishart (R1, kappa1)
	//MatrixXd _R1;
	//MatrixXd _Tau00;
	//MatrixXd _Sigma00;
	//VectorXd _mu00;

	//vector<MatrixXd> _R1_sparse;

	/// \brief The hyper parameter for dirichlet model
	double _alpha;

	/// \brief A vector of allocation variables
	vector<int> _z;

	/// \brief A vector containing the number of subjects in each cluster
	vector<unsigned int> _workNXInCluster;

	/// \brief A matrix containing a copy of X
	vector<vector<double> > _workContinuousX;

	/// \brief A matrix containing P(Xi|Zi)
	vector<double>  _workLogPXiGivenZi;

	/// \brief The actual number of cluster the individuals are individually
	/// allocated to. Note this is just needed for writing the log file.
	//unsigned int _workNClusInit;

	/// \brief Working vector of matrices containing matrix square root of Tau
	vector<MatrixXd> _workSqrtTau;

	/// \brief Working vector containing the log determinants of Tau
	vector<double> _workLogDetTau;

	/// \brief Working double containing the log determinants of R1
	//double _workLogDetR1;

	/// \brief Working  inverse of R1
	//MatrixXd _workInverseR1;

	//HyperParams _hyperParams;

};

