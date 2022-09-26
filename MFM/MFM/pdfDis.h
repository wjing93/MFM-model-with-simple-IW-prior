#pragma once
#include<cmath>
#include<string>

#include<boost/math/distributions/normal.hpp>
#include<boost/math/special_functions/gamma.hpp>
#include<boost/math/constants/constants.hpp>

#include<Eigen/Core>
#include<Eigen/Cholesky>
#include<Eigen/LU>

#include"MathFun.h"

using boost::math::lgamma;
using std::string;

using namespace Eigen;
using namespace boost::math::constants;

double logPdfBernoulli(const unsigned int& x, const double& p) {
	return x*log(p) + (1 - x)*log(1.0 - p);
}

double logPdfBeta(const double& x, const double& a, const double& b) {
	return lgamma(a + b) - lgamma(a) - lgamma(b) + (a - 1)*log(x) + (b - 1)*log(1 - x);
}

double logPdfBinomial(const unsigned int& x, const unsigned int& n, const double& p) {
	return (double)x*log(p) + (double)(n - x)*log(1.0 - p)
		+ lgamma((double)(n + 1)) - lgamma((double)(x + 1)) - lgamma((double)(n - x + 1));
}

double logPdfDirichlet(const vector<double>& x, const vector<double>& alpha, const bool& xOnLogScale) {

	double out = 0.0;
	unsigned int K = x.size();
	double sumAlpha = 0.0;
	for (unsigned int k = 0; k<K; k++) {
		double logxk = x[k];
		if (!xOnLogScale) {
			logxk = log(logxk);
		}

		out += (alpha[k] - 1)*logxk - lgamma(alpha[k]);
		sumAlpha += alpha[k];
	}
	out += lgamma(sumAlpha);

	return out;
}

double logPdfExponential(const double& x, const double& lambda) {
	return log(lambda) - lambda*x;
}

double logPdfLaplace(const double& x, const double& lambda) {
	return log(lambda / 2.0) - lambda*abs(x);
}

double logPdfGamma(const double& x, const double& shape, const double& rate) {
	return shape*log(rate) + (shape - 1)*log(x) - rate*x - lgamma(shape);
}

double logPdfInverseGamma(const double& x, const double& shape, const double& scale) {
	return shape*log(scale) - (shape + 1)*log(x) - scale / x - lgamma(shape);
}

double logPdfNormal(const double& x, const double& mu, const double& sigma) {
	return -0.5*log(2 * pi<double>()) - log(sigma) - 0.5*(x - mu)*(x - mu) / (sigma*sigma);
}

double logPdflogNormal(const double& x, const double& mu, const double& sigma) {
	return -0.5*log(2 * pi<double>()) - log(sigma)-log(x) - 0.5*(log(x) - mu)*(log(x) - mu) / (sigma*sigma);
}

double logPdfQuantile(const double& x, const double& mu, const double& sigma, const double& pQuantile) {
	return log(pQuantile) + log(1 - pQuantile) - log(sigma) - (abs((x - mu) / sigma) + (2 * pQuantile - 1)*(x - mu) / sigma) / 2;
}

double logPdfPoisson(const unsigned int &x, const double& mu) {
	return (double)x*log(mu) - mu - lgamma((double)x + 1);
}

double logPdfStudentsT(const double &x, const unsigned int& dof) {
	return -0.5*log((double)dof*pi<double>()) - (((double)dof + 1) / 2.0)*log(1.0 + x*x / (double)dof)
		+ lgamma((double)(dof + 1) / 2.0) - lgamma((double)dof / 2.0);
}

double logPdfLocationScaleT(const double &x, const double& mu,
	const double& sigma, const unsigned int& dof) {
	// http://www.mathworks.com/help/toolbox/stats/brn2ivz-145.html
	// location param=mu, scale param=sigma, shape param=dof
	return logPdfStudentsT((x - mu) / sigma, dof);

}

double logPdfTruncatedNormal(const double& x, const double& mean, const double& stdDev,
	const string& distType, const double& lower, const double& upper) {

	boost::math::normal_distribution<double> normDist(mean, stdDev);
	double pLower, pUpper;

	if (distType.compare("U") == 0) {
		pLower = 0.0000000001;
		pUpper = cdf(normDist, upper);
	}
	else if (distType.compare("L") == 0) {
		pLower = cdf(normDist, lower);
		pUpper = 1.0 - 0.0000000001;
	}
	else {
		pLower = cdf(normDist, lower);
		pUpper = cdf(normDist, upper);
	}

	return log(pdf(normDist, x)) - log(pUpper - pLower);

}

double logPdfMultinomialSizeOne(const unsigned int& x, const vector<double>& p) {
	return log(p[x]);
}


double logPdfMultivarNormal(const unsigned int& sizeX, const VectorXd& x, const VectorXd& meanVec, const MatrixXd& sqrtPrecMat, const double& logDetPrecMat) {
	// If S is the (upper triangular matrix) square root of the precision P, then S'S=P

	VectorXd work = sqrtPrecMat*(x - meanVec);
	double tmp = work.squaredNorm();
	return -0.5*((double)sizeX*log(2 * pi<double>()) - logDetPrecMat + tmp);
}


//used in the log prior case and the argument is SigmalogVec
double logPdfMultivarT(const unsigned int& sizeX, const VectorXd& x, const VectorXd& meanVec, const MatrixXd& sqrtPrecMat, const double& logDetPrecMat, const double df) {
	// If S is the (upper triangular matrix) square root of the precision P, then S'S=P
	double out;
	VectorXd work = sqrtPrecMat*(x - meanVec);
	double tmp = work.squaredNorm();
	out = lgamma((df + sizeX) / 2.0) - lgamma(df / 2.0) - sizeX / 2.0*log(df*pi<double>()) - 0.5*logDetPrecMat - ((df + sizeX) / 2.0)*log(1 + tmp / df);
	return out;
}

//when separation strategy is used to sample sigma_c
double logPdfMultivarNormalSS(const unsigned int& sizeX, const VectorXd& x, const VectorXd& meanVec, const MatrixXd& TauS, const double& logDetTauS, const MatrixXd& sqrtPrecMat, const double& logDetPrecMat) {
	// If S is the (upper triangular matrix) square root of the precision P, then S'S=P

	VectorXd work = sqrtPrecMat*TauS*(x - meanVec);
	double tmp = work.squaredNorm();
	return -0.5*((double)sizeX*log(2 * pi<double>()) - logDetPrecMat - 2 * logDetTauS + tmp);
}

double logPdfWishart(const unsigned int& dimA, const MatrixXd& A, const double& logDetA, const MatrixXd& inverseR, const double& logDetR, const double& kappa) {
	MatrixXd work = inverseR*A;
	return -0.5*(kappa*logDetR - (kappa - (double)dimA - 1)*logDetA + work.trace())
		- (0.5*kappa*(double)dimA*log(2.0) + logMultivarGammaFn(kappa / 2.0, dimA));
}

double logPdfInverseWishart(const unsigned int& dimA, const MatrixXd& inverseA, const double& logDetA, const MatrixXd& R, const double& logDetR, const double& kappa) {
	MatrixXd work = R*inverseA;
	return 0.5*(kappa*logDetR - (kappa + (double)dimA + 1)*logDetA - work.trace())
		- (0.5*kappa*(double)dimA*log(2.0) + logMultivarGammaFn(kappa / 2.0, dimA));
}

