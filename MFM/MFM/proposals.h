#pragma once


#include<Eigen/Dense>
#include<numeric>
#include<boost/math/distributions/normal.hpp>

#include"pdfDis.h"
#include"InputData.h"
#include"Params.h"
#include"randomDis.h"
#include"proposalClass.h"

using std::numeric_limits;

void gibbsForZ(Params& params, 
	HyperParams& hyperparams, inputData& dataset,
	baseGeneratorType& rndGenerator){

	unsigned int nSubjects = dataset.nSubjects();
	unsigned int nCovariates = dataset.nCovariates();
	unsigned int NComp = params.NComp();
	//number of components sampled added: NComp

	// Define a uniform random number generator
	randomUniform unifRand(0, 1);

	//vector that holds the number of subjects in each component
	vector<unsigned int> nMembers(NComp, 0);

	//the probability for each subject
	vector<double> rnd(nSubjects, 0.0);
	for (unsigned int i = 0; i<nSubjects; i++) {
		rnd[i] = unifRand(rndGenerator);
	}


	//likelihood times component weight for each cluster
	vector<vector<double>> logPZiGivenXi;
	logPZiGivenXi.resize(nSubjects);

	for (unsigned int i = 0; i<nSubjects; i++) {
		logPZiGivenXi[i].resize(NComp, 0.0);
		VectorXd xi = VectorXd::Zero(nCovariates);
		for (unsigned int c = 0; c<NComp; c++) {
			double logPsi=params.logPsi(c);
			if (params.z(i) == (int)c) {
				logPZiGivenXi[i][c] = params.workLogPXiGivenZi(i)+ logPsi;
			}
			else {
				for (unsigned int j = 0; j<nCovariates; j++) {
					xi(j) = params.workContinuousX(i, j);
				}
				
				logPZiGivenXi[i][c] = logPdfMultivarNormal(nCovariates, xi, params.mu(c),
					params.workSqrtTau(c), params.workLogDetTau(c))+ logPsi;
			}
		}			
	 }

	

	vector<double> meanVec(nSubjects, 0.0);


	for (unsigned int i = 0; i<nSubjects; i++) {

		vector<double> logDenom(NComp, 0.0);

		double maxlogDenom = -(numeric_limits<double>::max());

		
		for (unsigned int c = 0; c<NComp; c++) {
			logDenom[c] = logPZiGivenXi[i][c];

			if (logDenom[c]>maxlogDenom) {
				maxlogDenom = logDenom[c];
			}
		}
		//protection of very small log probability
		vector<double> Denom(NComp);
		double sumDenom = 0.0;
		for (unsigned int c = 0; c<NComp; c++) {
			double exponent= logDenom[c] - maxlogDenom;
			// Check for negative infinity (can only be negative)
			if (std::isinf(exponent) || std::isnan(exponent)) {
				exponent = -(numeric_limits<double>::max());
			}
			Denom[c] = exp(exponent);
			sumDenom += Denom[c];
		}

		//compute the cumsum
		vector<double> cumDenom(NComp);
		for (unsigned int c = 0; c<NComp; c++) {
			Denom[c] /= sumDenom;
			if (c == 0) {
				cumDenom[c] = Denom[c];
			}
			else {
				cumDenom[c] = cumDenom[c - 1] + Denom[c];
			}
			
		}
		
		unsigned int zi;

		if (NComp == 1) {
			zi = 0;
		}
		else {
			zi = 0;
			for (unsigned int c = 0; c<NComp; c++) {
				if (rnd[i]<cumDenom[c]) {
					zi = c;
					break;
				}
			}
		}


		params.z(i, zi, true);

		nMembers[zi]++;
	
	}

	params.workNXInCluster(nMembers);

	params.RelabelClust();

	
}




void gibbsForMuFilled(Params& params, HyperParams& hyperparams, inputData& dataset, baseGeneratorType& rndGenerator) {

	// Find the number of clusters
	unsigned int nk = params.NClust();
	// Find the number of covariates
	unsigned int nCovariates = dataset.nCovariates();
	// Find the number of subjects
	unsigned int nSubjects = dataset.nSubjects();


	// In the following it is useful to have the rows of X as
	// Eigen dynamic vectors
	vector<VectorXd> xi(nSubjects);
	for (unsigned int i = 0; i<nSubjects; i++) {
		xi[i].setZero(nCovariates);
		for (unsigned int j = 0; j<nCovariates; j++) {
			xi[i](j) = dataset.continuousX(i, j);
		}
	}



	// We begin by computing the mean X for individuals in each cluster
	vector<VectorXd> sumX(nk);
	for (unsigned int c = 0; c < nk; c++) {
		sumX[c].setZero(nCovariates);
	}

	for (unsigned int i = 0; i<nSubjects; i++) {
		sumX[params.z(i)] = sumX[params.z(i)] + xi[i];
	}


	for (unsigned int c = 0; c < nk; c++) {
		// Having computed this we can calcuate the posterior mean
		// and posterior covariance for each mu_c
		int nXInC = params.workNXInCluster(c);
		

		MatrixXd covMat(nCovariates, nCovariates);
		VectorXd meanVec(nCovariates);

		covMat = (hyperparams.Tau0() + nXInC* params.Tau(c)).inverse();
		meanVec = hyperparams.Tau0()*hyperparams.mu0() + params.Tau(c)*sumX[c];
		meanVec = covMat*meanVec;

		VectorXd mu(nCovariates);
		// We sample from this posterior
		mu = multivarNormalRand(rndGenerator, meanVec, covMat);

		// We store our sample
		params.mu(c, mu, true);
	}
}


void gibbsForTauFilled(Params& params, HyperParams& hyperparams, inputData& dataset, baseGeneratorType& rndGenerator) {

	// Find the number of clusters
	unsigned int nk = params.NClust();
	// Find the number of covariates
	unsigned int nCovariates = dataset.nCovariates();
	// Find the number of subjects
	unsigned int nSubjects = dataset.nSubjects();


	// In the following it is useful to have the rows of X as
	// Eigen dynamic vectors
	vector<VectorXd> xi(nSubjects);
	for (unsigned int i = 0; i<nSubjects; i++) {
		xi[i].setZero(nCovariates);
		for (unsigned int j = 0; j<nCovariates; j++) {
			xi[i](j) = dataset.continuousX(i, j);
		}
	}

	vector<MatrixXd> Rc(nk);
	for (unsigned int c = 0; c < nk; c++) {
		Rc[c].setZero(nCovariates, nCovariates);
	}

	for (unsigned int i = 0; i<nSubjects; i++) {
		unsigned int zi = params.z(i);
		Rc[zi] = Rc[zi] + (xi[i] - params.mu(zi))*((xi[i] - params.mu(zi)).transpose());
	}

	for (unsigned int c = 0; c < nk; c++) {
		Rc[c] = (Rc[c]+ hyperparams.R0()).inverse();
		MatrixXd Tau = wishartRand(rndGenerator, Rc[c], params.workNXInCluster(c) + hyperparams.kappa0());

		params.Tau(c, Tau);
	}

}

double logCompPrior(unsigned int k, HyperParams& hyperparams) {

	double out;

	double Comp_ck = hyperparams.Comp_ck();
	double Comp_ak = hyperparams.Comp_ak();
	double Comp_bk = hyperparams.Comp_bk();

	out = lgamma(Comp_ck + (double)k - 1.0) + lgamma(Comp_ck + Comp_ak) + lgamma((double)k - 1.0 + Comp_bk) - lgamma(Comp_ck + Comp_ak + (double)k - 1.0 + Comp_bk)
		- (lgamma(Comp_ck) + lgamma((double)k) + lgamma(Comp_ak) + lgamma(Comp_bk) - lgamma(Comp_ak+ Comp_bk));

	return out;


}

double logCompCondi(unsigned int k, HyperParams& hyperparams, Params& params) {

	double out;
	double alpha = params.alpha();
	unsigned int nk = params.NClust();

	double alphak = alpha / (double)k;

	out = (double)nk * log(alpha) + log(factorial<double>(k)) - ((double)nk*log((double)k) + log(factorial<double>(k-nk)));

	for (unsigned int c = 0; c < nk; c++) {
		int nc = params.workNXInCluster(c);
		out = out + lgamma((double)nc +alphak) - lgamma(1.0+alphak);
	}

	return out;



}


void gibbsForNComp(Params& params, HyperParams& hyperparams, inputData& dataset, baseGeneratorType& rndGenerator) {

	// Find the number of clusters
	unsigned int nk = params.NClust();
	// Find the number of covariates
	unsigned int nCovariates = dataset.nCovariates();
	// Find the number of subjects
	unsigned int nSubjects = dataset.nSubjects();              
	//Find the upper bound of the number of components
	unsigned int UpperCom = hyperparams.UpperNComp();

	//the range of the number of component 
	unsigned int NCompRange = UpperCom - nk + 1;
	VectorXd logCompProp=VectorXd::Zero(NCompRange);

	randomUniform unifRand(0, 1);

	double rnd = unifRand(rndGenerator);

	//numerical protection of very large log values
	double maxlogCompProp = -(numeric_limits<double>::max());

	for (unsigned int c = nk; c <= UpperCom; c++) {

		int indc = c - nk;

		logCompProp(indc) = logCompPrior(c, hyperparams) + logCompCondi(c, hyperparams, params);

		if (logCompProp(indc)>maxlogCompProp) {
			maxlogCompProp = logCompProp(indc);
		}
	}

	vector<double> CompProp(NCompRange);
	double sumCompProp = 0.0;

	for (unsigned int c = 0; c < NCompRange; c++) {
		double exponent = logCompProp[c] - maxlogCompProp;

		if (std::isinf(exponent) || std::isnan(exponent)) {
			exponent = -(numeric_limits<double>::max());
		}	

		CompProp[c] = exp(exponent);
		sumCompProp += exp(exponent);
	}

	vector<double> cumSum(NCompRange);
	for (unsigned int c = 0; c < NCompRange; c++) {
		CompProp[c] /= sumCompProp;
		if (c == 0) {
			cumSum[c] = CompProp[c];
		}
		else {
			cumSum[c] = cumSum[(c - 1)] + CompProp[c];
		}
	}

	unsigned int newComp=nk;
	for (unsigned int c = 0; c < NCompRange; c++) {
		if (cumSum[c]>= rnd) {
			newComp = newComp + c;
			break;
		}

	}

	params.NComp(newComp);


}

double logTargetDis(double alpha, HyperParams& hyperparams, Params& params, unsigned int samplesize) {
	double out;
	double f_df1 = hyperparams.f_df1();
	double f_df2 = hyperparams.f_df2();
	unsigned int nk = params.NClust();
	unsigned int ncomp = params.NComp();
	double alphak = alpha / (double)ncomp;

	out = (f_df1 / 2.0 - 1.0)*log(alpha) - (f_df1 + f_df2) / 2.0*log(1.0 + f_df1 / f_df2*alpha);
	out = out + (double)nk*log(alpha) + lgamma(alpha) - lgamma(alpha + (double)samplesize);
	for (unsigned int c = 0; c < nk; c++) {
		int nc = params.workNXInCluster(c);
		out = out + lgamma((double)nc + alphak) - lgamma(1.0 + alphak);
	}

	return out;

}


void AdaptiveMCMCForAlpha(Params& params, HyperParams& hyperparams, inputData& dataset, baseGeneratorType& rndGenerator,
	                      PropParams& propParams) {

	unsigned int nCovariates = dataset.nCovariates();
	unsigned int nSubjects = dataset.nSubjects();

	randomUniform unifRand(0, 1);
	double alpha_curr = params.alpha();
	double stdDev = propParams.alphaStdDev();
	double var_prop = stdDev*stdDev;
	double log_alpha_prop = NormalRand(rndGenerator, log(alpha_curr), var_prop);
	double alpha_prop = exp(log_alpha_prop);

	double logAcceptRatio = 0.0;



	//Add the target distribution
	logAcceptRatio += logTargetDis(alpha_prop, hyperparams, params, nSubjects);
	logAcceptRatio -= logTargetDis(alpha_curr, hyperparams, params, nSubjects);

	// Add the proposal contribution
	logAcceptRatio += logPdflogNormal(alpha_curr, log(alpha_prop), stdDev);
	logAcceptRatio -= logPdflogNormal(alpha_prop, log(alpha_curr), stdDev);

	propParams.alphaAddTry();
	if (unifRand(rndGenerator)<exp(logAcceptRatio)) {

		propParams.alphaAddAccept();
		// If the move was accepted update the state
		params.alpha(alpha_prop);

		// Also update the proposal standard deviation
		if (propParams.nTryAlpha() % propParams.alphaUpdateFreq() == 0) {
			stdDev += 10 * (propParams.alphaLocalAcceptRate() - propParams.alphaAcceptTarget()) /
				pow((double)(propParams.nTryAlpha() / propParams.alphaUpdateFreq()) + 2.0, 0.75);
			propParams.alphaAnyUpdates(true);
			if (stdDev>propParams.alphaStdDevUpper() || stdDev<propParams.alphaStdDevLower()) {
				propParams.alphaStdDevReset();
			}
			propParams.alphaLocalReset();
		}
	}
	else {
		// Otherwise update the proposal standard deviation
		if (propParams.nTryAlpha() % propParams.alphaUpdateFreq() == 0) {
			stdDev += 10 * (propParams.alphaLocalAcceptRate() - propParams.alphaAcceptTarget()) /
				pow((double)(propParams.nTryAlpha() / propParams.alphaUpdateFreq()) + 2.0, 0.75);
			propParams.alphaAnyUpdates(true);
			if (stdDev>propParams.alphaStdDevUpper() || stdDev<propParams.alphaStdDevLower()) {
				propParams.alphaStdDevReset();
			}
			propParams.alphaLocalReset();
		}

	}
}



//Only include if NComp is greater than NClust
void gibbsForMuEmpty(Params& params, HyperParams& hyperparams, inputData& dataset, baseGeneratorType& rndGenerator) {

	// Find the number of clusters
	unsigned int nk = params.NClust();
	unsigned int nComp = params.NComp();
	// Find the number of covariates
	unsigned int nCovariates = dataset.nCovariates();
	MatrixXd sigma0 = hyperparams.Tau0().inverse();
	VectorXd mu0 = hyperparams.mu0();

	for (unsigned int c = nk; c < nComp; c++) {

		VectorXd mu(nCovariates);
		// We sample from this posterior
		mu = multivarNormalRand(rndGenerator, mu0, sigma0);

		// We store our sample
		params.mu(c, mu, true);
	}
}

//Only include if NComp is greater than NClust
void gibbsForTauEmpty(Params& params, HyperParams& hyperparams, inputData& dataset, baseGeneratorType& rndGenerator) {

	// Find the number of clusters
	unsigned int nk = params.NClust();
	unsigned int nComp = params.NComp();
	// Find the number of covariates
	unsigned int nCovariates = dataset.nCovariates();
	double kappa0 = hyperparams.kappa0();
	MatrixXd R0Ivs = hyperparams.R0().inverse();

	for (unsigned int c = nk; c < nComp; c++) {

		MatrixXd Tau = wishartRand(rndGenerator, R0Ivs, kappa0);

		params.Tau(c, Tau);
	}

}

void gibbsForPsi(Params& params, HyperParams& hyperparams, inputData& dataset, baseGeneratorType& rndGenerator) {

	unsigned int nComp = params.NComp();
	unsigned int nc;
	double alpha = params.alpha();
	vector<double> alpha_psi(nComp);
	for (unsigned int c = 0; c < nComp; c++) {
		nc= params.workNXInCluster(c);

		alpha_psi[c] = alpha / (double)nComp + nc;

	}

	//sample from the dirichlet distribution
	vector<double> Psi = dirichletRand(rndGenerator, alpha_psi);

	params.Psi(Psi);
}
