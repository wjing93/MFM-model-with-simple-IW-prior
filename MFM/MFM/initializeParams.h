#pragma once

// Standard includes
#include<vector>
#include<string>

// Custom includes
#include"randomDis.h"
#include"Params.h"
#include"HyperParams.h"
#include"InputData.h"


using std::vector;
using std::string;

void initialiseChain(baseGeneratorType& rndGenerator,
	Params& params, HyperParams& hyperparams, inputData & dataset, unsigned int & nCompInit) {

	unsigned int nSubjects = dataset.nSubjects();
	unsigned int nCovariates = dataset.nCovariates();

	//bool useHyperpriorR1 = options.useHyperpriorR1();
	//bool useSparsePrior = options.useSparsePrior();

	//unsigned int maxNComp = hyperparams.UpperNComp();

	// Define a uniform random number generator

	randomUniform unifRand(0, 1);

	// Copy the dataset X matrix to a working object in params
	params.workContinuousX(dataset.continuousX());

	//initialise alpha
	ramdomFisherF FRand(hyperparams.f_df1(), hyperparams.f_df2());

	double alpha = FRand(rndGenerator);
	params.alpha(alpha);

	vector<unsigned int> nXInCluster(nCompInit, 0);
	//params.workNClusInit(nClusInit);

	//initialize cluster allocation
	for (unsigned int i = 0; i<nSubjects; i++) {
		int c = (int)nCompInit*unifRand(rndGenerator);
		params.z(i, c, false);
		nXInCluster[c]++;
	}

	params.workNXInCluster(nXInCluster);
	params.RelabelClust();


	//initialise psi
	unsigned int nc;
	vector<double> alpha_psi(nCompInit);
	for (unsigned int c = 0; c < nCompInit; c++) {
		nc = params.workNXInCluster(c);

		alpha_psi[c] = alpha / (double)nCompInit + nc;

	}

	//sample from the dirichlet distribution
	vector<double> Psi = dirichletRand(rndGenerator, alpha_psi);

	params.Psi(Psi);


	//retrieve the data
	vector<VectorXd> xi(nSubjects);
	for (unsigned int i = 0; i<nSubjects; i++) {
		xi[i].setZero(nCovariates);
		for (unsigned int j = 0; j<nCovariates; j++) {
			xi[i](j) = dataset.continuousX(i, j);
		}
	}


	//initialise mu_c

	//initialise mu_c for filled components
	unsigned int nk = params.NClust();
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

		covMat = (hyperparams.Tau0() + nXInC* hyperparams.Tau0()).inverse();
		meanVec = hyperparams.Tau0()*hyperparams.mu0() + hyperparams.Tau0()*sumX[c];
		meanVec = covMat*meanVec;

		VectorXd mu(nCovariates);
		// We sample from this posterior
		mu = multivarNormalRand(rndGenerator, meanVec, covMat);

		// We store our sample
		params.mu(c, mu, false);
	}


	//initialise mu_c for empty components
	MatrixXd sigma0 = hyperparams.Tau0().inverse();
	VectorXd mu0 = hyperparams.mu0();
	for (unsigned int c = nk; c < nCompInit; c++) {

		VectorXd mu(nCovariates);
		// We sample from this posterior
		mu = multivarNormalRand(rndGenerator, mu0, sigma0);

		// We store our sample
		params.mu(c, mu, false);
	}

	
	//initialise hyperparameters involved in the component density



	//initialize Tau_c for each component

	//initialize Tau_c for filled components
	vector<MatrixXd> Rc(nk);
	for (unsigned int c = 0; c < nk; c++) {
		Rc[c].setZero(nCovariates, nCovariates);
	}

	for (unsigned int i = 0; i<nSubjects; i++) {
		unsigned int zi = params.z(i);
		Rc[zi] = Rc[zi] + (xi[i] - params.mu(zi))*((xi[i] - params.mu(zi)).transpose());
	}

	for (unsigned int c = 0; c < nk; c++) {
		Rc[c] = (Rc[c] + hyperparams.R0()).inverse();
		MatrixXd Tau = wishartRand(rndGenerator, Rc[c], params.workNXInCluster(c) + hyperparams.kappa0());

		params.Tau(c, Tau);
	}


	//initialize Tau_c for empty components
	double kappa0 = hyperparams.kappa0();
	MatrixXd R0Ivs = hyperparams.R0().inverse();

	for (unsigned int c = nk; c < nCompInit; c++) {

		MatrixXd Tau = wishartRand(rndGenerator, R0Ivs, kappa0);

		params.Tau(c, Tau);
	}


	
	

}
