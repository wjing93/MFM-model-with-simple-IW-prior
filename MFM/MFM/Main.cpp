
#include <iostream>    
#include <algorithm>    
#include <vector> 


#include"InputData.h"
#include"ReadData.h"
#include"proposals.h"
#include"Params.h"
#include"HyperParams.h"
#include"writeOutputFiles.h"
#include"proposalClass.h"
#include"initializeParams.h"

void print(const vector<unsigned int> & vec) {
	for (int i = 0; i < vec.size(); i++) {
		cout << vec[i] << " ";
	}

	cout << endl;
}

void main()
{
	
	inputData mfmData;

	string inputFilename = "output_input.txt";

	importInputData(inputFilename, mfmData);

	
	unsigned int nCovariates = mfmData.nCovariates();
	unsigned int nSubjects = mfmData.nSubjects();
	unsigned int comp_init = 50;


	//set the seed
	string seedst = "12345";
	uint_fast32_t seed = (uint_fast32_t)atoi(seedst.c_str());

	//boost random generator	
	baseGeneratorType rndGenerator;
	rndGenerator.seed(seed);


	HyperParams hyperparams;
	hyperparams.setSizes(nCovariates);
	hyperparams.setDefaults(mfmData);

	Params params;
	params.setSizes(nSubjects, nCovariates,comp_init);
	//need to initialize parameters

	initialiseChain(rndGenerator, params, hyperparams, mfmData, comp_init);


	print(params.workNXInCluster());

	
	unsigned int burnin=100;
	unsigned int sweeps=100;
	unsigned int totalSweeps = burnin + sweeps;
	
	//name of the output files


	PropParams propParams;

    
	string filenames="output";

	
	outputFiles writeout(filenames);

	for (unsigned int i = 0; i < totalSweeps; i++) {

		//First update the partition
		gibbsForZ(params, hyperparams, mfmData, rndGenerator);

		cout << "1" << endl;

		//update filled parameters in the component density
		gibbsForMuFilled(params, hyperparams, mfmData, rndGenerator);
		gibbsForTauFilled(params, hyperparams, mfmData, rndGenerator);

		cout << "2" << endl;

		//update number of compoenent
		gibbsForNComp(params, hyperparams, mfmData, rndGenerator);

		cout << "3" << endl;
		//update alpha adaptive metropolis-hastings
		AdaptiveMCMCForAlpha(params, hyperparams, mfmData, rndGenerator, propParams);

		cout << "4" << endl;

		//update empty components
		if (params.NComp()>params.NClust()) {
			gibbsForMuEmpty(params, hyperparams, mfmData, rndGenerator);
			gibbsForTauEmpty(params, hyperparams, mfmData, rndGenerator);
		}

		cout << "5" << endl;

		//update component weight
		gibbsForPsi(params, hyperparams, mfmData, rndGenerator);

		//write outputs
		if (i >= burnin) {
			writeout.writeOutput(params, propParams);
		}

		cout << "iteration " << i << " finish" << endl;

	}

	cout << "finished" << endl;
	
	writeout.closeFiles();


	cout << "finish!!!!!!!!!!!" << endl;

	writeout.~outputFiles();

	std::cin.get();

}

