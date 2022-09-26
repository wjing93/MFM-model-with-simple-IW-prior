#pragma once

#include<iostream>
#include<fstream>
#include"Params.h"
#include"proposalClass.h"

using namespace std;

class outputFiles {

private:
		vector<ofstream*> outFiles;


public:

	outputFiles(string filenames){
		//file for number of clusters
		string fileName = filenames + "_nClusters.txt";
		outFiles.push_back(new ofstream(fileName.c_str()));
		//file for number of component
		fileName = filenames + "_nComponents.txt";
		outFiles.push_back(new ofstream(fileName.c_str()));

		//file for component weight
		fileName = filenames + "_psi.txt";
		outFiles.push_back(new ofstream(fileName.c_str()));

		//file for component mean and variance
		fileName = filenames + "_mu.txt";
		outFiles.push_back(new ofstream(fileName.c_str()));
		fileName = filenames + "_Sigma.txt";
		outFiles.push_back(new ofstream(fileName.c_str()));

		//file for cluster allocation
		fileName = filenames + "_z.txt";
		outFiles.push_back(new ofstream(fileName.c_str()));
		fileName = filenames + "_alpha.txt";
		outFiles.push_back(new ofstream(fileName.c_str()));
		//file for nc
		fileName = filenames + "_nMembers.txt";
		outFiles.push_back(new ofstream(fileName.c_str()));
		//acceptance probability for alpha
		fileName = filenames + "_alphaProp.txt";
		outFiles.push_back(new ofstream(fileName.c_str()));		
	};

	~outputFiles(){
		for (unsigned int i = 0; i<outFiles.size(); i++) {
			delete outFiles[i];
		}

		std::cout << "destroyed" << endl;

	};

	void writeOutput(Params& params, PropParams& propParams) {
		// File indices
		int nClustersInd = -1, nComponentsInd = -1, psiInd = -1, muInd = -1, SigmaInd = -1, zInd = -1, alphaInd = -1, nMembersInd = -1, alphaPropInd = -1;

		int r = 0;
		nClustersInd = r++;
		nComponentsInd = r++;
		psiInd = r++;
		muInd = r++;
		SigmaInd = r++;
		zInd = r++;
		alphaInd = r++;
		nMembersInd = r++;
		alphaPropInd = r++;
		unsigned int Nclust = params.NClust();
		unsigned int Ncomp = params.NComp();
		unsigned int nCovariates = params.nCovariates();
		unsigned int nSubjects = params.nSubjects();
		*(outFiles[nClustersInd]) << Nclust << endl;
		*(outFiles[nComponentsInd]) << Ncomp << endl;
		*(outFiles[alphaInd]) << params.alpha() << endl;

		for (unsigned int c = 0; c < Ncomp; c++) {
			unsigned int nXinC = params.workNXInCluster(c);
			*(outFiles[nMembersInd]) << nXinC;
			if (c < Ncomp - 1) {
				*(outFiles[nMembersInd]) << " ";
			}
			else {
				*(outFiles[nMembersInd]) << endl;
			}
		}


		for (unsigned int c = 0; c<Ncomp; c++) {
			// Print logPsi
			*(outFiles[psiInd]) << exp(params.logPsi(c));

			if (c<Ncomp - 1) {
				*(outFiles[psiInd]) << " ";
			}
			else {
				*(outFiles[psiInd]) << endl;
			}
		}

		for (unsigned int j = 0; j<nCovariates; j++) {
			for (unsigned int c = 0; c<Ncomp; c++) {
				*(outFiles[muInd]) << params.mu(c, j);
				if (c<(Ncomp - 1) || j<(nCovariates - 1)) {
					*(outFiles[muInd]) << " ";
				}
			}
		}
		*(outFiles[muInd]) << endl;

		for (unsigned int j1 = 0; j1<nCovariates; j1++) {
			for (unsigned int j2 = 0; j2<nCovariates; j2++) {
				for (unsigned int c = 0; c<Ncomp; c++) {
					*(outFiles[SigmaInd]) << params.Sigma(c, j1, j2);
					if (c<(Ncomp - 1) || j1<(nCovariates - 1) || j2<(nCovariates - 1)) {
						*(outFiles[SigmaInd]) << " ";
					}
				}
			}
		}
		*(outFiles[SigmaInd]) << endl;

		// Print the allocations
		for (unsigned int i = 0; i<nSubjects; i++) {
			*(outFiles[zInd]) << params.z(i);

			if (i<nSubjects - 1) {
				*(outFiles[zInd]) << " ";
			}
			else {

				*(outFiles[zInd]) << endl;
			}
		}

		//Print acceptance rate of alpha
		bool anyUpdates;
		anyUpdates = propParams.alphaAnyUpdates();
		if (anyUpdates) {
			*(outFiles[alphaPropInd]) << propParams.alphaAcceptRate() <<
				" " << propParams.alphaStdDev() << endl;
			propParams.alphaAnyUpdates(false);
		}

	}

	void closeFiles() {
		for (unsigned int i = 0; i<outFiles.size(); i++) {
			(*(outFiles[i])).close();
			//delete outFiles[i];
		}

		std::cout << "destroyed" << std::endl;
	}


};

