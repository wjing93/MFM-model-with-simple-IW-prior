#pragma once
#include<iostream>
#include<fstream>
#include<string>



#include"InputData.h"

using std::string;
using std::ifstream;


void importInputData(const string& fitFilename, inputData& dataset);

void importInputData(const string& fitFilename, inputData& dataset) {

	ifstream inputFile;
	inputFile.open(fitFilename.c_str());
	if (!inputFile.is_open()) {
		std::cout << "Error opening file" << std::endl;
		//	exit(-1);
	}

	unsigned int& nSubjects = dataset.nSubjects();
	unsigned int& nCovariates = dataset.nCovariates();

	vector<vector<double> >& continuousX = dataset.continuousX();
	vector<string>& covNames = dataset.covariateNames();


	// Get the number of subjects
	inputFile >> nSubjects;
	// Get the number of covariates
	inputFile >> nCovariates;
	covNames.resize(nCovariates);

	for (unsigned int i = 0; i<nCovariates; i++) {
		inputFile >> covNames[i];
	}

	// Get the data
	continuousX.resize(nSubjects);
	vector<double> meanX(nCovariates, 0);
	for (unsigned int i = 0; i<nSubjects; i++) {

		continuousX[i].resize(nCovariates);

		for (unsigned int j = 0; j<nCovariates; j++) {
			inputFile >> continuousX[i][j];
		}
	}


	inputFile.close();

}





