#pragma once

// Standard includes
#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>


using std::vector;
using std::ifstream;
using std::string;

/// \class pReMiuMData PReMiuMData.h "PReMiuMModel.h"
/// \brief A class for PReMiuMpp Data
class inputData {

public:
	/// \brief Default constructor
	inputData() : _nSubjects(0), _nCovariates(0) {};

	/// \brief Default destructor
	~inputData() {};

	/// \brief Return the number of subjects
	unsigned int nSubjects() const;

	/// \brief Return the number of subjects
	unsigned int& nSubjects();

	/// \brief Set the number of subjects
	void nSubjects(const unsigned int& nSubj);

	unsigned int size() const;

	/// \brief Return the number of covariates
	unsigned int nCovariates() const;

	/// \brief Return the number of covariates
	unsigned int& nCovariates();

	/// \brief Set the number of covariates
	void nCovariates(const unsigned int& nCov);

	/// \brief Return the vector of the covariate names
	vector<string> covariateNames() const;

	/// \brief Return the vector of the covariate names
	vector<string>& covariateNames();


	/// \brief Set the vector of the covariate names
	void covariateNames(const vector<string>& covNames);

	/// \brief Return name for covariate j
	string covariateNames(const unsigned int& j) const;

	/// \brief Return the covariate matrix
	const vector<vector<double> >& continuousX() const;

	/// \brief Return the covariate matrix
	vector<vector<double> >& continuousX();
	/// \brief Return the jth covariate for subject i
	double continuousX(const unsigned int& i, const unsigned int& j) const;

	/// \brief Set the jth covariate for subject i
	void continuousX(const unsigned int& i, const unsigned int& j, const double& x);


private:
	/// \brief The number of subjects
	unsigned int _nSubjects;

	/// \brief The number of covariates
	unsigned int _nCovariates;

	/// \brief A matrix (vector of vectors) of the covariate data
	/// \note this is a signed int because missing values are typically stored
	/// as negative values
	vector<vector<double> > _continuousX;

	/// \brief A vector of covariate names
	vector<string> _covariateNames;



};

unsigned int inputData::nSubjects() const {
	return _nSubjects;
}

/// \brief Return the number of subjects
unsigned int& inputData::nSubjects() {
	return _nSubjects;
}

/// \brief Set the number of subjects
void inputData::nSubjects(const unsigned int& nSubj) {
	_nSubjects = nSubj;
}

unsigned int inputData::size() const {
	return _nSubjects;
}

/// \brief Return the number of covariates
unsigned int inputData::nCovariates() const {
	return _nCovariates;
}

/// \brief Return the number of covariates
unsigned int& inputData::nCovariates() {
	return _nCovariates;
}

/// \brief Set the number of covariates
void inputData::nCovariates(const unsigned int& nCov) {
	_nCovariates = nCov;
}

/// \brief Return the vector of the covariate names
vector<string> inputData::covariateNames() const {
	return _covariateNames;
}

/// \brief Return the vector of the covariate names
vector<string>& inputData::covariateNames() {
	return _covariateNames;
}


/// \brief Set the vector of the covariate names
void inputData::covariateNames(const vector<string>& covNames) {
	_covariateNames.clear();
	_covariateNames.resize(covNames.size());
	_covariateNames.insert(_covariateNames.begin(), covNames.begin(), covNames.end());
}

/// \brief Return name for covariate j
string inputData::covariateNames(const unsigned int& j) const {
	return _covariateNames[j];
}

/// \brief Return the covariate matrix
const vector<vector<double> >& inputData::continuousX() const {
	return _continuousX;
}

/// \brief Return the covariate matrix
vector<vector<double> >& inputData::continuousX() {
	return _continuousX;
}

/// \brief Return the jth covariate for subject i
double inputData::continuousX(const unsigned int& i, const unsigned int& j) const {
	return _continuousX[i][j];
}

/// \brief Set the jth covariate for subject i
void inputData::continuousX(const unsigned int& i, const unsigned int& j, const double& x) {
	_continuousX[i][j] = x;
}

