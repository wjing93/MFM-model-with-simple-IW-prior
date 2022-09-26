#pragma once

#include<cmath>

class PropParams {

public:
	// Default constructor

	PropParams() {

		_alphaStdDev = 1.0;
		_alphaStdDevLower = 0.1;
		_alphaStdDevUpper = 99.9;
		_nTryAlpha = 0;
		_nAcceptAlpha = 0;
		_nLocalAcceptAlpha = 0;
		_nResetAlpha = 0;
		_alphaAcceptTarget = 0.44;
		_alphaUpdateFreq = 25;
		_alphaAnyUpdates = true;

	};

	~PropParams() {};
	

	unsigned int nTryAlpha() const {
		return _nTryAlpha;
	}

	unsigned int nAcceptAlpha() const {
		return _nAcceptAlpha;
	}

	double alphaAcceptRate() const {
		if (_nTryAlpha>0) {
			return (double)_nAcceptAlpha / (double)_nTryAlpha;
		}
		else {
			return 0.0;
		}
	}

	unsigned int alphaUpdateFreq() const {
		return _alphaUpdateFreq;
	}

	unsigned int nLocalAcceptAlpha() const {
		return _nLocalAcceptAlpha;
	}

	double alphaLocalAcceptRate() const {
		return (double)_nLocalAcceptAlpha / (double)_alphaUpdateFreq;
	}

	double alphaAcceptTarget() const {
		return _alphaAcceptTarget;
	}

	void alphaAddTry() {
		_nTryAlpha++;
	}

	void alphaAddAccept() {
		_nAcceptAlpha++;
		_nLocalAcceptAlpha++;
	}

	void alphaLocalReset() {
		_nLocalAcceptAlpha = 0;
	}

	unsigned int nResetAlpha() const {
		return _nResetAlpha;
	}

	void alphaStdDevReset() {
		_alphaStdDev = 1.0;
		_nResetAlpha++;
		_alphaStdDevLower = pow(10.0, -((double)_nResetAlpha + 1.0));
		_alphaStdDevUpper = 100.0 - pow(10.0, -((double)_nResetAlpha + 1.0));
	}


	const double alphaStdDev() const {
		return _alphaStdDev;
	}

	double& alphaStdDev() {
		return _alphaStdDev;
	}

	double alphaStdDevLower() const {
		return _alphaStdDevLower;
	}

	double alphaStdDevUpper() const {
		return _alphaStdDevUpper;
	}

	// Member function for setting the standard deviation for
	// proposal for beta for fixed effect j
	void alphaStdDev(const double& sd) {
		_alphaStdDev = sd;
	}

	bool alphaAnyUpdates() const {
		return _alphaAnyUpdates;
	}

	void alphaAnyUpdates(const bool& newStatus) {
		_alphaAnyUpdates = newStatus;
	}



	// Need to define a copy iterator
	PropParams& operator=(const PropParams& propParams) {
		_nTryAlpha = propParams.nTryAlpha();
		_nAcceptAlpha = propParams.nAcceptAlpha();
		_nLocalAcceptAlpha = propParams.nLocalAcceptAlpha();
		_nResetAlpha = propParams.nResetAlpha();
		_alphaStdDev = propParams.alphaStdDev();
		_alphaStdDevLower = propParams.alphaStdDevLower();
		_alphaStdDevUpper = propParams.alphaStdDevUpper();
		_alphaAcceptTarget = propParams.alphaAcceptTarget();
		_alphaUpdateFreq = propParams.alphaUpdateFreq();
		_alphaAnyUpdates = propParams.alphaAnyUpdates();
		return *this;

	}

private:
	unsigned int _nTryAlpha;
	unsigned int _nAcceptAlpha;
	unsigned int _nLocalAcceptAlpha;
	unsigned int _nResetAlpha;
	double _alphaStdDev;
	double _alphaStdDevLower;
	double _alphaStdDevUpper;
	double _alphaAcceptTarget;
	unsigned int _alphaUpdateFreq;
	bool _alphaAnyUpdates;

};