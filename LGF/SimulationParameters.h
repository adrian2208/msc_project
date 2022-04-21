#pragma once
#include <vector>
class SimulationParameters {
public:
	SimulationParameters(int NrDims, int extdofs, std::vector<int> shape, double beta, std::vector<int> cuts) {
		m_NrDims = NrDims;
		m_extdofs = extdofs;
		m_shape = shape;
		m_beta = beta;
		m_cuts = cuts;
	}

	int getNrDims() {
		return m_NrDims;
	}
	int getextdofs() {
		return m_extdofs;
	}
	int getshape(int i) {
		return m_shape[i];
	}
	double getbeta() {
		return m_beta;
	}
	int getcuts(int i) {
		return m_cuts[i];
	}
private:
	int m_NrDims;
	int m_extdofs;
	std::vector<int> m_shape;
	double m_beta;
	std::vector<int> m_cuts;
};

class FlowParameters {
public:
	FlowParameters(int ConfigurationStart, int ConfigurationStop, int flowSteps, double epsilon, int measure_every_nth_step, const std::string& dataFolder) {
		m_ConfigurationStart = ConfigurationStart;
		m_ConfigurationStop = ConfigurationStop;
		m_flowSteps = flowSteps;
		m_epsilon = epsilon;
		m_measure_every_nth_step = measure_every_nth_step;
		m_dataFolder = dataFolder;
	}

	int getConfigurationStart() {
		return m_ConfigurationStart;
	}
	int getConfigurationStop() {
		return m_ConfigurationStop;
	}
	int getflowSteps() {
		return m_flowSteps;
	}
	double getepsilon() {
		return m_epsilon;
	}
	int getMeasuringInterval() {
		return m_measure_every_nth_step;
	}
	std::string getdataFolder() {
		return m_dataFolder;
	}
private:
	int m_ConfigurationStart; 
	int m_ConfigurationStop;
	int m_flowSteps;
	double m_epsilon;
	int m_measure_every_nth_step;
	std::string m_dataFolder;
};

class HBParameters {
public:
	HBParameters(int ConfigurationStart, int ConfigurationStop, int ThermSteps, int OR_per_HB, int configSep, const std::string& dataFolder) {
		m_ConfigurationStart = ConfigurationStart;
		m_ConfigurationStop = ConfigurationStop;
		m_ThermSteps = ThermSteps;
		m_OR_per_HB = OR_per_HB;
		m_configSep = configSep;
		m_dataFolder = dataFolder;
	}

	int getConfigurationStart() {
		return m_ConfigurationStart;
	}
	int getConfigurationStop() {
		return m_ConfigurationStop;
	}
	int getThermSteps() {
		return m_ThermSteps;
	}
	int getOR_per_HB() {
		return m_OR_per_HB;
	}
	int getconfigSep() {
		return m_configSep;
	}
	std::string getdataFolder() {
		return m_dataFolder;
	}
private:
	int m_ConfigurationStart;
	int m_ConfigurationStop;
	int m_ThermSteps;
	int m_OR_per_HB;
	int m_configSep;
	std::string m_dataFolder;
};