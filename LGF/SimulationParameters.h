#pragma once
class SimulationParameters {
public:
	SimulationParameters(int NrDims, int extdofs, int shape[], double beta, int cuts[]) {
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
	int* m_shape;
	double m_beta;
	int* m_cuts;
};