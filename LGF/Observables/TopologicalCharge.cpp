#include "TopologicalCharge.h"
#include <fstream>

TopologicalCharge::TopologicalCharge(SU3_field& U) {
	m_U = &U;
	calculate(0);
}

void TopologicalCharge::calculate(double flowTime) {
	static const double PI = 4.0 * atan(1.0);
	static const double PI2 = PI * PI;
	double localSum = 0.0;
	for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) { //This loop, alongside the multiple of eight constitutes the 4D levi cevita, thanks to the cyclicity of the trace 
		localSum += (((*m_U).Improved_fieldStrengthTensor(i, 3, 2)) * ((*m_U).Improved_fieldStrengthTensor(i, 0, 1))).ReTr();
		localSum += (((*m_U).Improved_fieldStrengthTensor(i, 3, 1)) * ((*m_U).Improved_fieldStrengthTensor(i, 2, 0))).ReTr();
		localSum += (((*m_U).Improved_fieldStrengthTensor(i, 3, 0)) * ((*m_U).Improved_fieldStrengthTensor(i, 1, 2))).ReTr();
	}
	
	double totalSum = 0.0;
	MPI_Allreduce(&localSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, mpiWrapper::comm());
	/*totalSum *= 8.0 / (4.0*32.0 * PI2);*/
	totalSum *= 1.0 / (4.0 * PI2);
	if (mpiWrapper::id() == 0) {
		m_resultVector.push_back(totalSum);
		m_FlowMeasurementTimeVector.push_back(flowTime);
		std::cout << "Topological Charge = " << totalSum << "\n";
		std::cout.flush();
	}

}

void TopologicalCharge::saveTopologicalChargeToFile(double beta, const std::string& updateMethod,const std::string& identifier, const std::string& dataFolder) {
	if (mpiWrapper::id() == 0) {
		std::string beta_str = std::to_string(beta);
		std::replace(beta_str.begin(), beta_str.end(), '.', '_');
		std::filesystem::path fieldType("/Topological_Charge/beta" + beta_str + "/");
		//move to the ensembles directory
		std::filesystem::path outPath(dataFolder + "Observables");
		//if not existing, create a directory for the field type
		outPath += fieldType;
		int i;
		for (i = 0; i < (*m_U).getLatticePtr().getNdims() - 1; i++) {
			outPath += std::to_string((*m_U).getLatticePtr().getShape()[i]) + "X";
		}
		outPath += std::to_string((*m_U).getLatticePtr().getShape()[i]) + "/" + updateMethod +"/";
		std::filesystem::create_directories(outPath);
		//write the lattice shape to the filename

		//write the lattice type and the .bin extension to the filename
		outPath += (*m_U).getLatticePtr().getType() + "_extdof" + std::to_string((*m_U).getNrExtDOF());
		outPath += identifier;
		outPath += ".csv";
		//convert the filesystem path to a format suitable for the MPI_File_open argument
		std::string outPath_string = outPath.string();
		const char* outPath_pointer = outPath_string.c_str();

		//std::ofstream ofs(outPath_string, std::ios::out | std::ofstream::binary);
		//std::ostream_iterator<char> osi{ ofs };
		//const char* beginByte = (char*)&m_resultVector[0];

		//const char* endByte = (char*)&m_resultVector.back() + sizeof(double);
		//std::copy(beginByte, endByte, osi);

		std::ofstream outFile(outPath_string,std::ios_base::app);
		//for (const auto& e : m_resultVector) outFile << e << "\n";
		int size = m_resultVector.size();
		for (int i = 0; i < size; i++) {
			outFile << m_FlowMeasurementTimeVector[i] << "," << m_resultVector[i] << "\n";
		}
		outFile.close();
	}
	MPI_Barrier(mpiWrapper::comm());
}