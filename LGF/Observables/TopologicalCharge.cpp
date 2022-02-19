#include "TopologicalCharge.h"
#include <fstream>

TopologicalCharge::TopologicalCharge(SU3_field& U) {
	m_AutoCorrTime = calculate_AutoCorrTime();
	m_U = &U;
	calculate(0);
}

void TopologicalCharge::calculate(double flowTime) {
	static const double PI = 4.0 * atan(1.0);
	static const double PI2 = PI * PI;
	double localSum = 0.0;
	for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) { //This loop, alongside the multiple of eight constitutes the 4D levi cevita, thanks to the cyclicity of the trace 
		localSum += ((*m_U).clover_avg(i, 0, 1) * (*m_U).clover_avg(i, 2, 3)).ReTr();
		localSum -= ((*m_U).clover_avg(i, 0, 2) * (*m_U).clover_avg(i, 1, 3)).ReTr();
		localSum += ((*m_U).clover_avg(i, 0, 3) * (*m_U).clover_avg(i, 1, 2)).ReTr();
		
	}

	double totalSum = 0.0;
	MPI_Allreduce(&localSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, mpiWrapper::comm());
	totalSum *= 8.0 / (4.0*32.0 * PI2);
	if (mpiWrapper::id() == 0) {
		m_resultVector.push_back(totalSum);
		std::cout << "Topological Charge = " << totalSum << "\n";
		std::cout.flush();
	}

}






//void TopologicalCharge::calculate() {
//	static const double PI = 4.0 * atan(1.0);
//	static const double PI2 = PI * PI;
//	double localSum = 0.0;
//	//for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) {
//	//	localSum += ((*m_U).clover_avg(i, 0, 1) * (*m_U).clover_avg(i, 2, 3)).ReTr();
//	//	localSum -= ((*m_U).clover_avg(i, 0, 2) * (*m_U).clover_avg(i, 1, 3)).ReTr();
//	//	localSum += ((*m_U).clover_avg(i, 0, 3) * (*m_U).clover_avg(i, 1, 2)).ReTr();
//	//}
//
//	//for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) {
//	//	for (int mu = 0; mu < 4; mu++) {
//	//		for (int nu = 0; nu < 4; nu++) {
//	//			if (nu != mu) {
//	//				for (int rho = 0; rho < 4; rho++) {
//	//					if (rho != mu && rho != nu) {
//	//						for (int sig = 0; sig < 4; sig++) {
//	//							if (sig != mu && sig != nu && sig != rho) {
//	//								localSum += leviCivita(mu,nu,rho,sig)*((*m_U).clover_avg(i, mu, nu) * (*m_U).clover_avg(i, rho, sig)).ReTr();
//	//							}
//	//						}
//	//					}
//	//				}
//	//			}
//	//		}
//	//	}
//	//}
//	for (int i = (*m_U).Responsible_Start(); i < (*m_U).Responsible_Stop(); i++) {
//		int mu = 0;
//		for (int nu = 1; nu < 4; nu++) {
//			int rho = nu % 3;
//			rho++;
//			int sig = rho % 3;
//			sig++;
//			if (rho < sig)
//				localSum -= ((*m_U).clover_avg(i, mu, nu) * (*m_U).clover_avg(i, rho, sig)).ReTr();
//			else
//				localSum -= ((*m_U).clover_avg(i, mu, nu) * (*m_U).clover_avg(i, sig, rho).dagger()).ReTr();
//		}
//	}
//
//
//	double totalSum = 0.0;
//	MPI_Allreduce(&localSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, mpiWrapper::comm());
//	//totalSum *= 8.0 / (32.0 * PI2);
//	totalSum /= (16.0 * PI2);
//
//	if (mpiWrapper::id() == 0) {
//		m_resultVector.push_back(totalSum);
//		std::cout << totalSum << "\n";
//	}
//	//double avg_Plaquette = (*m_U).Avg_Plaquette();
//	//std::cout << "exp_val_Plaquette: " << avg_Plaquette << "\n\n";
//}

int TopologicalCharge::calculate_AutoCorrTime(){
	return 0;
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

		std::ofstream outFile(outPath_string);
		for (const auto& e : m_resultVector) outFile << e << "\n";
	}
	MPI_Barrier(mpiWrapper::comm());
}