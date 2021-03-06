#include "SU3_field.h"
#include <algorithm>

SU3_field::SU3_field(Lattice& lattice, int NrExtDOF) : Field(lattice, NrExtDOF){
	m_FieldArray_NrBytes = 8*2*9 * NrExtDOF * lattice.m_thisProc_Volume;//9 entries of C_double each containing 16 bytes
}

void SU3_field::saveSU3ToFileDEPRECATED(double beta, const std::string& updateMethod, const std::string& dataFolder, const std::string& identifier) {
	(*this).transfer_FieldValues();
	std::string beta_str = std::to_string(beta);
	std::replace(beta_str.begin(), beta_str.end(), '.', '_');
	n_fs::path fieldType("/su3/beta"+beta_str+"/");
	int result;
	MPI_File file;
	MPI_Status status;
	MPI_Offset displacement = 0;
	//move to the ensembles directory
	n_fs::path outPath(dataFolder+"ensembles");
	//if not existing, create a directory for the field type
	outPath += fieldType;
	int i;
	for (i = 0; i < m_lattice->getNdims() - 1; i++) {
		outPath += std::to_string(m_lattice->getShape()[i]) + "X";
	}
	outPath += std::to_string(m_lattice->getShape()[i]) + "/" + updateMethod + "/";
	n_fs::create_directories(outPath);
	//write the lattice shape to the filename

	//write the lattice type and the .bin extension to the filename
	outPath +=  m_lattice->getType() + "_extdof" + std::to_string(m_NrExtDOF) + "_Nprocs" + std::to_string(mpiWrapper::nProcs()) + "_";
	outPath += identifier;
	outPath += ".bin";
	//convert the filesystem path to a format suitable for the MPI_File_open argument
	std::string outPath_string = outPath.string();
	const char* outPath_pointer = outPath_string.c_str();

	//Divide the size of FieldArray in bytes by the size of a char
	//then tell MPI that the datatype you're using is char.
	// This way it allocates enough bytes regardless of the size of T
	MPI_Datatype dataType = MPI_CHAR;
	int NrDatatype = m_FieldArray_NrBytes / sizeof(char);

	result = MPI_File_open(mpiWrapper::comm(), outPath_pointer, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
	if (result != MPI_SUCCESS) {
		char Error_string[MPI_MAX_ERROR_STRING];
		int Error_class;
		int Error_len;
		MPI_Error_class(result, &Error_class);
		MPI_Error_string(result, Error_string, &Error_len);
		std::cout << Error_class << " Error: " << Error_string << std::endl;
	}
	char datarep[7] = "native";
	result = MPI_File_set_view(file, displacement, dataType, dataType, datarep,MPI_INFO_NULL);
	if (result != MPI_SUCCESS) {
		char Error_string[MPI_MAX_ERROR_STRING];
		int Error_class;
		int Error_len;
		MPI_Error_class(result, &Error_class);
		MPI_Error_string(result, Error_string, &Error_len);
		std::cout << Error_class << " Error: " << Error_string << std::endl;
	}

	//the displacement in bytes sets the offset between the field stored by 
	//process n and process n-1
	displacement = mpiWrapper::id();
	displacement *= m_FieldArray_NrBytes;
	if (mpiWrapper::id() == 0) {
		std::cout << "Saving SU3 Field onto " << std::to_string(mpiWrapper::nProcs()) << " processes\nTo path " << outPath_string << "\n";
		std::cout.flush();
	}
	//The file is written
	for (int site = 0; site < m_lattice->getthisProc_Volume(); site++) {
		for (int extDOF = 0; extDOF < m_NrExtDOF; extDOF++) {
			result = MPI_File_write_at(file, displacement, FieldArray[site*m_NrExtDOF+extDOF].getMemPointer(), 144, dataType, &status);
			displacement += 144;
			if (result != MPI_SUCCESS) {
				char Error_string[MPI_MAX_ERROR_STRING];
				int Error_class;
				int Error_len;
				MPI_Error_class(result, &Error_class);
				MPI_Error_string(result, Error_string, &Error_len);
				std::cout << Error_class << " Error: " << Error_string << std::endl;
			}
		}
	}

	MPI_File_close(&file);
}
void SU3_field::saveSU3ToFile(double beta, const std::string& updateMethod, const std::string& dataFolder, const std::string& identifier) {
	(*this).transfer_FieldValues();
	std::string beta_str = std::to_string(beta);
	std::replace(beta_str.begin(), beta_str.end(), '.', '_');
	n_fs::path fieldType("/su3/beta" + beta_str + "/");
	int result;
	MPI_File file;
	MPI_Status status;
	MPI_Offset displacement = 0;
	//move to the ensembles directory
	n_fs::path outPath(dataFolder + "ensembles");
	//if not existing, create a directory for the field type
	outPath += fieldType;
	int i;
	for (i = 0; i < m_lattice->getNdims() - 1; i++) {
		outPath += std::to_string(m_lattice->getShape()[i]) + "X";
	}
	outPath += std::to_string(m_lattice->getShape()[i]) + "/" + updateMethod + "/";
	n_fs::create_directories(outPath);
	//write the lattice shape to the filename

	//write the lattice type and the .bin extension to the filename
	outPath += m_lattice->getType() + "_extdof" + std::to_string(m_NrExtDOF) + "_";
	outPath += identifier;
	outPath += ".bin";
	//convert the filesystem path to a format suitable for the MPI_File_open argument
	std::string outPath_string = outPath.string();
	const char* outPath_pointer = outPath_string.c_str();

	//Divide the size of FieldArray in bytes by the size of a char
	//then tell MPI that the datatype you're using is char.
	// This way it allocates enough bytes regardless of the size of T
	MPI_Datatype dataType = MPI_CHAR;
	int NrDatatype = 8 * 2 * 9 * m_NrExtDOF * m_lattice->m_responsible_Volume / sizeof(char);

	result = MPI_File_open(mpiWrapper::comm(), outPath_pointer, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
	if (result != MPI_SUCCESS) {
		char Error_string[MPI_MAX_ERROR_STRING];
		int Error_class;
		int Error_len;
		MPI_Error_class(result, &Error_class);
		MPI_Error_string(result, Error_string, &Error_len);
		std::cout << Error_class << " Error: " << Error_string << std::endl;
	}
	char datarep[7] = "native";
	result = MPI_File_set_view(file, displacement, dataType, dataType, datarep, MPI_INFO_NULL);
	if (result != MPI_SUCCESS) {
		char Error_string[MPI_MAX_ERROR_STRING];
		int Error_class;
		int Error_len;
		MPI_Error_class(result, &Error_class);
		MPI_Error_string(result, Error_string, &Error_len);
		std::cout << Error_class << " Error: " << Error_string << std::endl;
	}

	//the displacement in bytes sets the offset between the field stored by 
	//process n and process n-1
	//displacement = mpiWrapper::id();
	//displacement *= 8 * 2 * 9 * m_NrExtDOF * m_lattice->m_responsible_Volume;
	if (mpiWrapper::id() == 0) {
		std::cout << "Saving SU3 Field to path " << outPath_string << "\n";
		std::cout.flush();
	}
	//The file is written
	for (int site = (*this).Responsible_Start(); site < (*this).Responsible_Stop(); site++) {
		for (int extDOF = 0; extDOF < m_NrExtDOF; extDOF++) {
			//std::cout << "id: " << mpiWrapper::id() << " location: " << m_lattice->m_InternalToTotal_idx[site] << std::endl;
			displacement = (m_lattice->m_InternalToTotal_idx[site] * m_NrExtDOF + extDOF)*144;
			result = MPI_File_write_at(file, displacement, FieldArray[site * m_NrExtDOF + extDOF].getMemPointer(), 144, dataType, &status);
			//result = MPI_File_write_at(file, displacement, FieldArray[site * m_NrExtDOF + extDOF].getMemPointer(), 144, dataType, &status);
			//displacement += 144;
			if (result != MPI_SUCCESS) {
				char Error_string[MPI_MAX_ERROR_STRING];
				int Error_class;
				int Error_len;
				MPI_Error_class(result, &Error_class);
				MPI_Error_string(result, Error_string, &Error_len);
				std::cout << Error_class << " Error: " << Error_string << std::endl;
			}
		}
	}

	MPI_File_close(&file);
}
void SU3_field::loadSU3FromFileDEPRECATED(double beta, const std::string& updateMethod, const std::string& dataFolder, const std::string& identifier) {

	std::string beta_str = std::to_string(beta);
	std::replace(beta_str.begin(), beta_str.end(), '.', '_');
	n_fs::path fieldType("/su3/beta" + beta_str + "/");
	int result;
	MPI_File file;
	MPI_Status status;
	MPI_Offset displacement = 0;
	//move to the ensembles directory
	n_fs::path outPath(dataFolder+"ensembles");
	//if not existing, create a directory for the field type
	outPath += fieldType;
	int i;
	for (i = 0; i < m_lattice->getNdims() - 1; i++) {
		outPath += std::to_string(m_lattice->getShape()[i]) + "X";
	}
	outPath += std::to_string(m_lattice->getShape()[i]) + "/" + updateMethod + "/";
	n_fs::create_directories(outPath);
	//write the lattice shape to the filename

	//write the lattice type and the .bin extension to the filename
	outPath += m_lattice->getType() + "_extdof" + std::to_string(m_NrExtDOF) + "_Nprocs" + std::to_string(mpiWrapper::nProcs()) + "_";
	outPath += identifier;
	outPath += ".bin";
	//convert the filesystem path to a format suitable for the MPI_File_open argument
	std::string outPath_string = outPath.string();
	const char* outPath_pointer = outPath_string.c_str();
	if (mpiWrapper::id() == 0) {
		std::cout << "Loading SU3 Field onto " << std::to_string(mpiWrapper::nProcs()) << " processes\nFrom path " << outPath_string << "\n";
		std::cout.flush();
	}
	
	//Divide the size of FieldArray in bytes by the size of a char
	//then tell MPI that the datatype you're using is char.
	// This way it allocates enough bytes regardless of the size of T
	MPI_Datatype dataType = MPI_CHAR;
	int NrDatatype = m_FieldArray_NrBytes / sizeof(char);

	result = MPI_File_open(mpiWrapper::comm(), outPath_pointer, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	if (result != MPI_SUCCESS) {
		char Error_string[MPI_MAX_ERROR_STRING];
		int Error_class;
		int Error_len;
		MPI_Error_class(result, &Error_class);
		MPI_Error_string(result, Error_string, &Error_len);
		std::cout << Error_class <<" Error: " << Error_string << std::endl;
	}
	char datarep[7] = "native";
	result = MPI_File_set_view(file, displacement, dataType, dataType, datarep, MPI_INFO_NULL);
	if (result != MPI_SUCCESS) {
		char Error_string[MPI_MAX_ERROR_STRING];
		int Error_class;
		int Error_len;
		MPI_Error_class(result, &Error_class);
		MPI_Error_string(result, Error_string, &Error_len);
		std::cout << Error_class << " Error: " << Error_string << std::endl;
	}

	//the displacement in bytes sets the offset between the field stored by 
	//process n and process n-1
	displacement = mpiWrapper::id();
	displacement *= m_FieldArray_NrBytes;

	//The file is written
	for (int site = 0; site < m_lattice->getthisProc_Volume(); site++) {
		for (int extDOF = 0; extDOF < m_NrExtDOF; extDOF++) {
			result = MPI_File_read_at(file, displacement, FieldArray[site * m_NrExtDOF + extDOF].getMemPointer(), 144, dataType, &status);
			displacement += 144;
			if (result != MPI_SUCCESS) {
				char Error_string[MPI_MAX_ERROR_STRING];
				int Error_class;
				int Error_len;
				MPI_Error_class(result, &Error_class);
				MPI_Error_string(result, Error_string, &Error_len);
				std::cout << Error_class << " Error: " << Error_string << std::endl;
			}
		}
	}

	MPI_File_close(&file);
	if (mpiWrapper::id() == 0) {
		std::cout << "File loaded \n";
		std::cout.flush();
	}
	(*this).transfer_FieldValues();
}
void SU3_field::loadSU3FromFile(double beta, const std::string& updateMethod, const std::string& dataFolder, const std::string& identifier) {

	std::string beta_str = std::to_string(beta);
	std::replace(beta_str.begin(), beta_str.end(), '.', '_');
	n_fs::path fieldType("/su3/beta" + beta_str + "/");
	int result;
	MPI_File file;
	MPI_Status status;
	MPI_Offset displacement = 0;
	//move to the ensembles directory
	n_fs::path outPath(dataFolder + "ensembles");
	//if not existing, create a directory for the field type
	outPath += fieldType;
	int i;
	for (i = 0; i < m_lattice->getNdims() - 1; i++) {
		outPath += std::to_string(m_lattice->getShape()[i]) + "X";
	}
	outPath += std::to_string(m_lattice->getShape()[i]) + "/" + updateMethod + "/";
	n_fs::create_directories(outPath);
	//write the lattice shape to the filename

	//write the lattice type and the .bin extension to the filename
	outPath += m_lattice->getType() + "_extdof" + std::to_string(m_NrExtDOF) + "_";
	outPath += identifier;
	outPath += ".bin";
	//convert the filesystem path to a format suitable for the MPI_File_open argument
	std::string outPath_string = outPath.string();
	const char* outPath_pointer = outPath_string.c_str();
	if (mpiWrapper::id() == 0) {
		std::cout << "Loading SU3 Field onto " << std::to_string(mpiWrapper::nProcs()) << " processes\nFrom path " << outPath_string << "\n";
		std::cout.flush();
	}

	//Divide the size of FieldArray in bytes by the size of a char
	//then tell MPI that the datatype you're using is char.
	// This way it allocates enough bytes regardless of the size of T
	MPI_Datatype dataType = MPI_CHAR;
	int NrDatatype = 8 * 2 * 9 * m_NrExtDOF * m_lattice->m_responsible_Volume / sizeof(char);

	result = MPI_File_open(mpiWrapper::comm(), outPath_pointer, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	if (result != MPI_SUCCESS) {
		char Error_string[MPI_MAX_ERROR_STRING];
		int Error_class;
		int Error_len;
		MPI_Error_class(result, &Error_class);
		MPI_Error_string(result, Error_string, &Error_len);
		std::cout << Error_class << " Error: " << Error_string << std::endl;
		int term;
		std::cin >> term;
	}
	char datarep[7] = "native";
	result = MPI_File_set_view(file, displacement, dataType, dataType, datarep, MPI_INFO_NULL);
	if (result != MPI_SUCCESS) {
		char Error_string[MPI_MAX_ERROR_STRING];
		int Error_class;
		int Error_len;
		MPI_Error_class(result, &Error_class);
		MPI_Error_string(result, Error_string, &Error_len);
		std::cout << Error_class << " Error: " << Error_string << std::endl;
		int term;
		std::cin >> term;
	}

	//the displacement in bytes sets the offset between the field stored by 
	//process n and process n-1
	//displacement = mpiWrapper::id();
	//displacement *= 8 * 2 * 9 * m_NrExtDOF * m_lattice->m_responsible_Volume;

	//The file is written
	for (int site = (*this).Responsible_Start(); site < (*this).Responsible_Stop(); site++) {
		for (int extDOF = 0; extDOF < m_NrExtDOF; extDOF++) {
			displacement = (m_lattice->m_InternalToTotal_idx[site] * m_NrExtDOF + extDOF) * 144;
			result = MPI_File_read_at(file, displacement, FieldArray[site * m_NrExtDOF + extDOF].getMemPointer(), 144, dataType, &status);
			//result = MPI_File_read_at(file, displacement, FieldArray[site * m_NrExtDOF + extDOF].getMemPointer(), 144, dataType, &status);
			//displacement += 144;
			if (result != MPI_SUCCESS) {
				char Error_string[MPI_MAX_ERROR_STRING];
				int Error_class;
				int Error_len;
				MPI_Error_class(result, &Error_class);
				MPI_Error_string(result, Error_string, &Error_len);
				std::cout << Error_class << " Error: " << Error_string << std::endl;
			}
		}
	}

	MPI_File_close(&file);
	if (mpiWrapper::id() == 0) {
		std::cout << "File loaded \n";
		std::cout.flush();
	}
	(*this).transfer_FieldValues();
}
void SU3_field::transfer_FieldValues(){
	MPI_Request request;
	MPI_Status status;
	int result;
	int Send_id;
	int Recv_id;
	int Nr_sites_toSend;
	int Nr_sites_toRecv;
	bool Sending = false;
	bool Receiving = false;
	C_double* Package;
	C_double* Package_recv;

	for (int i = 1; i < mpiWrapper::nProcs(); i++) {
		Send_id = (mpiWrapper::id() + i) % mpiWrapper::nProcs();
		Recv_id = (mpiWrapper::id() + mpiWrapper::nProcs() - i) % mpiWrapper::nProcs();
		/// <summary>
		/// If the processor <Send_id> shares boundary sites with <id>, Nr_sites_toSend is greater than 0:
		/// Package the field values into an array of C_doubles and send it to the Process in question
		/// </summary>
		Nr_sites_toSend = m_lattice->m_BufferSize[Send_id];

		if (Nr_sites_toSend > 0) {
			Sending = true;
			int Nr_C_doubles_toSend = Nr_sites_toSend * m_NrExtDOF*9;
			Package = new C_double[Nr_C_doubles_toSend];
			for (int i = 0; i < Nr_sites_toSend;i++) {
				for (int mu = 0; mu < m_NrExtDOF; mu++) {
					for (int matEntry = 0; matEntry < 9; matEntry++) {
						Package[i * 9 * m_NrExtDOF + 9 * mu + matEntry] = FieldArray[m_lattice->m_Buffer_receive[Send_id][i] * m_NrExtDOF + mu][matEntry];
					}
				}
			}

			result = MPI_Isend(Package, Nr_C_doubles_toSend * sizeof(C_double) / sizeof(char), MPI_CHAR, Send_id, mpiWrapper::id() * mpiWrapper::nProcs() + Send_id, mpiWrapper::comm(), &request);
			if (result != MPI_SUCCESS) {
				std::cout << "Error: " << result << "-> MPI_Isend	[FIELDVALUES]  From: "<<mpiWrapper::id() << " To: "<<Send_id << "\n";
				char Error_string[MPI_MAX_ERROR_STRING];
				int Error_class;
				int Error_len;
				MPI_Error_class(result, &Error_class);
				MPI_Error_string(result, Error_string, &Error_len);
				std::cout << Error_class << " Error: " << Error_string << std::endl;
			}
		}
		

		Nr_sites_toRecv = m_lattice->m_InternalIdx_stop[Recv_id][1] - m_lattice->m_InternalIdx_start[Recv_id][0];
		if (Nr_sites_toRecv > 0) {
			Receiving = true;
			int Nr_C_doubles_toRecv = Nr_sites_toRecv * m_NrExtDOF * 9;
			Package_recv = new C_double[Nr_C_doubles_toRecv];
			result = MPI_Recv(Package_recv, Nr_C_doubles_toRecv*sizeof(C_double)/sizeof(char),MPI_CHAR,Recv_id, Recv_id * mpiWrapper::nProcs() + mpiWrapper::id(), mpiWrapper::comm(), &status);
			if (result != MPI_SUCCESS) {
				std::cout << "Error: " << result << "-> MPI_Recv	[FIELDVALUES]  From: " << Send_id  << " To: " << mpiWrapper::id() << "\n";
				char Error_string[MPI_MAX_ERROR_STRING];
				int Error_class;
				int Error_len;
				MPI_Error_class(result, &Error_class);
				MPI_Error_string(result, Error_string, &Error_len);
				std::cout << Error_class << " Error: " << Error_string << std::endl;
			}
			int idx = 0;
			for (int intern_idx = m_lattice->m_InternalIdx_start[Recv_id][0]; intern_idx < m_lattice->m_InternalIdx_stop[Recv_id][1];intern_idx++) {
				for (int mu = 0; mu < m_NrExtDOF; mu++) {
					for (int matEntry = 0; matEntry < 9; matEntry++) {
						FieldArray[intern_idx * m_NrExtDOF + mu][matEntry] = Package_recv[idx * 9 * m_NrExtDOF + 9 * mu + matEntry];
					}
				}
				idx++;
			}
		}
		if (Sending) {
			MPI_Status status_waiting;
			MPI_Wait(&request, &status_waiting);
			Sending = false;
			delete[] Package;
			
		}
		if (Receiving) {
			Receiving = false;
			delete[] Package_recv;
		}
		
	}
}

void SU3_field::InitializeHotStart(){
	Random rng;
	su3_mat LinCombGen;
	SU3_gen generators;
	for (int i = (*this).Responsible_Start(); i < (*this).Responsible_Stop(); i++) {
		for (int mu = 0; mu < m_NrExtDOF; mu++) {
			/*for (int j = 0; j < 8; j++) {
				LinCombGen = LinCombGen + rng.Gaussian_Double(0.0, 1.0) * generators(j);
			}
			LinCombGen = LinCombGen.timesI();
			(*this)(i, mu) = HermTrLessExp(LinCombGen);
			LinCombGen.setToZeros();*/
			(*this)(i, mu) = rng.rnd_SU3_Group_elem();

		}
	}
	(*this).transfer_FieldValues();
}
void SU3_field::InitializeColdStart() {
	su3_mat unit;
	unit.setToIdentity();
	for (int i = 0; i < m_lattice->m_thisProc_Volume; i++) {
		for (int mu = 0; mu < m_NrExtDOF; mu++) {
			(*this)(i, mu) = unit;
		}
	}
}

/// <summary>
/// IMPORTANT: This function will not update the boundary field values before calculating!
/// Make sure that the boundary field values are up to date before calling this function.
/// </summary>
double SU3_field::total_PlaquetteSum(){
	double localSum = 0.0;
	double totalSum = 0.0;
	for (int i = Responsible_Start(); i < Responsible_Stop(); i++) {
		for (int mu = 0; mu < m_NrExtDOF-1; mu++) {
			for (int nu = mu + 1; nu < m_NrExtDOF; nu++) {
				localSum += 1.0-1.0/3.0*(plaquette(i, mu, nu).ReTr());
			}
		}
	}
	MPI_Allreduce(&localSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, mpiWrapper::comm());
	return totalSum;
}
double SU3_field::Avg_Plaquette() {
	double localSum = 0.0;
	double totalSum = 0.0;
	int div = m_lattice->m_totalVolume;
	for (int i = Responsible_Start(); i < Responsible_Stop(); i++) {
		for (int mu = 0; mu < m_NrExtDOF - 1; mu++) {
			for (int nu = mu + 1; nu < m_NrExtDOF; nu++) {
				localSum += plaquette(i, mu, nu).ReTr();
			}
		}
	}
	localSum /= 3.0*6.0*(double)div;
	MPI_Allreduce(&localSum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, mpiWrapper::comm());
	return totalSum;
}
/// <summary>
/// Calculates the 
/// </summary>
/// <param name="internal_index"></param>
/// <param name="mu"></param>
/// <returns></returns>
//su3_mat SU3_field::staple(int internal_index, int mu){
//	su3_mat out;
//	int displacedIdx;
//	for (int nu = 0; nu < m_NrExtDOF; nu++) {
//		if (mu != nu) {
//			displacedIdx = m_lattice->m_back[internal_index][nu];
//			out = out + this->fwd_fieldVal(internal_index, mu, nu) * (this->fwd_fieldVal(internal_index, nu, mu)).dagger() * ((*this)(internal_index, nu)).dagger()
//				+ (this->fwd_fieldVal(displacedIdx, mu, nu)).dagger() * ((*this)(displacedIdx, mu)).dagger() * (*this)(displacedIdx, nu);
//		}
//	}
//	return out;
//}

su3_mat SU3_field::staple(int internal_index, int mu){
	su3_mat out;
	int displacedIdx;
	for (int nu = 0; nu < m_NrExtDOF; nu++) {
		if (mu != nu) {
			displacedIdx = m_lattice->m_back[internal_index][nu];
			out = out + this->fwd_fieldVal(internal_index, mu, nu) * ((*this)(internal_index, nu)*this->fwd_fieldVal(internal_index, nu, mu)).dagger()
				+ ((*this)(displacedIdx, mu)*this->fwd_fieldVal(displacedIdx, mu, nu)).dagger() * (*this)(displacedIdx, nu);
		}
	}
	return out;
}


su3_mat SU3_field::clover_avg(int internal_index, int mu, int nu) {
	su3_mat out;
	int displacedIdx;

	//displacedIdx = m_lattice->m_fwd[internal_index][nu];
	//out = (*this)(internal_index, nu) * (*this)(displacedIdx, mu) * this->fwd_fieldVal(internal_index, mu, nu).dagger() * (*this)(internal_index, mu).dagger();
	//
	//displacedIdx = m_lattice->m_back[internal_index][nu];
	//out = out + (*this)(internal_index, mu) * this->fwd_fieldVal(displacedIdx, mu, nu).dagger() * (*this)(displacedIdx, mu).dagger() * (*this)(displacedIdx, nu);
	//
	//displacedIdx = m_lattice->m_back[internal_index][mu];
	//out = out + (*this)(displacedIdx, mu).dagger() * (*this)(displacedIdx, nu) * this->fwd_fieldVal(displacedIdx, nu, mu) * (*this)(internal_index, nu).dagger();

	//displacedIdx = m_lattice->m_back[displacedIdx][nu];
	//out = out + this->back_fieldVal(internal_index, nu, nu).dagger() * (*this)(displacedIdx, mu).dagger() * (*this)(displacedIdx, nu) * this->back_fieldVal(internal_index, mu, mu);

	//By exploiting the property that A.dagger()*B.dagger() = (B*A).dagger(), we can reduce the number of evaluations needed here
	displacedIdx = m_lattice->m_fwd[internal_index][nu];
	out = (*this)(internal_index, nu) * (*this)(displacedIdx, mu) * ((*this)(internal_index, mu) * this->fwd_fieldVal(internal_index, mu, nu)).dagger();

	displacedIdx = m_lattice->m_back[internal_index][nu];
	out = out + (*this)(internal_index, mu) * ((*this)(displacedIdx, mu) * this->fwd_fieldVal(displacedIdx, mu, nu)).dagger() * (*this)(displacedIdx, nu);

	displacedIdx = m_lattice->m_back[internal_index][mu];
	out = out + (*this)(displacedIdx, mu).dagger() * (*this)(displacedIdx, nu) * this->fwd_fieldVal(displacedIdx, nu, mu) * (*this)(internal_index, nu).dagger();

	displacedIdx = m_lattice->m_back[displacedIdx][nu];
	out = out + ((*this)(displacedIdx, mu) * this->back_fieldVal(internal_index, nu, nu)).dagger() * (*this)(displacedIdx, nu) * this->back_fieldVal(internal_index, mu, mu);

	out = 0.25 * out;

	// // RESTORE
	//su3_mat identity;
	//identity.setToIdentity();
	//double trace = 1.0 / 3.0 * out.Tr().I();
	//out = out - out.dagger();
	//out[0].Im -= trace;
	//out[4].Im -= trace;
	//out[8].Im -= trace;
	// // RESTORE

	// // TESTING
	su3_mat identity;
	identity.setToIdentity();
	su3_mat temp;
	temp = out - out.dagger();
	out = (1.0/2.0)*(temp - (1.0 / 3.0) * temp.Tr() * identity).timesMinusI();
	//out = (19.0 / 9.0) * out;// FOR IMPROVED FIELD STRENGTH TENSOR
	// // TESTING

	return out;
}
su3_mat SU3_field::FieldStrengthTensor(int internal_index, int mu, int nu) {
	return clover_avg(internal_index, mu, nu);
}
su3_mat SU3_field::Improved_fieldStrengthTensor(int internal_index, int mu, int nu) {
	//double k1 = 1.7341;
	//double k3 = 0.0251;
	//double k2 = -0.4179;
	double k1 = 19.0/9.0;
	double k2 = 1.0/36.0;
	double k3 = -32.0/45.0;
	double k4 = 1.0 / 15.0;
	return k1 * clover_avg(internal_index, mu, nu) + k2 * clover2x2_avg(internal_index, mu, nu) + k3 * RectangleClover_avg(internal_index, mu, nu) + k4* RectangleClover1x3_avg(internal_index, mu, nu);
}
su3_mat SU3_field::clover2x2_avg(int internal_index, int mu, int nu) {
	su3_mat out;
	int displacedIdx;
	int displacedIdx2;
	int displacedIdx3;

	

	displacedIdx = m_lattice->m_fwd[internal_index][nu];
	displacedIdx2 = m_lattice->m_fwd[displacedIdx][mu];
	displacedIdx3 = m_lattice->m_fwd[internal_index][mu];
	//out = (*this)(internal_index, nu) * (*this)(displacedIdx, nu) * this->fwd_fieldVal(displacedIdx, nu, mu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * this->fwd_fieldVal(displacedIdx2, mu, nu).dagger() * this->fwd_fieldVal(displacedIdx3, mu, nu).dagger() * (*this)(displacedIdx3, mu).dagger()*(*this)(internal_index, mu).dagger();
	out = (*this)(internal_index, nu) * (*this)(displacedIdx, nu) * this->fwd_fieldVal(displacedIdx, nu, mu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * ( (*this)(internal_index, mu) * (*this)(displacedIdx3, mu) * this->fwd_fieldVal(displacedIdx3, mu, nu) * this->fwd_fieldVal(displacedIdx2, mu, nu)).dagger();
	
	
	displacedIdx = m_lattice->m_back[internal_index][nu];
	displacedIdx2 = m_lattice->m_fwd[displacedIdx][mu];
	displacedIdx3 = m_lattice->m_back[displacedIdx2][nu];
	//out = out + (*this)(internal_index, mu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * this->fwd_fieldVal(displacedIdx2, mu, nu).dagger() * this->fwd_fieldVal(displacedIdx3, mu, nu).dagger() * (*this)(displacedIdx3, mu).dagger() * this->back_fieldVal(displacedIdx, nu, mu).dagger() * this->back_fieldVal(displacedIdx, nu, nu) * (*this)(displacedIdx, nu);
	out = out + (*this)(internal_index, mu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * (this->back_fieldVal(displacedIdx, nu, mu) * (*this)(displacedIdx3, mu) * this->fwd_fieldVal(displacedIdx3, mu, nu)  * this->fwd_fieldVal(displacedIdx2, mu, nu)).dagger() * this->back_fieldVal(displacedIdx, nu, nu) * (*this)(displacedIdx, nu);
	
	
	displacedIdx = m_lattice->m_back[internal_index][mu];
	displacedIdx2 = m_lattice->m_fwd[displacedIdx][nu];
	displacedIdx3 = m_lattice->m_back[displacedIdx2][mu];
	//out = out + (*this)(displacedIdx, mu).dagger() * this->back_fieldVal(displacedIdx, mu, mu).dagger() * this->back_fieldVal(displacedIdx, mu, nu) * (*this)(displacedIdx3, nu)	* this->fwd_fieldVal(displacedIdx3, nu, mu) * this->fwd_fieldVal(displacedIdx2, nu, mu)	* this->fwd_fieldVal(displacedIdx2, mu, nu).dagger() * (*this)(internal_index, nu).dagger();
	out = out + (this->back_fieldVal(displacedIdx, mu, mu)* (*this)(displacedIdx, mu)).dagger() * this->back_fieldVal(displacedIdx, mu, nu) * (*this)(displacedIdx3, nu) * this->fwd_fieldVal(displacedIdx3, nu, mu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * ( (*this)(internal_index, nu)* this->fwd_fieldVal(displacedIdx2, mu, nu)).dagger();


	displacedIdx = m_lattice->m_back[displacedIdx][nu];
	displacedIdx2 = m_lattice->m_back[displacedIdx][nu];
	displacedIdx3 = m_lattice->m_back[displacedIdx][mu];
	//out = out + this->fwd_fieldVal(displacedIdx, mu, nu).dagger() * this->fwd_fieldVal(displacedIdx2, mu, nu).dagger() * (*this)(displacedIdx2, mu).dagger() * this->back_fieldVal(displacedIdx3, nu, mu).dagger() * this->back_fieldVal(displacedIdx3, nu, nu) * (*this)(displacedIdx3, nu) * this->fwd_fieldVal(displacedIdx3, nu, mu) * this->fwd_fieldVal(displacedIdx, nu, mu);
	out = out + (this->back_fieldVal(displacedIdx3, nu, mu) * (*this)(displacedIdx2, mu) * this->fwd_fieldVal(displacedIdx2, mu, nu) * this->fwd_fieldVal(displacedIdx, mu, nu)).dagger() * this->back_fieldVal(displacedIdx3, nu, nu) * (*this)(displacedIdx3, nu) * this->fwd_fieldVal(displacedIdx3, nu, mu) * this->fwd_fieldVal(displacedIdx, nu, mu);
	
	
	out = 0.25 * out;

	// // TESTING
	su3_mat identity;
	identity.setToIdentity();
	su3_mat temp;
	temp = out - out.dagger();
	out = (1.0 / 2.0) * (temp - (1.0 / 3.0) * temp.Tr() * identity).timesMinusI();
	//out = (1.0 / 36.0) * out;// FOR IMPROVED FIELD STRENGTH TENSOR
	// // TESTING

	return out;
}
su3_mat SU3_field::RectangleClover_avg(int internal_index, int mu, int nu) {
	su3_mat tall_rectangle;
	su3_mat wide_rectangle;
	su3_mat out;
	int displacedIdx;
	int displacedIdx2;

	displacedIdx = m_lattice->m_fwd[internal_index][nu];
	//tall_rectangle = (*this)(internal_index, nu) * (*this)(displacedIdx, nu) * this->fwd_fieldVal(displacedIdx, nu,mu) * this->fwd_fieldVal(displacedIdx, mu, nu).dagger() * this->fwd_fieldVal(internal_index, mu, nu).dagger()* (*this)(internal_index, mu).dagger();
	tall_rectangle = (*this)(internal_index, nu) * (*this)(displacedIdx, nu) * this->fwd_fieldVal(displacedIdx, nu, mu) * ((*this)(internal_index, mu)* this->fwd_fieldVal(internal_index, mu, nu) * this->fwd_fieldVal(displacedIdx, mu, nu)).dagger();
	
	displacedIdx = m_lattice->m_back[internal_index][nu];
	displacedIdx2 = m_lattice->m_back[displacedIdx][nu];
	//tall_rectangle = tall_rectangle + (*this)(internal_index, mu) * this->fwd_fieldVal(displacedIdx, mu, nu).dagger() * this->fwd_fieldVal(displacedIdx2, mu, nu).dagger() * (*this)(displacedIdx2, mu).dagger() * (*this)(displacedIdx2, nu) * (*this)(displacedIdx, nu);
	tall_rectangle = tall_rectangle + (*this)(internal_index, mu) * ((*this)(displacedIdx2, mu) * this->fwd_fieldVal(displacedIdx2, mu, nu) * this->fwd_fieldVal(displacedIdx, mu, nu)).dagger() * (*this)(displacedIdx2, nu) * (*this)(displacedIdx, nu);
	
	displacedIdx = m_lattice->m_back[internal_index][mu];
	displacedIdx2 = m_lattice->m_fwd[displacedIdx][nu];
	//tall_rectangle = tall_rectangle + (*this)(displacedIdx, mu).dagger() * (*this)(displacedIdx, nu) * (*this)(displacedIdx2, nu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * this->fwd_fieldVal(displacedIdx2, mu, nu).dagger() * (*this)(internal_index, nu).dagger();
	tall_rectangle = tall_rectangle + (*this)(displacedIdx, mu).dagger() * (*this)(displacedIdx, nu) * (*this)(displacedIdx2, nu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * ((*this)(internal_index, nu) *this->fwd_fieldVal(displacedIdx2, mu, nu) ).dagger();

	displacedIdx = m_lattice->m_back[displacedIdx][nu];
	displacedIdx2 = m_lattice->m_back[displacedIdx][nu];
	//tall_rectangle = tall_rectangle + this->back_fieldVal(internal_index, nu, nu).dagger() * this->fwd_fieldVal(displacedIdx2, mu, nu).dagger() * (*this)(displacedIdx2, mu).dagger() * (*this)(displacedIdx2, nu) * (*this)(displacedIdx, nu) * this->back_fieldVal(internal_index, mu, mu);
	tall_rectangle = tall_rectangle + ((*this)(displacedIdx2, mu) * this->fwd_fieldVal(displacedIdx2, mu, nu) * this->back_fieldVal(internal_index, nu, nu)).dagger() * (*this)(displacedIdx2, nu) * (*this)(displacedIdx, nu) * this->back_fieldVal(internal_index, mu, mu);
	tall_rectangle = (1.0/8.0) * tall_rectangle;




	displacedIdx = m_lattice->m_back[internal_index][nu];
	displacedIdx2 = m_lattice->m_fwd[displacedIdx][mu];
	//wide_rectangle = (*this)(internal_index, mu) * this->fwd_fieldVal(internal_index, mu, mu) * this->fwd_fieldVal(displacedIdx2, mu, nu).dagger() * (*this)(displacedIdx2, mu).dagger() *(*this)(displacedIdx, mu).dagger() * (*this)(displacedIdx, nu);
	wide_rectangle = (*this)(internal_index, mu) * this->fwd_fieldVal(internal_index, mu, mu) * ( (*this)(displacedIdx, mu) * (*this)(displacedIdx2, mu) * this->fwd_fieldVal(displacedIdx2, mu, nu)).dagger() * (*this)(displacedIdx, nu);

	displacedIdx = m_lattice->m_back[internal_index][mu];
	displacedIdx2 = m_lattice->m_back[displacedIdx][mu];
	//wide_rectangle = wide_rectangle + this->back_fieldVal(internal_index, nu, nu).dagger() * this->back_fieldVal(displacedIdx, nu, mu).dagger() * this->back_fieldVal(displacedIdx2, nu, mu).dagger() * this->back_fieldVal(displacedIdx2, nu, nu) * (*this)(displacedIdx2, mu) * (*this)(displacedIdx, mu);
	wide_rectangle = wide_rectangle + (this->back_fieldVal(displacedIdx2, nu, mu) * this->back_fieldVal(displacedIdx, nu, mu) * this->back_fieldVal(internal_index, nu, nu)).dagger() * this->back_fieldVal(displacedIdx2, nu, nu) * (*this)(displacedIdx2, mu) * (*this)(displacedIdx, mu);

	displacedIdx = m_lattice->m_back[internal_index][mu];
	displacedIdx2 = m_lattice->m_back[displacedIdx][mu];
	//wide_rectangle = wide_rectangle + (*this)(displacedIdx, mu).dagger() * (*this)(displacedIdx2, mu).dagger() * (*this)(displacedIdx2, nu) * this->fwd_fieldVal(displacedIdx2, nu, mu)* this->fwd_fieldVal(displacedIdx, nu, mu) * (*this)(internal_index, nu).dagger();
	wide_rectangle = wide_rectangle + ((*this)(displacedIdx2, mu)*(*this)(displacedIdx, mu)).dagger() * (*this)(displacedIdx2, nu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * this->fwd_fieldVal(displacedIdx, nu, mu) * (*this)(internal_index, nu).dagger();

	displacedIdx = m_lattice->m_fwd[internal_index][mu];
	//wide_rectangle = wide_rectangle + (*this)(internal_index, nu) * this->fwd_fieldVal(internal_index, nu, mu) * this->fwd_fieldVal(displacedIdx, nu, mu) * this->fwd_fieldVal(displacedIdx, mu, nu).dagger() * (*this)(displacedIdx, mu).dagger() * (*this)(internal_index, mu).dagger();
	wide_rectangle = wide_rectangle + (*this)(internal_index, nu) * this->fwd_fieldVal(internal_index, nu, mu) * this->fwd_fieldVal(displacedIdx, nu, mu) * ((*this)(internal_index, mu) * (*this)(displacedIdx, mu) * this->fwd_fieldVal(displacedIdx, mu, nu)).dagger();

	wide_rectangle = (1.0 / 8.0) * wide_rectangle;

	

	out = tall_rectangle + wide_rectangle;

	// // TESTING
	su3_mat identity;
	identity.setToIdentity();
	su3_mat temp;
	temp = out - out.dagger();
	out = (1.0 / 2.0) * (temp - (1.0 / 3.0) * temp.Tr() * identity).timesMinusI();
	//out = -(32.0 / 45.0) * out;
	// // TESTING
	return out;
}

su3_mat SU3_field::RectangleClover1x3_avg(int internal_index, int mu, int nu) {
	su3_mat tall_rectangle;
	su3_mat wide_rectangle;
	su3_mat out;
	int displacedIdx;
	int displacedIdx2;
	int displacedIdx3;

	

	displacedIdx = m_lattice->m_fwd[internal_index][nu];
	displacedIdx2 = m_lattice->m_fwd[displacedIdx][nu];
	//tall_rectangle = (*this)(internal_index, nu) * (*this)(displacedIdx, nu)* (*this)(displacedIdx2, nu) * this->fwd_fieldVal(displacedIdx2, nu, mu) *	this->fwd_fieldVal(displacedIdx2, mu, nu).dagger() * this->fwd_fieldVal(displacedIdx, mu, nu).dagger() *  this->fwd_fieldVal(internal_index, mu, nu).dagger() * (*this)(internal_index, mu).dagger();
	tall_rectangle = (*this)(internal_index, nu) * (*this)(displacedIdx, nu) * (*this)(displacedIdx2, nu) * this->fwd_fieldVal(displacedIdx2, nu, mu)  * ((*this)(internal_index, mu)* this->fwd_fieldVal(internal_index, mu, nu)* this->fwd_fieldVal(displacedIdx, mu, nu)* this->fwd_fieldVal(displacedIdx2, mu, nu)).dagger();

	displacedIdx = m_lattice->m_back[internal_index][nu];
	displacedIdx2 = m_lattice->m_back[displacedIdx][nu];
	displacedIdx3 = m_lattice->m_back[displacedIdx2][nu];
	//tall_rectangle = tall_rectangle + (*this)(internal_index, mu) * this->fwd_fieldVal(displacedIdx, mu, nu).dagger() * this->fwd_fieldVal(displacedIdx2, mu, nu).dagger() * this->fwd_fieldVal(displacedIdx3, mu, nu).dagger()	* (*this)(displacedIdx3, mu).dagger() * (*this)(displacedIdx3, nu) * (*this)(displacedIdx2, nu) * (*this)(displacedIdx, nu);
	tall_rectangle = tall_rectangle + (*this)(internal_index, mu) * ((*this)(displacedIdx3, mu) * this->fwd_fieldVal(displacedIdx3, mu, nu) * this->fwd_fieldVal(displacedIdx2, mu, nu)* this->fwd_fieldVal(displacedIdx, mu, nu)).dagger() * (*this)(displacedIdx3, nu) * (*this)(displacedIdx2, nu) * (*this)(displacedIdx, nu);
	

	displacedIdx = m_lattice->m_fwd[internal_index][nu];
	displacedIdx2 = m_lattice->m_fwd[displacedIdx][nu];
	displacedIdx3 = m_lattice->m_fwd[displacedIdx2][nu];
	//tall_rectangle = tall_rectangle + this->back_fieldVal(internal_index, mu, mu).dagger() * this->back_fieldVal(internal_index, mu, nu) * this->back_fieldVal(displacedIdx, mu, nu) * this->back_fieldVal(displacedIdx2, mu, nu) * this->back_fieldVal(displacedIdx3, mu, mu) * (*this)(displacedIdx2, nu).dagger() * (*this)(displacedIdx, nu).dagger() * (*this)(internal_index, nu).dagger();
	tall_rectangle = tall_rectangle + this->back_fieldVal(internal_index, mu, mu).dagger() * this->back_fieldVal(internal_index, mu, nu) * this->back_fieldVal(displacedIdx, mu, nu) * this->back_fieldVal(displacedIdx2, mu, nu) * this->back_fieldVal(displacedIdx3, mu, mu) * ((*this)(internal_index, nu) * (*this)(displacedIdx, nu) * (*this)(displacedIdx2, nu)).dagger();
	

	displacedIdx = m_lattice->m_back[internal_index][nu];
	displacedIdx2 = m_lattice->m_back[displacedIdx][nu];
	displacedIdx3 = m_lattice->m_back[displacedIdx2][nu];
	//tall_rectangle = tall_rectangle + (*this)(displacedIdx, nu).dagger() * (*this)(displacedIdx2, nu).dagger() * (*this)(displacedIdx3, nu).dagger() * this->back_fieldVal(displacedIdx3, mu, mu).dagger() * this->back_fieldVal(displacedIdx3, mu, nu) * this->back_fieldVal(displacedIdx2, mu, nu) * this->back_fieldVal(displacedIdx, mu, nu) * this->back_fieldVal(internal_index, mu, mu);
	tall_rectangle = tall_rectangle + (this->back_fieldVal(displacedIdx3, mu, mu) * (*this)(displacedIdx3, nu) * (*this)(displacedIdx2, nu)* (*this)(displacedIdx, nu)).dagger() * this->back_fieldVal(displacedIdx3, mu, nu) * this->back_fieldVal(displacedIdx2, mu, nu) * this->back_fieldVal(displacedIdx, mu, nu) * this->back_fieldVal(internal_index, mu, mu);
	
	tall_rectangle = (1.0 / 8.0) * tall_rectangle;

	


	displacedIdx = m_lattice->m_fwd[internal_index][mu];
	displacedIdx2 = m_lattice->m_fwd[displacedIdx][mu];
	displacedIdx3 = m_lattice->m_fwd[displacedIdx2][mu];
	//wide_rectangle = (*this)(internal_index, nu) * this->fwd_fieldVal(internal_index, nu, mu) * this->fwd_fieldVal(displacedIdx, nu, mu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * (*this)(displacedIdx3, nu).dagger() * (*this)(displacedIdx2, mu).dagger() * (*this)(displacedIdx, mu).dagger() * (*this)(internal_index, mu).dagger();
	wide_rectangle = (*this)(internal_index, nu) * this->fwd_fieldVal(internal_index, nu, mu) * this->fwd_fieldVal(displacedIdx, nu, mu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * ((*this)(internal_index, mu) * (*this)(displacedIdx, mu) * (*this)(displacedIdx2, mu) * (*this)(displacedIdx3, nu)).dagger();

	displacedIdx = m_lattice->m_fwd[internal_index][mu];
	displacedIdx2 = m_lattice->m_fwd[displacedIdx][mu];
	displacedIdx3 = m_lattice->m_fwd[displacedIdx2][mu];
	//wide_rectangle = wide_rectangle + (*this)(internal_index, mu) * (*this)(displacedIdx, mu) * (*this)(displacedIdx2, mu) * this->back_fieldVal(displacedIdx3, nu, nu).dagger() * this->back_fieldVal(displacedIdx2, nu, mu).dagger() * this->back_fieldVal(displacedIdx, nu, mu).dagger() * this->back_fieldVal(internal_index, nu, mu).dagger() * this->back_fieldVal(internal_index, nu, nu);
	wide_rectangle = wide_rectangle + (*this)(internal_index, mu) * (*this)(displacedIdx, mu) * (*this)(displacedIdx2, mu) * (this->back_fieldVal(internal_index, nu, mu) * this->back_fieldVal(displacedIdx, nu, mu) * this->back_fieldVal(displacedIdx2, nu, mu) * this->back_fieldVal(displacedIdx3, nu, nu)).dagger() * this->back_fieldVal(internal_index, nu, nu);

	displacedIdx = m_lattice->m_back[internal_index][mu];
	displacedIdx2 = m_lattice->m_back[displacedIdx][mu];
	displacedIdx3 = m_lattice->m_back[displacedIdx2][mu];
	//wide_rectangle = wide_rectangle + (*this)(displacedIdx, mu).dagger() * (*this)(displacedIdx2, mu).dagger() * (*this)(displacedIdx3, mu).dagger() * (*this)(displacedIdx3, nu) * this->fwd_fieldVal(displacedIdx3, nu, mu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * this->fwd_fieldVal(displacedIdx, nu, mu) * (*this)(internal_index, nu).dagger();
	wide_rectangle = wide_rectangle + ((*this)(displacedIdx3, mu) * (*this)(displacedIdx2, mu)* (*this)(displacedIdx, mu)).dagger() * (*this)(displacedIdx3, nu) * this->fwd_fieldVal(displacedIdx3, nu, mu) * this->fwd_fieldVal(displacedIdx2, nu, mu) * this->fwd_fieldVal(displacedIdx, nu, mu) * (*this)(internal_index, nu).dagger();

	displacedIdx = m_lattice->m_back[internal_index][mu];
	displacedIdx2 = m_lattice->m_back[displacedIdx][mu];
	displacedIdx3 = m_lattice->m_back[displacedIdx2][mu];
	//wide_rectangle = wide_rectangle + this->back_fieldVal(internal_index, nu, nu).dagger() * this->back_fieldVal(displacedIdx, nu, mu).dagger() * this->back_fieldVal(displacedIdx2, nu, mu).dagger() * this->back_fieldVal(displacedIdx3, nu, mu).dagger() * this->back_fieldVal(displacedIdx3, nu, nu) * (*this)(displacedIdx3, mu) * (*this)(displacedIdx2, mu) * (*this)(displacedIdx, mu);
	wide_rectangle = wide_rectangle + (this->back_fieldVal(displacedIdx3, nu, mu) * this->back_fieldVal(displacedIdx2, nu, mu) * this->back_fieldVal(displacedIdx, nu, mu)*this->back_fieldVal(internal_index, nu, nu)).dagger() * this->back_fieldVal(displacedIdx3, nu, nu) * (*this)(displacedIdx3, mu) * (*this)(displacedIdx2, mu) * (*this)(displacedIdx, mu);
	

	wide_rectangle = (1.0 / 8.0) * wide_rectangle;


	out = tall_rectangle + wide_rectangle;

	// // TESTING
	su3_mat identity;
	identity.setToIdentity();
	su3_mat temp;
	temp = out - out.dagger();
	out = (1.0 / 2.0) * (temp - (1.0 / 3.0) * temp.Tr() * identity).timesMinusI();
	//out = -(32.0 / 45.0) * out;
	// // TESTING
	return out;
}

inline su3_mat SU3_field::plaquette(int internal_index, int mu, int nu){
	su3_mat out;
	out = (*this)(internal_index, mu) * this->fwd_fieldVal(internal_index, mu, nu) * ((*this)(internal_index, nu)*this->fwd_fieldVal(internal_index, nu, mu)).dagger();
	return out;
}

su3_mat SU3_field::Plaquette_FieldStrengthTensor(int internal_index, int mu, int nu) {
	su3_mat temp;
	su3_mat identity;
	identity.setToIdentity();
	temp = (*this).plaquette(internal_index, mu, nu) - (*this).plaquette(internal_index, mu, nu).dagger();
	return (1.0 / 2.0) * (temp - (1.0 / 3.0) * temp.Tr() * identity).timesMinusI();
}

void SU3_field::operator=(const SU3_field& field){
	for (int i = 0; i < field.m_NrExtDOF*field.m_lattice->getthisProc_Volume(); i++) {
		FieldArray[i] = field.FieldArray[i];
	}
}


