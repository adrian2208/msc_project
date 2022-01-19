#include "SU3_field.h"

SU3_field::SU3_field(Lattice& lattice, int NrExtDOF) : Field(lattice, NrExtDOF){
	m_FieldArray_NrBytes = 8*2*9 * NrExtDOF * lattice.m_thisProc_Volume;//9 entries of C_double each containing 16 bytes
}

void SU3_field::saveSU3ToFile(const std::filesystem::path& identifier){
	std::filesystem::path fieldType("/su3/");
	int result;
	MPI_File file;
	MPI_Status status;
	MPI_Offset displacement = 0;
	//move to the ensembles directory
	std::filesystem::path outPath("../../../../data/ensembles");
	//if not existing, create a directory for the field type
	outPath += fieldType;
	std::filesystem::create_directories(outPath);
	//write the lattice shape to the filename
	int i;
	for (i = 0; i < m_lattice->getNdims() - 1; i++) {
		outPath += std::to_string(m_lattice->getShape()[i]) + "X";
	}
	//write the lattice type and the .bin extension to the filename
	outPath += std::to_string(m_lattice->getShape()[i]) + "_" + m_lattice->getType() + "_extdof" + std::to_string(m_NrExtDOF) + "_" + std::to_string(mpiWrapper::nProcs()) + "_";
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
		std::cout << "Error: " << result << "-> MPI_File_open\n";
	}
	result = MPI_File_set_view(file, displacement, dataType, dataType, (char*)NULL, MPI_INFO_NULL);
	if (result != MPI_SUCCESS) {
		std::cout << "Error: " << result << "-> MPI_File_set_view\n";
	}

	//the displacement in bytes sets the offset between the field stored by 
	//process n and process n-1
	displacement = mpiWrapper::id();
	displacement *= m_FieldArray_NrBytes;

	//The file is written
	for (int site = 0; site < m_lattice->getthisProc_Volume(); site++) {
		for (int extDOF = 0; extDOF < m_NrExtDOF; extDOF++) {
			result = MPI_File_write_at(file, displacement, FieldArray[site*m_NrExtDOF+extDOF].getMemPointer(), 144, dataType, &status);
			displacement += 144;
		}
	}
	if (result != MPI_SUCCESS) {
		std::cout << "Error: " << result << "-> MPI_File_write_at\n";
	}
	MPI_File_close(&file);
}

void SU3_field::loadSU3FromFile(const std::filesystem::path& identifier){
	std::filesystem::path fieldType("/su3/");
	int result;
	MPI_File file;
	MPI_Status status;
	MPI_Offset displacement = 0;
	//move to the ensembles directory
	std::filesystem::path outPath("../../../../data/ensembles");
	//if not existing, create a directory for the field type
	outPath += fieldType;
	std::filesystem::create_directories(outPath);
	//write the lattice shape to the filename
	int i;
	for (i = 0; i < m_lattice->getNdims() - 1; i++) {
		outPath += std::to_string(m_lattice->getShape()[i]) + "X";
	}
	//write the lattice type and the .bin extension to the filename
	outPath += std::to_string(m_lattice->getShape()[i]) + "_" + m_lattice->getType() + "_extdof" + std::to_string(m_NrExtDOF) + "_" + std::to_string(mpiWrapper::nProcs()) + "_";
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

	result = MPI_File_open(mpiWrapper::comm(), outPath_pointer, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	if (result != MPI_SUCCESS) {
		std::cout << "Error: " << result << "-> MPI_File_open\n";
	}
	result = MPI_File_set_view(file, displacement, dataType, dataType, (char*)NULL, MPI_INFO_NULL);
	if (result != MPI_SUCCESS) {
		std::cout << "Error: " << result << "-> MPI_File_set_view\n";
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
		}
	}
	if (result != MPI_SUCCESS) {
		std::cout << "Error: " << result << "-> MPI_File_write_at\n";
	}
	MPI_File_close(&file);
}

void SU3_field::CommunicateFieldValues(){
}

double SU3_field::total_PlaquetteSum(){
	double localSum = 0.0;
	double totalSum = 0.0;
	for (int i = 0; i < m_lattice->m_responsible_Volume; i++) {
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
	int div = m_lattice->m_responsible_Volume;
	for (int i = 0; i < m_lattice->m_responsible_Volume; i++) {
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
su3_mat SU3_field::staple(int internal_index, int mu){
	su3_mat out;
	//NOT NESCESSARY -- SEE DEFAULT CONSTRUCTOR!
	out[0] = C_double(0.0, 0.0);
	out[1] = C_double(0.0, 0.0);
	out[2] = C_double(0.0, 0.0);
	out[3] = C_double(0.0, 0.0);
	out[4] = C_double(0.0, 0.0);
	out[5] = C_double(0.0, 0.0);
	out[6] = C_double(0.0, 0.0);
	out[7] = C_double(0.0, 0.0);
	out[8] = C_double(0.0, 0.0);
	//-----DELETE ME ----------------
	int displacedIdx;
	for (int nu = 0; nu < m_NrExtDOF; nu++) {
		if (mu != nu) {
			displacedIdx = m_lattice->m_back[internal_index][nu];
			out = out + this->fwd_fieldVal(internal_index, mu, nu) * (this->fwd_fieldVal(internal_index, nu, mu)).dagger() * ((*this)(internal_index, nu)).dagger()
				+ (this->fwd_fieldVal(displacedIdx, mu, nu)).dagger() * ((*this)(displacedIdx, mu)).dagger() * (*this)(displacedIdx, nu);
		}
	}
	return out;
}

inline su3_mat SU3_field::plaquette(int internal_index, int mu, int nu){
	su3_mat out;
	out = (*this)(internal_index, mu) * this->fwd_fieldVal(internal_index, mu, nu) * this->fwd_fieldVal(internal_index, nu, mu).dagger() * (*this)(internal_index, nu).dagger();
	isSpecialUnitary(out);
	return out;
}

void SU3_field::operator=(const SU3_field& field){
	for (int i = 0; i < field.m_NrExtDOF*field.m_lattice->getthisProc_Volume(); i++) {
		FieldArray[i] = field.FieldArray[i];
	}
}

//void  SU3_field::copyFieldVals(const SU3_field& field) {
//	for (int i = 0; i < field.m_NrExtDOF * field.m_lattice->getthisProc_Volume(); i++) {
//				FieldArray[i] = field.FieldArray[i];
//		}
//}