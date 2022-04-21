#include "Double_field.h"
#include <algorithm>
Double_field::Double_field(Lattice& lattice, int NrExtDOF) : Field(lattice, NrExtDOF) {
	for (int i = (*this).Responsible_Start(); i < (*this).Responsible_Stop(); i++) {
		for (int mu = 0; mu < m_NrExtDOF; mu++) {
			(*this)(i, mu) = 0.0;
		}
	}
}

void Double_field::operator=(const Double_field& field) {
	for (int i = 0; i < field.m_NrExtDOF * field.m_lattice->getthisProc_Volume(); i++) {
		FieldArray[i] = field.FieldArray[i];
	}
}

void Double_field::saveDoubleToFile(double beta, const std::string& updateMethod, const std::string& dataFolder, const std::string& identifier) {

	std::string beta_str = std::to_string(beta);
	std::replace(beta_str.begin(), beta_str.end(), '.', '_');
	n_fs::path fieldType("/Qdensity/beta" + beta_str + "/");
	int result;
	MPI_File file;
	MPI_Status status;
	MPI_Offset displacement = 0;
	//move to the ensembles directory
	n_fs::path outPath(dataFolder + "Observables");
	//if not existing, create a directory for the field type
	outPath += fieldType;
	int i;
	for (i = 0; i < (*this).getLatticePtr().getNdims() - 1; i++) {
		outPath += std::to_string((*this).getLatticePtr().getShape()[i]) + "X";
	}
	outPath += std::to_string((*this).getLatticePtr().getShape()[i]) + "/" + updateMethod + "/";
	n_fs::create_directories(outPath);
	//write the lattice shape to the filename

	//write the lattice type and the .bin extension to the filename
	outPath += (*this).getLatticePtr().getType() + "_extdof" + std::to_string((*this).getNrExtDOF()) + "_";
	outPath += identifier;
	outPath += ".bin";
	//convert the filesystem path to a format suitable for the MPI_File_open argument
	std::string outPath_string = outPath.string();
	const char* outPath_pointer = outPath_string.c_str();

	//Divide the size of FieldArray in bytes by the size of a char
	//then tell MPI that the datatype you're using is char.
	// This way it allocates enough bytes regardless of the size of T
	MPI_Datatype dataType = MPI_CHAR;
	int NrDatatype = 8 * m_NrExtDOF * m_lattice->m_responsible_Volume / sizeof(char);

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
		std::cout << "Saving Double Field to path " << outPath_string << "\n";
		std::cout.flush();
	}
	//The file is written
	for (int site = (*this).Responsible_Start(); site < (*this).Responsible_Stop(); site++) {
		for (int extDOF = 0; extDOF < m_NrExtDOF; extDOF++) {
			//std::cout << "id: " << mpiWrapper::id() << " location: " << m_lattice->m_InternalToTotal_idx[site] << std::endl;
			displacement = (m_lattice->m_InternalToTotal_idx[site] * m_NrExtDOF + extDOF) * 8;
			result = MPI_File_write_at(file, displacement, &FieldArray[site * m_NrExtDOF + extDOF], 8, dataType, &status);
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