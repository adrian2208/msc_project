#pragma once
#include "../Lattice/Lattice.h"
#include "../Vertex.h"
#include <filesystem>
#include <string>

template<class T>
class Field {
public:
	Field(Lattice& lattice, int NrExtDOF) {
		m_lattice = &lattice;
		m_NrExtDOF = NrExtDOF;
		int NrTs = NrExtDOF * lattice.m_thisProc_Volume;
		FieldArray = new T[NrTs];
		m_FieldArray_NrBytes = sizeof(T)* NrExtDOF * lattice.m_thisProc_Volume;
	}
	
	// overloads () operator on Field instance to return T situated on Vertex x and ext. DOF mu by reference
	inline T& operator() (Vertex x, int mu) {
		return FieldArray[x.index * m_NrExtDOF + mu];
	}
	inline T& operator() (int i, int mu) {
		return FieldArray[i * m_NrExtDOF + mu];
	}
	
	//Writes a binary file containing the FieldArray of each process
	//ordered by their rank. Contained in data- directory
	void saveToFile(const std::filesystem::path& fieldType = "/generic/", const std::filesystem::path& identifier = "") {
		int result;
		MPI_File file;
		MPI_Status status;
		MPI_Offset displacement = 0;
		//move to the ensembles directory
		std::filesystem::path outPath("../../../../data/ensembles");
		//if not existing, create a directory for the field type
		outPath += fieldType;
		std::filesystem::create_directory(outPath);
		//write the lattice shape to the filename
		int i;
		for (i = 0; i < m_lattice->getNdims()-1; i++) {
			outPath += std::to_string(m_lattice->getShape()[i])+"X";
		}
		//write the lattice type and the .bin extension to the filename
		outPath += std::to_string(m_lattice->getShape()[i]) + "_" + m_lattice->getType() + "_" + std::to_string(mpiWrapper::nProcs()) + "_";
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
		result = MPI_File_write_at(file, displacement, FieldArray, NrDatatype, dataType, &status);
		if (result != MPI_SUCCESS) {
			std::cout << "Error: " << result << "-> MPI_File_write_at\n";
		}
		MPI_File_close(&file);
	}


	//Field configuration must only be loaded To a program using the 
	//same process allocation used when that file was created.
	//The file is not configured by spatial position, but by
	//the rank of the process that contained the field
	void loadFromFile(const std::filesystem::path& fieldType = "/generic/", const std::filesystem::path& identifier = "") {
		int result;
		MPI_File file;
		MPI_Status status;
		MPI_Offset displacement = 0;
		//move to the ensembles directory
		std::filesystem::path outPath("../../../../data/ensembles");
		//if not existing, create a directory for the field type
		outPath += fieldType;
		std::filesystem::create_directory(outPath);
		//write the lattice shape to the filename
		int i;
		for (i = 0; i < m_lattice->getNdims() - 1; i++) {
			outPath += std::to_string(m_lattice->getShape()[i]) + "X";
		}
		//write the lattice type and the .bin extension to the filename
		outPath += std::to_string(m_lattice->getShape()[i]) + "_" + m_lattice->getType() + "_" + std::to_string(mpiWrapper::nProcs()) + "_";
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
		result = MPI_File_read_at(file, displacement, FieldArray, NrDatatype, dataType, &status);
		if (result != MPI_SUCCESS) {
			std::cout << "Error: " << result << "-> MPI_File_read_at\n";
		}
		MPI_File_close(&file);

	};

protected:
	T* FieldArray;
	int m_NrExtDOF;//the number of external DOF e.g. 4: mu=0,1,2,3
private:
	Lattice* m_lattice;
	int m_FieldArray_NrBytes;
};