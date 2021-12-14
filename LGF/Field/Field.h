#pragma once
#include "../Lattice/Lattice.h"
#include "../Vertex.h"
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
	void saveToFile() {
		MPI_File file;
		MPI_Status status;
		const char filename[128] = "../../../../data/ensembles/trial.bin";
		MPI_Offset displacement=0;
		int result;
		//Divide the size of FieldArray in bytes by the size of a char
		//then tell MPI that the datatype you're using is char.
		// This way it allocates enough bytes regardless of the size of T
		MPI_Datatype dataType = MPI_CHAR;
		int NrDatatype = m_FieldArray_NrBytes / sizeof(char);

		result = MPI_File_open(mpiWrapper::comm(), filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
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
	void loadFromFile() {};

protected:
	T* FieldArray;
	int m_NrExtDOF;//the number of external DOF e.g. 4: mu=0,1,2,3
private:
	Lattice* m_lattice;
	int m_FieldArray_NrBytes;
};