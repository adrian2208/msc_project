#include "Field.h"
#include "../Math/SU3_mat.h"
class SU3_field : public Field<su3_mat> {
public:
	SU3_field(Lattice& lattice, int NrExtDOF);//NrExtDOF corresponds to e.g. mu={0,1,2,3}
	void saveSU3ToFile(const std::filesystem::path& identifier = "");
	void loadSU3FromFile(const std::filesystem::path& identifier = "");

};