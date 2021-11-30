#include "Field.h"
class SU3_field : public Field<int> {
public:
	SU3_field(Lattice& lattice, int NrExtDOF);//NrExtDOF corresponds to e.g. mu={0,1,2,3}

};