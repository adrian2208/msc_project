#pragma once
#include "Field.h"
#include "../Math/Random.h"
#pragma once
class RNG_field : public Field<Random> {
public:
	RNG_field(Lattice& lattice);

};