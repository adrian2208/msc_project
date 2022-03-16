#pragma once
#include "../Field/Field.h"
#include "../Field/SU3_field.h"
class Action {
public:
	Action();
	virtual double calculate_Action(SU3_field& field) = 0;
	virtual void calculate_Force(SU3_field& field, SU3_field& F) = 0;
	virtual void calculate_FlowGradient(SU3_field& field, SU3_field& F) = 0;
	virtual void calculate_FlowGradient(SU3_field& field, SU3_field& F, int parity, int mu) = 0;
};