#pragma once
#include "../Field/SU3_field.h"
#include "../Action/Action.h"
#include "../Math/SU3_mat.h"
#include "../Observables/TopologicalCharge.h"
class GradientFlow {
public:
	GradientFlow(Action& GaugeAction, SU3_field& U, double epsilon);
	void flow();

	void Include_TopCharge(TopologicalCharge& TopCharge);

	void MakeMeasurements();

private: 
	Action* m_GaugeAction;
	SU3_field* m_U;
	SU3_field* m_Z0;
	SU3_field* m_Z1;
	SU3_field* m_Z2;

	double m_epsilon;
	double m_flowTime;

	TopologicalCharge* m_topologicalCharge = nullptr;
};