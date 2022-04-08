#pragma once
#include "../Field/SU3_field.h"
#include "../Action/Action.h"
#include "../Math/SU3_mat.h"
#include "../Observables/TopologicalCharge.h"
#include "../Observables/EnergyDensity.h"
class GradientFlow {
public:
	GradientFlow(Action& GaugeAction, SU3_field& U, double epsilon, int FlowStepsPerMeasurement=1);
	~GradientFlow();
	void flow();

	void Include_TopCharge(TopologicalCharge& TopCharge);
	void Include_EnergyDensity(EnergyDensity& Edensity);
	
	void MakeMeasurements();
	std::string getupdateMethod() const;
	double GetFlowTime() const;
private: 
	std::string m_updateMethod;
	Action* m_GaugeAction;
	SU3_field* m_U;
	SU3_field* m_Z0;
	SU3_field* m_Z1;
	SU3_field* m_Z2;
	int m_NrSteps;
	int m_FlowStepsPerMeasurement;
	double m_epsilon;
	double m_flowTime;

	TopologicalCharge* m_topologicalCharge = nullptr;
	EnergyDensity* m_EnergyDensity = nullptr;
};