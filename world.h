#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

#include <time.h>

// include elastic rod class
#include "elasticRod.h"

// include force classes
#include "elasticStretchingForce.h"
#include "elasticBendingForce.h"
#include "elasticTwistingForce.h"
#include "externalGravityForce.h"
#include "inertialForce.h"

// include external force
#include "dampingForce.h"

// include time stepper
#include "timeStepper.h"

// include input file and option
#include "setInput.h"

class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	void setRodStepper();
	void updateTimeStep();
	int simulationRunning();
	int numPoints();
	double getScaledCoordinate(int i);
	double getCurrentTime();
	double getTotalTime();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);
		
private:

	// Physical parameters
	double RodLength;
	double helixradius;
	double helixpitch;
	double rodRadius;
	int numVertices;
	double youngM;
	double Poisson;
	double shearM;
	double deltaTime;
	double totalTime;
	double density;
	Vector3d gVector;
	double viscosity;

	// Viscous drag coefficients
    double eta_per, eta_par;
    
	double tol, stol;
	int maxIter; // maximum number of iterations
	double characteristicForce;
	double forceTol;
	
	// Geometry
	MatrixXd vertices;
	
	// Rod
	elasticRod *rod;
	
	// set up the time stepper
	timeStepper *stepper;
	double *totalForce;
	double currentTime;
	
	// declare the forces
	elasticStretchingForce *m_stretchForce;
	elasticBendingForce *m_bendingForce;
	elasticTwistingForce *m_twistingForce;
	inertialForce *m_inertialForce;	
	externalGravityForce *m_gravityForce;
	dampingForce *m_dampingForce;
	
	int Nstep;
	int timeStep;
	int iter;

	void rodGeometry();
	void rodBoundaryCondition();

	// Variables about angular velocity
	double deltaTheta;
	string inputName; // name of file containing omega
	std::vector<double> timeSeries, omegaSeries;
    int numOmegaPoints, currentOmegaIndex;
    void ReadOmegaData();
    
	bool render; // should the OpenGL rendering be included?
	bool saveData; // should data be written to a file?

	double axisLengthInput;
	double deltaLengthInput;

	int start_s;
	int stop_s;
};

#endif
