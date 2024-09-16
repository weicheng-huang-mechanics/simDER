#include "dampingForce.h"
#include <iostream>

dampingForce::dampingForce(elasticRod &m_rod, timeStepper &m_stepper, 
	double m_viscosity, double m_eta_per, double m_eta_par)
{
	rod = &m_rod;
	stepper = &m_stepper;
	viscosity = m_viscosity;
	dt = rod->dt;
	eta_per = m_eta_per;
	eta_par = m_eta_par;
		
	Id3<<1,0,0,
		0,1,0,
        0,0,1;
         
}

dampingForce::~dampingForce()
{
	;
}

void dampingForce::computeFd()
{
	for (int i=0;i<rod->ne;i++)
	{
		u = rod->getVelocity(i);
		// f = - viscosity * u *  rod->voronoiLen(i);
        t = rod->getTangent(i);
        
        f = -((eta_par - eta_per) * (t.dot(u)) * t + eta_per * u) * rod->voronoiLen(i);
		
		for (int k = 0; k < 3; k++)
		{
			ind = 4*i + k;
			stepper->addForce(ind, - f[k]); // subtracting external force
		}
	}
}

void dampingForce::computeJd()
{
	// Remember that dF/dx = 1/dt * dF/dv 
	for (int i=0;i<rod->ne;i++)
	{
		u = rod->getVelocity(i);
		t = rod->getTangent(i);
        // jac = - viscosity/dt * rod->voronoiLen(i) * u;
		
		jac = -((eta_par - eta_per) * t * t.transpose() + eta_per * Id3) * rod->voronoiLen(i) / dt;
		
		for (int kx = 0; kx < 3; kx++)
		{
			indx = 4 * i + kx;
			for (int ky = 0; ky < 3; ky++)
			{
				indy = 4 * i + ky;
				stepper->addJacobian(indx, indy, - jac(kx,ky)); // subtracting external force
			}
		}
	}
}
