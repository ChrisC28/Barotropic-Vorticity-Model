#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "modelGrid.h"
#include "doubleArray.h"

//#include "modelGrid.cpp"

doubleArray Laplacian(doubleArray &inputField, modelGrid &grid);
doubleArray JacobianDeterminant(doubleArray &variableOne, doubleArray &variableTwo, modelGrid &grid);
doubleArray ForwardTimeStep(doubleArray &vorticityStepNmin1, doubleArray &vorticityStepN, doubleArray &streamFunction, modelGrid &grid, double deltaT, const double beta);
void SolvePoissonSOR(doubleArray &variable, doubleArray &rhs, modelGrid &grid, double omega, double relTol);

const double pi   = 3.1415926535;
int main (int argc, char * const argv[]) {

    using namespace std;
    // Set up input grid
	
    double deltaX = 200e3;
    double deltaY = 200e3;
    int nX	      = 101;
    int nY        = 51;
    const double beta = +1.0e-11; 
    //const double beta = 0.0;
    double relTol = 1.0e-7;
	
    const double deltaT = 15 * 60;
    char inputStreamFunctionFileName[256];
    strcpy(inputStreamFunctionFileName, "/home/cchlod/Barotropic_Vorticity_Model/BarotropicVorticityModel/barotropic_streamfunction_input.dat");
	
    char outputFileNameBase[256]; 
    strcpy(outputFileNameBase, "/home/cchlod/Barotropic_Vorticity_Model/BarotropicVorticityModel/output/streamFunctionOut_");
    
    double forecastTime = 144 * 3600;
    int numTimeSteps = forecastTime / deltaT;
	
    /////////////////////////////////////////////////////
    // Instantiation
    /////////////////////////////////////////////////////
    modelGrid grid(nX, nY, deltaX, deltaY);
    cout << "Read input" << endl; 	
    doubleArray streamFunction(inputStreamFunctionFileName, grid.nXgrid(), grid.nYgrid());
	
	
    cout << "input file name:........ " << inputStreamFunctionFileName << endl;
    cout << "num x grid points:...... " << grid.nXgrid() << endl;
    cout << "num y grid points:...... " << grid.nYgrid() << endl;
    cout << "time step: .............." << deltaT << " seconds" << endl;
    cout << "delta x: ............... " << grid.getDeltaX() << endl;
    cout << "delta y: ............... " << grid.getDeltaY() << endl;
    cout << "number of time steps.... " << numTimeSteps << endl;
	
	
    /////////////////////////////////////////////////////
    // Initialisation
    /////////////////////////////////////////////////////
	
    //Initialise vorticity field
    cout << "Initialising Fields" << endl;
    doubleArray vorticity(grid.nXgrid(), grid.nYgrid());
    doubleArray vorticityStep(vorticity);
    
    vorticity = Laplacian(streamFunction, grid);

    double omega;
	
    double r = 0.5 * (cos(pi/grid.nXgrid()) + cos(pi/grid.nXgrid()));
    omega = 2.0 / (1 + sqrt(1-(r*r)));
    omega = 1.5;
    cout << "Using omega value: " << omega << endl;
	
		
    cout << "========================" << endl;
    cout << "Completed Initialisation" << endl;
    cout << "========================" << endl;
    
    /////////////////////////////////////////////////////
    // Take first time step
    /////////////////////////////////////////////////////
		
    vorticityStep = ForwardTimeStep(vorticity, vorticity, streamFunction, grid, deltaT, beta);
    SolvePoissonSOR(streamFunction, vorticity, grid, omega, relTol);
	
    /////////////////////////////////////////////////////
    // Solve for new stream function
    /////////////////////////////////////////////////////
	

    //SolvePoissonSOR(streamFunction, vorticityStep, grid, omega,  relTol);
	
    // Begin Main Loop	
    doubleArray vorticityStepPlus(vorticityStep);
    //Main timing loop
    for(int iT=0; iT <= numTimeSteps; iT++)	{
        
        vorticityStepPlus = ForwardTimeStep(vorticity, vorticityStep, streamFunction, grid, 2*deltaT, beta);
	SolvePoissonSOR(streamFunction, vorticityStepPlus, grid, omega,  relTol);
		
	vorticity = vorticityStep;
	vorticityStep = vorticityStepPlus;
		
	cout << "Completed: " << 100 * (double) iT/ (double) numTimeSteps << "\% \r";
	
	if(iT % 12 == 0 )	{
            char outputFileName[128];
	    sprintf(outputFileName, "%s%d.dat",outputFileNameBase, iT/12);
	    cout << "Writing data to: " << outputFileName << endl;
	    streamFunction.WriteData(outputFileName);
			
	    }  //End write routine
	} //END main time steping loop
	
    return 0;
} //END: main model function

doubleArray Laplacian(doubleArray &inputField, modelGrid &grid)	{
    using namespace std;
    double laplacianUpdate;
    doubleArray lhs(grid.nXgrid(), grid.nYgrid());
	
    for(int iY=1; iY < grid.nYgrid()-1; iY++)	{
        for(int iX=1; iX < grid.nXgrid()-1; iX++)	{
	    laplacianUpdate  = ((inputField(iX+1,iY) + inputField(iX-1,iY) - (2.0 * inputField(iX,iY))) / (grid.getDeltaX() * grid.getDeltaX()));
	    laplacianUpdate += ((inputField(iX,iY+1) + inputField(iX,iY-1) - (2.0 * inputField(iX,iY))) / (grid.getDeltaY() * grid.getDeltaY()));
	    lhs(iX,iY) = laplacianUpdate;
	    }
    }
    
    return lhs;
	
} //END Laplacian


doubleArray ForwardTimeStep(doubleArray &vorticityStepNmin1, doubleArray &vorticityStepN, doubleArray &streamFunction, modelGrid &grid, double deltaT, const double beta) {
    
    doubleArray vorticityStepNplus1(vorticityStepN);
    doubleArray detJ(grid.nXgrid(), grid.nYgrid());
    detJ = JacobianDeterminant(streamFunction, vorticityStepN, grid);
    double vorticityStep;
    
    for(int iY=1; iY < grid.nYgrid()-1; iY++)	{
	    
        vorticityStep  = vorticityStepNmin1(0,iY);
	vorticityStep -= deltaT * (detJ(0,iY) + (beta*(streamFunction(1,iY)-streamFunction(grid.nXgrid()-1,iY))/(2*grid.getDeltaX())));
	vorticityStepNplus1(0,iY) = vorticityStep;
		
	for(int iX=1; iX < grid.nXgrid()-1; iX++)	{
	    vorticityStep  = vorticityStepNmin1(iX,iY);
	    vorticityStep -= deltaT * (detJ(iX,iY) + (beta*(streamFunction(iX+1,iY)-streamFunction(iX-1,iY))/(2*grid.getDeltaX())));
	    vorticityStepNplus1(iX,iY) = vorticityStep;
	}
	vorticityStep  = vorticityStepNmin1(grid.nXgrid()-1,iY);
	vorticityStep -= deltaT * (detJ(grid.nXgrid()-1,iY) + (beta*(streamFunction(0,iY)-streamFunction(grid.nXgrid()-2,iY))/(2*grid.getDeltaX())));
	vorticityStepNplus1(grid.nXgrid()-1,iY) = vorticityStep;
	
	
    }
	
    return vorticityStepNplus1;
} //END ForwardTimeStep


doubleArray JacobianDeterminant(doubleArray &variableOne, doubleArray &variableTwo, modelGrid &grid)	{
	
    doubleArray detJ(grid.nXgrid(), grid.nYgrid());
	
    double J1, J2, J3, Jtotal;
	
    for(int iY=1; iY < grid.nYgrid()-1; iY++)	{
        //Compute the Jacobian at the boundaries (using periodic BCs)
	J1  = (variableOne(1,iY) - variableOne(grid.nXgrid()-1,iY)) * (variableTwo(0,iY+1) - variableTwo(0,iY-1));
	J1 -= (variableTwo(1,iY) - variableTwo(grid.nXgrid()-1,iY)) * (variableOne(0,iY+1) - variableOne(0,iY-1));
		
	J2  = (variableOne(1,iY) * (variableTwo(1,iY+1) - variableTwo(1,iY-1))) - (variableOne(grid.nXgrid()-1,iY) * (variableTwo(grid.nXgrid()-1,iY+1) - variableTwo(grid.nXgrid()-1,iY-1)));
	J2 -= (variableOne(0,iY+1) * (variableTwo(1,iY+1) - variableTwo(grid.nXgrid()-1,iY+1))) - (variableOne(0,iY-1) * (variableTwo(1,iY-1) - variableTwo(grid.nXgrid()-1,iY-1)));
		
		
	J3  = (variableTwo(0,iY+1) * (variableOne(1,iY+1) - variableOne(grid.nXgrid()-1,iY+1))) - (variableTwo(0,iY-1) * (variableOne(1,iY-1) - variableOne(grid.nXgrid()-1,iY-1)));
	J3 -= (variableTwo(1,iY) * (variableOne(1,iY+1) - variableOne(1,iY-1))) - (variableTwo(grid.nXgrid()-1,iY) * (variableOne(grid.nXgrid()-1,iY+1) - variableOne(grid.nXgrid()-1,iY-1)));
		
		
	Jtotal = (1.0/3.0) * (J1 + J2 + J3);
	Jtotal /= 4.0 *(grid.getDeltaX() * grid.getDeltaY());
	detJ(0,iY) = Jtotal;
		
		
	for(int iX=1; iX < grid.nXgrid()-1; iX++)	{
			
	    J1  = (variableOne(iX+1,iY) - variableOne(iX-1,iY)) * (variableTwo(iX,iY+1) - variableTwo(iX,iY-1));
	    J1 -= (variableTwo(iX+1,iY) - variableTwo(iX-1,iY)) * (variableOne(iX,iY+1) - variableOne(iX,iY-1));
			
			
            J2  = (variableOne(iX+1,iY) * (variableTwo(iX+1,iY+1) - variableTwo(iX+1,iY-1))) - (variableOne(iX-1,iY) * (variableTwo(iX-1,iY+1) - variableTwo(iX-1,iY-1)));
	    J2 -= (variableOne(iX,iY+1) * (variableTwo(iX+1,iY+1) - variableTwo(iX-1,iY+1))) - (variableOne(iX,iY-1) * (variableTwo(iX+1,iY-1) - variableTwo(iX-1,iY-1)));
			
			
	    J3  = (variableTwo(iX,iY+1) * (variableOne(iX+1,iY+1) - variableOne(iX-1,iY+1))) - (variableTwo(iX,iY-1) * (variableOne(iX+1,iY-1) - variableOne(iX-1,iY-1)));
	    J3 -= (variableTwo(iX+1,iY) * (variableOne(iX+1,iY+1) - variableOne(iX+1,iY-1))) - (variableTwo(iX-1,iY) * (variableOne(iX-1,iY+1) - variableOne(iX-1,iY-1)));
			
			
	    Jtotal = (1.0/3.0) * (J1 + J2 + J3);
	    Jtotal /= 4.0 *(grid.getDeltaX() * grid.getDeltaY());
	    detJ(iX,iY) = Jtotal;
			
	} //END for (iX=1;iX<grid.nXgrid()-1;iX++)

        J1  = (variableOne(0,iY) - variableOne(grid.nXgrid()-2,iY)) * (variableTwo(grid.nXgrid()-1,iY+1) - variableTwo(grid.nXgrid()-1,iY-1));
	J1 -= (variableTwo(0,iY) - variableTwo(grid.nXgrid()-2,iY)) * (variableOne(grid.nXgrid()-1,iY+1) - variableOne(grid.nXgrid()-1,iY-1));
		
		
	J2  = (variableOne(0,iY) * (variableTwo(0,iY+1) - variableTwo(0,iY-1))) - (variableOne(grid.nXgrid()-2,iY) * (variableTwo(grid.nXgrid()-2,iY+1) - variableTwo(grid.nXgrid()-2,iY-1)));
	J2 -= (variableOne(grid.nXgrid()-1,iY+1) * (variableTwo(0,iY+1) - variableTwo(grid.nXgrid()-2,iY+1))) - (variableOne(grid.nXgrid()-1,iY-1) * (variableTwo(0,iY-1) - variableTwo(grid.nXgrid()-2,iY-1)));
		
		
	J3  = (variableTwo(grid.nXgrid()-1,iY+1) * (variableOne(0,iY+1) - variableOne(grid.nXgrid()-2,iY+1))) - (variableTwo(grid.nXgrid()-1,iY-1) * (variableOne(0,iY-1) - variableOne(grid.nXgrid()-2,iY-1)));
	J3 -= (variableTwo(0,iY) * (variableOne(0,iY+1) - variableOne(0,iY-1))) - (variableTwo(grid.nXgrid()-2,iY) * (variableOne(grid.nXgrid()-2,iY+1) - variableOne(grid.nXgrid()-2,iY-1)));
		
		
	Jtotal = (1.0/3.0) * (J1 + J2 + J3);
	Jtotal /= 4.0 *(grid.getDeltaX() * grid.getDeltaY());
	detJ(grid.nXgrid()-1,iY) = Jtotal;
	
	}
	
     return detJ;
}


void SolvePoissonSOR(doubleArray &variable, doubleArray &rhs, modelGrid &grid, double omega, double relTol)	{
    
    ///////////////////////////////////////////////////////////////////
    //Solve a poisson equation using Successive Over Relaxiation
    //The equation is of the form: div . grad(variable) = rhs
    //////////////////////////////////////////////////////////////////

    double dTau = 0.5 * grid.getDeltaX()*grid.getDeltaX()*grid.getDeltaY()*grid.getDeltaY();
    dTau /= (grid.getDeltaX()*grid.getDeltaX() + grid.getDeltaY()*grid.getDeltaY());
		  
    int maxIterations = grid.nXgrid() * grid.nYgrid();
	
    double variableStep;
    for(int iT=0; iT < maxIterations; iT++)	{
	double changeSum=0;
	//Loop over grid, using an implicit procedure 
	for(int iY=1; iY < grid.nYgrid()-1; iY++)	{
	    variableStep  = variable(0,iY);

	    variableStep +=omega * (dTau/(grid.getDeltaX()*grid.getDeltaX()))*(variable(1,iY) + variable(grid.nXgrid()-1,iY)-(2.0*variable(0,iY)));
	    variableStep +=omega * (dTau/(grid.getDeltaX()*grid.getDeltaX()))*(variable(0,iY+1) + variable(0,iY-1)-(2.0*variable(0,iY)));
	    variableStep -= dTau * omega * rhs(0,iY);
			
	    variable(0,iY) = variableStep;
			
	    for(int iX=1; iX < grid.nXgrid()-1; iX++)	{
	        variableStep  = variable(iX,iY);
		variableStep += omega * (dTau/(grid.getDeltaX()*grid.getDeltaX()))*(variable(iX+1,iY) + variable(iX-1,iY)-(2.0*variable(iX,iY)));
		variableStep += omega * (dTau/(grid.getDeltaY()*grid.getDeltaY()))*(variable(iX,iY+1) + variable(iX,iY-1)-(2.0*variable(iX,iY)));
		variableStep -= dTau * omega * rhs(iX,iY);
				
		//Compute total residual over domain
				
		changeSum += fabs(1-variable(iX,iY)/variableStep);
		variable(iX,iY) = variableStep;
				
		} //END for iX
	
        variableStep  = variable(grid.nXgrid()-1,iY);
	variableStep +=omega * (dTau/(grid.getDeltaX()*grid.getDeltaX()))*(variable(0,iY) + variable(grid.nXgrid()-2,iY)-(2.0*variable(grid.nXgrid()-1,iY)));
	variableStep +=omega * (dTau/(grid.getDeltaX()*grid.getDeltaX()))*(variable(grid.nXgrid()-1,iY+1) + variable(grid.nXgrid()-1,iY-1)-(2.0*variable(grid.nXgrid()-1,iY)));
	variableStep -= dTau * omega * rhs(grid.nXgrid()-1,iY);
			
	variable(grid.nXgrid()-1,iY) = variableStep;
		
		
	} //END for iY 
		
	changeSum /= ((grid.nXgrid()-2) * (grid.nYgrid()-2));
	//	if(((iT % 10) < 1 ) && (iT != 0))	{
	//		std::cout<< "After " << iT << "iterations, fractional change = " << changeSum << std::endl; 
	//	}
		
	if(changeSum < relTol)	{
	    std::cout << "Convergence after: " << iT << "iterations" << std::endl;
	    break;
	}
    
    } //End outer loops
	
} //END SolvePoissonSOR
 
