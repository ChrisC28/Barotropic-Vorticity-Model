/*
 *  modelGrid.cpp
 *  BarotropicVorticityModel
 *
 *  Created by Christopher Chapman on 14/02/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <stdlib.h>
#include "modelGrid.h"

modelGrid::modelGrid(int nX, int nY, double dX, double dY)	
{
	//////////////////////////////////////////////////////
	//Instantiate class variables
	//////////////////////////////////////////////////////
	
	deltaX = dX;
	deltaY = dY;
	
	numX = nX;
	numY = nY;
	
	//////////////////////////////////////////////////////
	//Set up model grid
	//////////////////////////////////////////////////////
	
	
	xArray = new double [numX];
	yArray = new double [numY];
	
	for(int iX=0; iX < numX; iX++)	{
		xArray[iX] = ((double) iX )* deltaX;
	}
	
	for(int iY=0; iY < numY; iY++)	{
		yArray[iY] = ((double) iY ) * deltaY;
	}
	
	
	
	
	
} //End: modelGrid(int nX, int nY, double dX, double dY);

modelGrid::modelGrid(int nX, int nY, double delta)	
{
	//////////////////////////////////////////////////////
	//Instantiate class variables
	//////////////////////////////////////////////////////
	
	deltaX = delta;
	deltaY = delta;
	
	numX = nX;
	numY = nY;
	
	//////////////////////////////////////////////////////
	//Set up model grid
	//////////////////////////////////////////////////////
	
	
	xArray = new double [numX];
	yArray = new double [numY];
	
	for(int iX=0; iX < numX; iX++)	{
		xArray[iX] = ((double) iX )* deltaX;
	}
	
	for(int iY=0; iY < numY; iY++)	{
		yArray[iY] = ((double) iY ) * deltaY;
	}
	
	
	
	
	
} //End: modelGrid(int nX, int nY, double delta);



modelGrid::~modelGrid()	
{

	delete [] yArray;
	delete [] xArray;
	
} //End: destructor


double modelGrid::getX(int iX)	{
	
	if((0 <= iX) && (numX > iX))
		return xArray[iX];
	else {
		return NULL;
	}


}  //End getX

double modelGrid::getY(int iY)	{
	
	if((0 <= iY) && (numY > iY))
		return xArray[iY];
	else {
		return NULL;
	}
	
	
}
double modelGrid::nXgrid()	{
	return numX;
}

double modelGrid::nYgrid()	{
	return numY;
}

double modelGrid::getDeltaX()	{
	return deltaX;
}

double modelGrid::getDeltaY()	{
	return deltaY;
}

