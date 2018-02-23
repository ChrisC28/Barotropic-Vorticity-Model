/*
 *  modelGrid.h
 *  BarotropicVorticityModel
 *
 *  Created by Christopher Chapman on 14/02/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MODELGRID_H
#define MODELGRID_H

class modelGrid	{

	
public:
	//////////////////////////////////////////////////////
	//Constructors
	//////////////////////////////////////////////////////
	
	modelGrid(int nX, int nY, double dX, double dY);  //Constructor: 
	modelGrid(int nX, int nY, double delta);  //Constructor:
	modelGrid(const modelGrid&);
	
	//////////////////////////////////////////////////////
	//Destructors
	//////////////////////////////////////////////////////
	
	~modelGrid();
	
	//////////////////////////////////////////////////////
	//Getters
	//////////////////////////////////////////////////////
	
	double getX(int iX);
	double getY(int iY);
	double getDeltaX();
	double getDeltaY();
	double nXgrid();
	double nYgrid();
	
private:
	
	int numX;
	int numY;
	double deltaX;
	double deltaY;
	
	double *xArray;
	double *yArray;
	
	

	
};

#endif