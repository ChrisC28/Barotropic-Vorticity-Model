/*
 *  doubleArray.h
 *  BarotropicVorticityModel
 *
 *  Created by Christopher Chapman on 14/02/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
class doubleArray  {

public:
    
    //////////////////////////////////////////////////////
    //Constructors
    //////////////////////////////////////////////////////
	
    doubleArray(int nX, int nY);
    doubleArray(char * inputFileName, int nX, int nY);
    doubleArray(const doubleArray &a);  
    doubleArray();
	
    //////////////////////////////////////////////////////
    //Destructor
    //////////////////////////////////////////////////////
	
     ~doubleArray();
	
     //////////////////////////////////////////////////////
     //Setters
     //////////////////////////////////////////////////////
     double& operator()(int iX, int iY);
     doubleArray& operator=(const doubleArray& rhs);
	
     void setSize(int nX, int nY);
     //////////////////////////////////////////////////////
     //Getters
     //////////////////////////////////////////////////////

     int getNumX();
     int getNumY();
	
     //////////////////////////////////////////////////////
     //Utilities
     //////////////////////////////////////////////////////
	
      bool WriteData(char * fileName);
	
private:
	
    int numCols;
    int numRows;

    double * dataArray;
    void copyArray(const doubleArray& a);
	
};
