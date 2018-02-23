/*
 *  doubleArray.cpp
 *  BarotropicVorticityModel
 *
 *  Created by Christopher Chapman on 14/02/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "doubleArray.h"

//////////////////////////////////////////////////////
//Constructors
//////////////////////////////////////////////////////


doubleArray::doubleArray(int nX, int nY)	{
	
    numCols = nX;
    numRows = nY;
	
    dataArray = new double [numRows * numCols];
	
     for(int index=0; index < numRows*numCols; index++)	{
		dataArray[index] = 0; 
      }
	
	
	
} //End: doubleArray(int nX, int nY);

doubleArray::doubleArray()	{

    numRows = 1;
    numCols = 1;
	
    dataArray = new double [1];
    dataArray[0] = 0;

} //End doubleArray()


doubleArray::doubleArray(char * inputFileName, int nX, int nY)	{
    using namespace std;
	
    FILE *inputFILE = fopen(inputFileName,"r");
    float inputNumber;
    numCols = nX;
    numRows = nY;
    cout << "declare array" << endl;
    dataArray = new double [numRows * numCols];
    cout << inputFileName << endl;


	
    for(int iY=0; iY < numRows; iY++)	{
        for(int iX=0; iX < numCols; iX++)	{
            cout << "iX: " << iX << endl;		
	     
            fscanf(inputFILE, "%f", &inputNumber);	
	    cout << "number: " << inputNumber;
            dataArray[((iY*numCols) + iX)] = inputNumber;
	}
    }

    fclose(inputFILE);
    return;
} //End doubleArray(char * inputFileName)


doubleArray::doubleArray(const doubleArray &a)	{
	
    numRows = a.numRows;
    numCols = a.numCols;
    dataArray = new double [numRows * numCols];
	
    this->copyArray(a);
	
}  //End doubleArray(const doubleArray& a)

//////////////////////////////////////////////////////
//Destructor
//////////////////////////////////////////////////////


doubleArray:: ~doubleArray()	{
	
    delete [] dataArray;
	
} //End destructor

//////////////////////////////////////////////////////
//Setters
//////////////////////////////////////////////////////


double& doubleArray :: operator()(int iX, int iY)	{
	
	 
    return dataArray[((iY*numCols) + iX)];
	
}

doubleArray& doubleArray :: operator=(const doubleArray& rhs)	{
	
    if(dataArray != rhs.dataArray)	{
        setSize(rhs.numCols, rhs.numRows);
	copyArray(rhs);
    }
	
    return *this;
} //End 


void doubleArray :: setSize(int nX, int nY)	{
    
    delete [] dataArray;
    numRows = nY;
    numCols = nX;
    dataArray = new double [numRows * numCols];
}

//////////////////////////////////////////////////////
//Getters
//////////////////////////////////////////////////////


int doubleArray :: getNumX()	{

    return numCols;

}

int doubleArray :: getNumY()	{
		
    return numRows;

}

void doubleArray :: copyArray(const doubleArray& a)	{

		
    double *temp = dataArray + (numCols * numRows);
    double *input = a.dataArray + (numCols * numRows); 
	
    while(temp>dataArray) {
	*--temp = *--input;
    }
}

//////////////////////////////////////////////////////
//Utilities
//////////////////////////////////////////////////////


bool doubleArray :: WriteData(char * fileName)	{
	
    bool outputGood = true;
    using namespace std;
	
    ofstream dataFile;
    dataFile.open(fileName, ios::out);
    if(dataFile.good())	{
        for(int iY =0; iY < numRows; iY++)	{
	    for(int iX=0; iX < numCols; iX++)	{
		dataFile << dataArray[(iY*numCols)+iX];
		    if(iX != numCols-1)	{
			dataFile << ", "; 
		    }
		}
            dataFile << "\n";
	    }
    }
	
    else {
        outputGood = false;
    }

    dataFile.close();
    return outputGood;
	
}// End WriteData
