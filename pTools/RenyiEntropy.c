/*******************************************************************************
** RenyiEntropy.c
** Part of the mutual information toolbox
**
** Contains functions to calculate the Renyi alpha entropy of a single variable 
** H_\alpha(X), the Renyi joint entropy of two variables H_\alpha(X,Y), and the 
** conditional Renyi entropy H_\alpha(X|Y)
** 
** Author: Adam Pocock
** Created 26/3/2010
** Updated - 22/02/2014 - Added checking on calloc.
**
** Copyright 2010-2017 Adam Pocock, The University Of Manchester
** www.cs.manchester.ac.uk
**
** This file is part of MIToolbox, licensed under the 3-clause BSD license.
*******************************************************************************/

#include "../pTools/MIToolbox.h"
#include "../pTools/ArrayOperations.h"
#include "../pTools/CalculateProbability.h"
#include "../pTools/RenyiEntropy.h"

double renyiEntropy(ProbabilityState state, double alpha) {
  double entropy = 0.0;
  double tempValue = 0.0;
  int i;
  
  /*H_\alpha(X) = 1/(1-alpha) * \log(\sum_x p(x)^alpha)*/
  printf("Printing entropy values:\n");
  for (i = 0; i < state.numStates; i++) {
    tempValue = state.probabilityVector[i];
    if (tempValue > 0) {
      entropy += pow(tempValue,alpha);
      printf("Pij ^ q Value at i=%d, is %f\n",i, pow(tempValue,alpha));
    }
  }
  
  printf("Sum of Pij ^ q = %f",entropy);
  entropy = log(entropy);
  printf(" after log = %f",entropy);
  entropy /= log(LOG_BASE);
  printf(" after base conversion = %f",entropy);
  entropy /= (1.0-alpha);
  printf(" divided by q-1 = %f\n",entropy);
  return entropy;
}

double calcRenyiEntropy(double alpha, uint *dataVector, int vectorLength) {
  ProbabilityState state = calculateProbability(dataVector,vectorLength);
  double h = renyiEntropy(state,alpha);
  
  freeProbabilityState(state);
  
  return h;
}/*calcRenyiEntropy(double,uint*,int)*/

double calcRenyiEntropyFromCoV(double alpha, uint *dataVector, int vectorLength) {
    printf("Calculating Renyi entropy from Co-occurrence vector for alpha (q) =  %f\n",alpha);

    double *stateProbs;
    double stateLength = vectorLength;
    ProbabilityState state;
    printf("Calculating probability :\n");
    stateProbs = (double *) checkedCalloc(stateLength,sizeof(double));

    int sumOfCo = sumState(dataVector, vectorLength);
    int i;
    for (i = 0; i < vectorLength; i++) {
       stateProbs[i] = dataVector[i]/(double)sumOfCo;
    }

    printf("Printing stateProbs\n");
    printDoubleVector(stateProbs,stateLength);

    state.probabilityVector = stateProbs;
    state.numStates = stateLength;

    double h = renyiEntropy(state,alpha);
    printf("Entropy Sum = %f\n",h);
    freeProbabilityState(state);
    
    return h;
}


double discAndCalcRenyiEntropy(double alpha, double *dataVector, int vectorLength) {
  ProbabilityState state = discAndCalcProbability(dataVector,vectorLength);
  double h = renyiEntropy(state,alpha);
  
  freeProbabilityState(state);
  
  return h;
}/*discAndCalcRenyiEntropy(double,double*,int)*/

double jointRenyiEntropy(JointProbabilityState state, double alpha) {
  double jointEntropy = 0.0;  
  double tempValue = 0.0;
  int i;
  
  /*H_\alpha(XY) = 1/(1-alpha) * log(2)(sum p(xy)^alpha)*/
  for (i = 0; i < state.numJointStates; i++) {
    tempValue = state.jointProbabilityVector[i];
    if (tempValue > 0) {
      jointEntropy += pow(tempValue,alpha);
    }
  }
  
  jointEntropy = log(jointEntropy);
  jointEntropy /= log(LOG_BASE);
  jointEntropy /= (1.0-alpha);

  return jointEntropy;
}

double calcJointRenyiEntropy(double alpha, uint *firstVector, uint *secondVector, int vectorLength) {
  JointProbabilityState state = calculateJointProbability(firstVector,secondVector,vectorLength);
  double h = jointRenyiEntropy(state,alpha);
  
  freeJointProbabilityState(state);
  
  return h;
}/*calcJointRenyiEntropy(double,uint*,uint*,int)*/

double discAndCalcJointRenyiEntropy(double alpha, double *firstVector, double *secondVector, int vectorLength) {
  JointProbabilityState state = discAndCalcJointProbability(firstVector,secondVector,vectorLength);
  double h = jointRenyiEntropy(state,alpha);
  
  freeJointProbabilityState(state);
  
  return h;
}/*discAndCalcJointRenyiEntropy(double,double*,double*,int)*/
