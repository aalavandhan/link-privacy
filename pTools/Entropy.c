/*******************************************************************************
 ** Entropy.c
 ** Part of the mutual information toolbox
 **
 ** Contains functions to calculate the entropy of a single variable H(X), 
 ** the joint entropy of two variables H(X,Y), and the conditional entropy
 ** H(X|Y)
 ** 
 ** Author: Adam Pocock
 ** Created 19/2/2010
 **
 ** Copyright 2010-2017 Adam Pocock, The University Of Manchester
 ** www.cs.manchester.ac.uk
 **
 ** This file is part of MIToolbox, licensed under the 3-clause BSD license.
 ******************************************************************************/

#include "../pTools/MIToolbox.h"
#include "../pTools/ArrayOperations.h"
#include "../pTools/CalculateProbability.h"
#include "../pTools/Entropy.h"

double entropy(ProbabilityState state) {
    double entropy = 0.0;
    double tempValue = 0.0;
    int i;

    /*H(X) = - \sum p(x) \log p(x)*/
    printf("Printing entropy values:\n");
    for (i = 0; i < state.numStates; i++) {
        tempValue = state.probabilityVector[i];

        if (tempValue > 0) {
            entropy -= tempValue * log(tempValue);
            printf("Entropy Value at i=%d, is %f\n",i,tempValue * log(tempValue));
        }
    }

    entropy /= log(LOG_BASE);

    return entropy;
}

double discAndCalcEntropy(double* dataVector, int vectorLength) {
    ProbabilityState state = discAndCalcProbability(dataVector, vectorLength);
    double h = entropy(state);

    freeProbabilityState(state);

    return h;
}/*discAndCalcEntropy(double* ,int)*/

double calcEntropy(uint* dataVector, int vectorLength) {
    printf("Calculating entropy: \n");
    ProbabilityState state = calculateProbability(dataVector, vectorLength);
    double h = entropy(state);
    printf("Entropy Sum = %f\n",h);
    freeProbabilityState(state);
    
    return h;
}/*calcEntropy(uint* ,int)*/

double calcEntropyFromCoV(uint* dataVector, int vectorLength) {
    printf("Calculating entropy from Co-occurrence vector: \n");

    double *stateProbs;
    double stateLength = vectorLength;
    ProbabilityState state;
    printf("Calculating probability : \n");
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

    double h = entropy(state);
    printf("Entropy Sum = %f\n",h);
    freeProbabilityState(state);
    
    return h;
}/*calcEntropy(uint* ,int)*/

double calcEntropyFromLocationVector(uint* dataVector, int vectorLength) {
    printf("Calculating entropy from Co-occurrence vector: \n");

    double *stateProbs;
    double stateLength = vectorLength;
    ProbabilityState state;
    printf("Calculating probability : \n");
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

    double h = entropy(state);
    printf("Entropy Sum = %f\n",h);
    freeProbabilityState(state);
    
    return h;
}/*calcEntropy(uint* ,int)*/




double jointEntropy(JointProbabilityState state) {
    double jointEntropy = 0.0;
    double tempValue = 0.0;
    int i;
    
    /*H(XY) = - \sum_x \sum_y p(xy) \log p(xy)*/
    for (i = 0; i < state.numJointStates; i++) {
        tempValue = state.jointProbabilityVector[i];
        if (tempValue > 0) {
            jointEntropy -= tempValue * log(tempValue);
        }
    }

    jointEntropy /= log(LOG_BASE);

    return jointEntropy;
}

double discAndCalcJointEntropy(double *firstVector, double *secondVector, int vectorLength) {
    JointProbabilityState state = discAndCalcJointProbability(firstVector, secondVector, vectorLength);
    double h = jointEntropy(state);

    freeJointProbabilityState(state);

    return h;
}/*discAndCalcJointEntropy(double *, double *, int)*/

double calcJointEntropy(uint *firstVector, uint *secondVector, int vectorLength) {
    JointProbabilityState state = calculateJointProbability(firstVector, secondVector, vectorLength);
    double h = jointEntropy(state);

    freeJointProbabilityState(state);

    return h;
}/*calcJointEntropy(uint *, uint *, int)*/

double condEntropy(JointProbabilityState state) {
    double condEntropy = 0.0;
    double jointValue = 0.0;
    double condValue = 0.0;
    int i;

    /*H(X|Y) = - \sum_x \sum_y p(x,y) \log p(x,y)/p(y)*/
    /** Indexing by numFirstStates use modulus of i
     ** Indexing by numSecondStates use integer division of i by numFirstStates*/
    for (i = 0; i < state.numJointStates; i++) {
        jointValue = state.jointProbabilityVector[i];
        condValue = state.secondProbabilityVector[i / state.numFirstStates];
        if ((jointValue > 0) && (condValue > 0)) {
            condEntropy -= jointValue * log(jointValue / condValue);
        }
    }

    condEntropy /= log(LOG_BASE);

    return condEntropy;
}

double discAndCalcConditionalEntropy(double *dataVector, double *conditionVector, int vectorLength) {
    JointProbabilityState state = discAndCalcJointProbability(dataVector, conditionVector, vectorLength);
    double h = condEntropy(state);

    freeJointProbabilityState(state);

    return h;
}/*discAndCalcConditionalEntropy(double *, double *, int)*/

double calcConditionalEntropy(uint *dataVector, uint *conditionVector, int vectorLength) {
    JointProbabilityState state = calculateJointProbability(dataVector, conditionVector, vectorLength);
    double h = condEntropy(state);

    freeJointProbabilityState(state);

    return h;
}/*calcConditionalEntropy(uint *, uint *, int)*/

