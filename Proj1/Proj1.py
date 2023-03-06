# libraries
import pandas as pd
import math
import random
import os

# function that generates posterior probabilities for a strictly biallelic model, assuming P(S = A) = P(S = B)
def posteriorProbs(df, A, B, pAA, pBB, pAB):
    # generate log liklihoods for each genotype
    llAA = 0
    llBB = 0
    llAB = 0
    for i in range(len(df)):
        O = df.at[i, 'observations']
        E = df.at[i, 'probability_of_error']
        if O == A:
            E = 1 - E
        llAA += math.log(E)
        llBB += math.log(1 - E) # all error probabilities are just the complement of the AA genotype's
        llAB += math.log(0.5) # read's likelihood simplifies to 0.5 under assumption that P(S = A) = P(S = B)
    # calculate denominator of Bayes' theorem equation
    d = math.exp(llAA + math.log(pAA)) + math.exp(llBB + math.log(pBB)) + math.exp(llAB + math.log(pAB))
    # calculate posterior probabilities for each genotype
    ppAA = math.exp(llAA + math.log(pAA) - math.log(d))
    ppBB = math.exp(llBB + math.log(pBB) - math.log(d))
    ppAB = math.exp(llAB + math.log(pAB) - math.log(d))
    return [ppAA, ppBB, ppAB]
  
# function that generates posterior probabilities for a model with/without transition probabilities, assuming P(S = A) = P(S = B)
def posteriorProbs2(df, A, B, pAA, pBB, pAB):
    # generate log liklihoods for each genotype
    llAA = 0
    llBB = 0
    llAB = 0
    bases = {"A" : 0, "C" : 1, "G" : 2, "T" : 3} # dictionary to help with indexing into transitions
    for i in range(len(df)):
        O = df.at[i, 'observations']
        E = df.at[i, 'probability_of_error']
        T_A = transitions.iat[bases[A], bases[O]] # probability of a transition from the major allele to the observation
        T_B = transitions.iat[bases[B], bases[O]] # probability of a transition from the minor allele to the observation
        # if the base is the major allele
        if O == A:
            llAA += math.log(1 - E)
            llBB += math.log(E * T_B)
            llAB += math.log((1 - E) * 0.5 + (E * T_B) * 0.5)
        # if the base is the minor allele
        elif O == B:
            llAA += math.log(E * T_A)
            llBB += math.log(1 - E)
            llAB += math.log((E * T_A) * 0.5 + (1 - E) * 0.5)
        # if the base is neither
        else:
            llAA += math.log(E * T_A)
            llBB += math.log(E * T_B)
            llAB += math.log((E * T_A) * 0.5 + (E * T_B) * 0.5)
    # calculate denominator of Bayes' theorem equation
    d = math.exp(llAA + math.log(pAA)) + math.exp(llBB + math.log(pBB)) + math.exp(llAB + math.log(pAB))
    # calculate posterior probabilities for each genotype
    ppAA = math.exp(llAA + math.log(pAA) - math.log(d))
    ppBB = math.exp(llBB + math.log(pBB) - math.log(d))
    ppAB = math.exp(llAB + math.log(pAB) - math.log(d))
    return [ppAA, ppBB, ppAB]
