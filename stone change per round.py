# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 14:37:58 2021

@author: FastRun
"""
import csv, pickle, copy, importlib
import numpy as np
import matplotlib.pyplot as plt
import wall_funcs as w


T = 5000
stones = 3000

# Run a set of simulations where the changes in space matrix is tracked round by round
outList= list()
reps = 20
for seed in range(reps):
    # record data every round
    out = w.FD(seed,0,18, T = 5000, stones = stones, tint = 1, recordStoneMovement = 'on')
    outList.append(out) # list of list - each sublist is the sequence of space matrices for a single simulation

simDtList = list()
stoneMovedMat = np.zeros([reps, T-1])
for sim in range(reps): # for each simulation
    matList = outList[sim][6] # extract space_mat list
    
    stoneMoved = 0
    previous = matList[0] # we do not have the space_mat at initialisation, so we use the end of round 0 as t0
    for tx in range(1,T):
        current = matList[tx] # space_mat for each year, starting from round 1
        for x in range(81):
            for y in range(81): # for each spatial location in the space_mat
                if current[x][y] != previous[x][y]: # check if the number of stones is identical between the two rounds
                    stoneMoved += 1
        stoneMovedMat [sim][tx-1] = stoneMoved/2 # each stone moved changes the value of two space locations (stones are always deposited within the same simulation round as they are picked up)
        stoneMoved = 0
        previous = copy.deepcopy(current) # save a true copy of the current round to become the comparison for the next round
        
TMat = np.transpose(stoneMovedMat) # transpose the matrix so that each row contains the values of one round for all simulations
yearlyAve = [ np.mean(y) for y in TMat]
yearlySD = [ np.std(y) for y in TMat]

csvfile = "Path/To/File/stone_movement.csv"
np.savetxt(csvfile, TMat, delimiter = ',', header = "Seed0,Seed1,Seed2,Seed3,Seed4,Seed5,Seed6,Seed7,Seed8,Seed9,Seed10,Seed11,Seed12,Seed13,Seed14,Seed15,Seed16,Seed17,Seed18,Seed19")

filepath = "Path/To/File/deposition_events_per_round_FD_stones3000T5000"
objToSave = [TMat, outList, yearlyAve, yearlySD]
file = open(filepath, 'wb')
pickle.dump(objToSave, file)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
file.close()

plt.plot(range(1,T), yearlyAve)
plt.plot(range(1,T), yearlySD)

plt.plot(range(1,T), stoneMovedMat[0])



    