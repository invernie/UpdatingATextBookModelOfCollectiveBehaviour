# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 10:06:44 2020

@author: FastRun
"""

from __future__ import division
import math, copy, importlib  # importlib is for module reloading, if necessary
import numpy as np
import matplotlib.pyplot as plt


# save global functions into variables for speed
unif = np.random.uniform
randint = np.random.randint
xp = math.exp
    

def templateOnly(seed,q,R2,stones = 3000, T = 5000, tint = 100, printing = "off", recordStoneMovement = 'off', showIm = "off"):

    np.random.seed(seed)
    
    #creates space and sets simulation variables
    N = 30
    a = 81
    
    space_matrix = np.zeros([a,a], dtype = 'int8')
    r_mat = np.zeros([a,a], dtype = 'float64')
    cos_mat = np.zeros([a,a], dtype = 'float64')
    sin_mat = np.zeros([a,a], dtype = 'float64')
    for x in range(a):
        for y in range(a):
            ##calculate cartesian coordinates from matrix coordinates. Space matrix refers the position of nodes and not cells.
            xcart = x - 40
            ycart = 40 - y
            r = math.sqrt(xcart*xcart+ycart*ycart)
            r_mat[x,y] = r
            if r > 0:
                cos_mat[x,y] = xcart/r
                sin_mat[x,y] = ycart/r #at 0;0 , sine and cosine have default values 0 and 0 (trigonometrically impossible)
    
    # Franks and Deneubourg (1997) model parameters, without effect of stone density
    R1 = 18
    Pmax = 0.35
    Dmax = 0.5
    tau = 0.025
    Qmax = 3

    # fitness function parameters
    beta = 1
    gamma = 1
    
    #Tot Ft calculated as average of Ft at t itervals
    Ft_list = np.zeros(int(T/tint), dtype = np.float)
    rbar_list = np.zeros(int(T/tint), dtype = np.float)
    CoV_list = np.zeros(int(T/tint), dtype = np.float)
    Rbar_list = np.zeros(int(T/tint), dtype = np.float)
    if recordStoneMovement == 'on':
        space_mat_list = list()

    #creates colony population
    ant_mat = np.zeros([N], dtype = 'float64')
    qi = int(round(N*q))
    ant_mat [:qi] = R2
    ant_mat [qi:] = R1
    np.random.shuffle(ant_mat)

    #pellet generation: 1000 stone items are generated at initialisation. The number of items does not change throughout the simulation.
    for _ in range(stones):
        assign = 0
        while (assign ==0):
            ax = randint(0,a)
            ay = randint(0,a)
            if space_matrix[ax,ay]<Qmax:
                space_matrix[ax,ay] += 1
                assign = 1



    #SIMULATION
    t = 0
    while t<T:
    
        for n in range(N):
            #at each cycle, each ant has one chance to pick up a pellet; if it picks up, it then repeats its action until the pellet is deposited
            Ropt = ant_mat[n]
            pickedup = 0
            deposited = 0
            while not pickedup:
                px = randint(0,a)
                py = randint(0,a)
                if  space_matrix[px,py] > 0:
                    rp = r_mat[px,py]
                    Pr = Pmax*(1-(1/(1+tau*(rp-Ropt)**2)))
                    if unif(0,1)<=Pr:
                        space_matrix[px,py]+= -1
                        pickedup = 1
                    else:
                        pickedup = 2
                    
            if pickedup == 1:
                while not deposited:
                    dx = randint(0,a)
                    dy = randint(0,a)
                    Qd = space_matrix[dx,dy]
                    if Qd < Qmax:
                        rd = r_mat[dx,dy]
                        Dr = (Dmax/(1+tau*(rd-Ropt)**2))
                        if unif(0,1)<=Dr:
                            space_matrix[dx,dy] += 1
                            deposited = 1
                        
        #calculate structure fitness at time step t
        ##calculate rbar
        stone_v = r_mat[space_matrix[:,:] > 0] #slice given as a list of extracted values
        rbar = sum(stone_v)/len(stone_v)
        ##calculate CoV
        CoV = np.std(stone_v)/rbar
        ##calculate mean resultant length Rbar
        l = len(stone_v)
        cos_v = cos_mat[space_matrix[:,:] > 0]
        sin_v = sin_mat[space_matrix[:,:] > 0]
        Rbar = math.sqrt(sum(cos_v)**2+sum(sin_v)**2)/l  # if for every stone there is a corresponding stone at the opposite coordinates on the circunference, the sum of both cos and sin is zero. i.e., a perfect circumference has an Rbar of zero.

        Ft = xp(-beta*CoV)*xp(-gamma*Rbar)
    
        t += 1
    
        if t%tint == 0:
            fi = int(t/tint)
            Ft_list [fi-1] = Ft 
            rbar_list [fi-1] = rbar
            CoV_list [fi-1] = CoV
            Rbar_list [fi-1] = Rbar
            if recordStoneMovement == 'on':
                space_mat_list.append(copy.deepcopy(space_matrix))

        if printing == "on":
            print(t)

    #plot matrix
    if showIm == "on":
        plt.imshow(space_matrix, cmap = 'gray');
        plt.colorbar()
    #Final Ft as average Ft across the last 1000 rounds of simulation
    fsteps = 10
    Ft = sum(Ft_list[-fsteps:])/len(Ft_list[-fsteps:])
    rbar = sum(rbar_list[-fsteps:])/len(rbar_list[-fsteps:])
    
    if recordStoneMovement == 'on':
        outObj = [Ft, rbar, Ft_list,rbar_list,CoV_list,Rbar_list,space_mat_list]
    else:
        outObj = [Ft, rbar, Ft_list,rbar_list,CoV_list,Rbar_list,space_matrix]
    
    return(outObj)



def FD(seed,q,R2, Pmax = 0.35, Fmin = 0.01, Fmax = 0.55, Dmax = 0.5, Gmin = 0.01, Gmax = 0.55, tau = 0.025, stones = 3000, T = 5000, tint = 100, printing = "off", delayedPopInc = "off", lateRopt = 24, recordStoneMovement = 'off', showIm = "off"):

    np.random.seed(seed)
    
    # save local functions into variables for speed
    def calcDr(rd,Ropt,Sd,Sc):
        if Sd < Sc:
            G = Gmin
        else:
            G = Gmax
        Dr = (1/(1+tau*(rd-Ropt)**2))*Dmax*G
        return(Dr)
    def calcPr(rp,Ropt,Sp,Sc):
        if Sp < Sc:
            F = Fmax
        else:
            F = Fmin
        Pr = Pmax*(1-(1/(1+tau*(rp-Ropt)**2)))*F
        return(Pr)
    
    
    #creates space and sets simulation variables
    N = 30    # n of workers
    a = 81    # 80 mm by 80 mm matrix; distance node-node = 1 mm
    
    R1 = 18   # allele carried by (1-q) proportion of workers
    # additional Franks and Deneubourg (1997) model parameters
    Qmax = 3  # max n of stones per node
    Sc = 6    # critical stone density
    
    # fitness function parameters
    beta = 1
    gamma = 1
    
    # record data at time intervals
    Ft_list = np.zeros(int(T/tint), dtype = np.float)
    rbar_list = np.zeros(int(T/tint), dtype = np.float)
    CoV_list = np.zeros(int(T/tint), dtype = np.float)
    Rbar_list = np.zeros(int(T/tint), dtype = np.float)
    if recordStoneMovement == 'on':
        space_mat_list = list()
    
    # space matrix and conversion to cartesian coordinates
    space_matrix = np.zeros([a,a], dtype = 'int8')
    r_mat = np.zeros([a,a], dtype = 'float64')
    cos_mat = np.zeros([a,a], dtype = 'float64')
    sin_mat = np.zeros([a,a], dtype = 'float64')
    for x in range(a):
        for y in range(a):
            ##calculate cartesian coordinates from matrix coordinates. Space matrix refers the position of nodes and not cells.
            xcart = x - 40
            ycart = 40 - y
            r = math.sqrt(xcart*xcart+ycart*ycart)
            r_mat[x,y] = r
            if r > 0:
                cos_mat[x,y] = xcart/r
                sin_mat[x,y] = ycart/r #at 0;0 , sine and cosine have default values 0 and 0 (trigonometrically impossible)
    
    #creates colony population and assigns alleles
    ant_mat = np.zeros([N], dtype = 'float64')
    qi = int(round(N*q))
    ant_mat [:qi] = R2
    ant_mat [qi:] = R1
    np.random.shuffle(ant_mat)

    #pellet generation: 1000 stone items are generated at initialisation. The number of items does not change throughout the simulation.
    for _ in range(stones):
        assign = 0
        while (assign ==0):
            ax = randint(0,a)
            ay = randint(0,a)
            if space_matrix[ax,ay]<Qmax:
                space_matrix[ax,ay] += 1
                assign = 1



    #SIMULATION
    t=0
    while t<T:
    
        for n in range(N):
            #at each cycle, each ant has one chance to pick up a pellet; if it picks up, it then repeats its action until the pellet is deposited
            Ropt = ant_mat[n]
            if delayedPopInc == "on" and t>=5000:  # for sims where we check for the effect of a delayed increase in worker pop size
                Ropt = lateRopt
            pickedup = 0
            deposited = 0
            while not pickedup:
                px = randint(0,a)
                py = randint(0,a)
                if  space_matrix[px,py] > 0:
                    rp = r_mat[px,py]
                    if px == 80:
                        xpr = 0
                        xpl = space_matrix[px-1,py]
                    elif px == 0:
                        xpl = 0
                        xpr = space_matrix[px+1,py]
                    else:
                        xpr = space_matrix[px+1,py]
                        xpl = space_matrix[px-1,py]
                    if py == 80:
                        ypd = 0
                        ypu = space_matrix[px,py-1]
                    elif py == 0:
                        ypu = 0
                        ypd = space_matrix[px,py+1]
                    else:
                        ypu = space_matrix[px,py-1]
                        ypd = space_matrix[px,py+1] 
                    Sp = xpl+xpr+ypd+ypu
                    Pr = calcPr(rp, Ropt, Sp, Sc)
                    if unif(0,1)<=Pr:
                        space_matrix[px,py]+= -1
                        pickedup = 1
                    else:
                        pickedup = 2
                    
            if pickedup == 1:
                while not deposited:
                    dx = randint(0,a)
                    dy = randint(0,a)
                    Qd = space_matrix[dx,dy]
                    if Qd < Qmax:
                        rd = r_mat[dx,dy]
                        if dx == 80:
                            xdr = 0
                            xdl = space_matrix[dx-1,dy]
                        elif dx == 0:
                            xdl = 0
                            xdr = space_matrix[dx+1,dy]
                        else:
                            xdr = space_matrix[dx+1,dy]
                            xdl = space_matrix[dx-1,dy]
                        if dy == 80:
                            ydd = 0
                            ydu = space_matrix[dx,dy-1] 
                        elif dy == 0:
                            ydu = 0
                            ydd = space_matrix[dx,dy+1]
                        else:
                            ydd = space_matrix[dx,dy+1]
                            ydu = space_matrix[dx,dy-1] 
                        Sd = xdl+xdr+ydd+ydu+Qd
                        Dr = calcDr(rd, Ropt, Sd, Sc)
                        if unif(0,1)<=Dr:
                            space_matrix[dx,dy] += 1
                            deposited = 1
                        
        #calculate structure fitness at time step t
        ##calculate rbar
        stone_v = r_mat[space_matrix[:,:] > 0] # locations where stones are present; note: slice given as a list of extracted values
        rbar = sum(stone_v)/len(stone_v)       # ave radius of those locations
        ##calculate CoV
        CoV = np.std(stone_v)/rbar
        ##calculate mean resultant length/circular spread Rbar
        l = len(stone_v)
        cos_v = cos_mat[space_matrix[:,:] > 0]
        sin_v = sin_mat[space_matrix[:,:] > 0]
        Rbar = math.sqrt(sum(cos_v)**2+sum(sin_v)**2)/l

        Ft = math.exp(-beta*CoV)*math.exp(-gamma*Rbar)
    
        t += 1
    
        if t%tint == 0:
            fi = int(t/tint)
            Ft_list [fi-1] = Ft 
            rbar_list [fi-1] = rbar
            CoV_list [fi-1] = CoV
            Rbar_list [fi-1] = Rbar
            if recordStoneMovement == 'on':
                space_mat_list.append(copy.deepcopy(space_matrix))

        if printing == "on":
            print(t)

    #plot matrix
    if showIm == "on":
        plt.imshow(space_matrix, cmap = 'gray');
        plt.colorbar()
    #Final Ft as average Ft across the last 1000 rounds of simulation
    fsteps = 10
    Ft = sum(Ft_list[-fsteps:])/len(Ft_list[-fsteps:])
    rbar = sum(rbar_list[-fsteps:])/len(rbar_list[-fsteps:])
    
    if recordStoneMovement == 'on':
        outObj = [Ft, rbar, Ft_list,rbar_list,CoV_list,Rbar_list,space_mat_list]
    else:
        outObj = [Ft, rbar, Ft_list,rbar_list,CoV_list,Rbar_list,space_matrix]
        
    return(outObj)  # output will be added to a list and then saved as a binary file


def FD_delayed(seed,q,R2, Pmax = 0.35, Fmin = 0.01, Fmax = 0.55, Dmax = 0.5, Gmin = 0.01, Gmax = 0.55, tau = 0.025, stones = 3000, printing = "off", delayedPopInc = "on", lateRopt = 24, showIm = "off"):

    np.random.seed(seed)
    
    # save local functions into variables for speed
    def calcDr(rd,Ropt,Sd,Sc):
        if Sd < Sc:
            G = Gmin
        else:
            G = Gmax
        Dr = (1/(1+tau*(rd-Ropt)**2))*Dmax*G
        return(Dr)
    def calcPr(rp,Ropt,Sp,Sc):
        if Sp < Sc:
            F = Fmax
        else:
            F = Fmin
        Pr = Pmax*(1-(1/(1+tau*(rp-Ropt)**2)))*F
        return(Pr)
    
    
    #creates space and sets simulation variables
    N = 30    # n of workers
    a = 81    # 80 mm by 80 mm matrix; distance node-node = 1 mm
    T = 75000
    
    R1 = 18   # allele carried by (1-q) proportion of workers
    # additional Franks and Deneubourg (1997) model parameters
    Qmax = 3  # max n of stones per node
    Sc = 6    # critical stone density
    
    # fitness function parameters
    beta = 1
    gamma = 1
    
    # record data at time intervals
    tint = 100
    Ft_list = np.zeros(int(T/tint), dtype = np.float)
    rbar_list = np.zeros(int(T/tint), dtype = np.float)
    CoV_list = np.zeros(int(T/tint), dtype = np.float)
    Rbar_list = np.zeros(int(T/tint), dtype = np.float)

    
    # space matrix and conversion to cartesian coordinates
    space_matrix = np.zeros([a,a], dtype = 'int8')
    r_mat = np.zeros([a,a], dtype = 'float64')
    cos_mat = np.zeros([a,a], dtype = 'float64')
    sin_mat = np.zeros([a,a], dtype = 'float64')
    for x in range(a):
        for y in range(a):
            ##calculate cartesian coordinates from matrix coordinates. Space matrix refers the position of nodes and not cells.
            xcart = x - 40
            ycart = 40 - y
            r = math.sqrt(xcart*xcart+ycart*ycart)
            r_mat[x,y] = r
            if r > 0:
                cos_mat[x,y] = xcart/r
                sin_mat[x,y] = ycart/r #at 0;0 , sine and cosine have default values 0 and 0 (trigonometrically impossible)
    
    #creates colony population and assigns alleles
    ant_mat = np.zeros([N], dtype = 'float64')
    qi = int(round(N*q))
    ant_mat [:qi] = R2
    ant_mat [qi:] = R1
    np.random.shuffle(ant_mat)

    #pellet generation: 1000 stone items are generated at initialisation. The number of items does not change throughout the simulation.
    for _ in range(stones):
        assign = 0
        while (assign ==0):
            ax = randint(0,a)
            ay = randint(0,a)
            if space_matrix[ax,ay]<Qmax:
                space_matrix[ax,ay] += 1
                assign = 1



    #SIMULATION
    t=0
    while t<T:
    
        for n in range(N):
            #at each cycle, each ant has one chance to pick up a pellet; if it picks up, it then repeats its action until the pellet is deposited
            Ropt = ant_mat[n]
            
            if t >= 5000:  # for sims where we check for the effect of a delayed increase in worker pop size
               Ropt = lateRopt
            
            pickedup = 0
            deposited = 0
            while not pickedup:
                px = randint(0,a)
                py = randint(0,a)
                if  space_matrix[px,py] > 0:
                    rp = r_mat[px,py]
                    if px == 80:
                        xpr = 0
                        xpl = space_matrix[px-1,py]
                    elif px == 0:
                        xpl = 0
                        xpr = space_matrix[px+1,py]
                    else:
                        xpr = space_matrix[px+1,py]
                        xpl = space_matrix[px-1,py]
                    if py == 80:
                        ypd = 0
                        ypu = space_matrix[px,py-1]
                    elif py == 0:
                        ypu = 0
                        ypd = space_matrix[px,py+1]
                    else:
                        ypu = space_matrix[px,py-1]
                        ypd = space_matrix[px,py+1] 
                    Sp = xpl+xpr+ypd+ypu
                    Pr = calcPr(rp, Ropt, Sp, Sc)
                    if unif(0,1)<=Pr:
                        space_matrix[px,py]+= -1
                        pickedup = 1
                    else:
                        pickedup = 2
                    
            if pickedup == 1:
                while not deposited:
                    dx = randint(0,a)
                    dy = randint(0,a)
                    Qd = space_matrix[dx,dy]
                    if Qd < Qmax:
                        rd = r_mat[dx,dy]
                        if dx == 80:
                            xdr = 0
                            xdl = space_matrix[dx-1,dy]
                        elif dx == 0:
                            xdl = 0
                            xdr = space_matrix[dx+1,dy]
                        else:
                            xdr = space_matrix[dx+1,dy]
                            xdl = space_matrix[dx-1,dy]
                        if dy == 80:
                            ydd = 0
                            ydu = space_matrix[dx,dy-1] 
                        elif dy == 0:
                            ydu = 0
                            ydd = space_matrix[dx,dy+1]
                        else:
                            ydd = space_matrix[dx,dy+1]
                            ydu = space_matrix[dx,dy-1] 
                        Sd = xdl+xdr+ydd+ydu+Qd
                        Dr = calcDr(rd, Ropt, Sd, Sc)
                        if unif(0,1)<=Dr:
                            space_matrix[dx,dy] += 1
                            deposited = 1
                        
        #calculate structure fitness at time step t
        ##calculate rbar
        stone_v = r_mat[space_matrix[:,:] > 0] # locations where stones are present; note: slice given as a list of extracted values
        rbar = sum(stone_v)/len(stone_v)       # ave radius of those locations
        ##calculate CoV
        CoV = np.std(stone_v)/rbar
        ##calculate mean resultant length/circular spread Rbar
        l = len(stone_v)
        cos_v = cos_mat[space_matrix[:,:] > 0]
        sin_v = sin_mat[space_matrix[:,:] > 0]
        Rbar = math.sqrt(sum(cos_v)**2+sum(sin_v)**2)/l

        Ft = math.exp(-beta*CoV)*math.exp(-gamma*Rbar)
    
        t += 1
    
        if t%tint == 0:
            fi = int(t/tint)
            Ft_list [fi-1] = Ft 
            rbar_list [fi-1] = rbar
            CoV_list [fi-1] = CoV
            Rbar_list [fi-1] = Rbar
            
        if t == 5000: 
            out5000 = [Ft, rbar, copy.deepcopy(Ft_list[0:50]),copy.deepcopy(rbar_list[0:50]),copy.deepcopy(CoV_list[0:50]),copy.deepcopy(Rbar_list[0:50]),copy.deepcopy(space_matrix)]
        if t == 6000:
            out6000 = [Ft, rbar, copy.deepcopy(Ft_list[0:60]),copy.deepcopy(rbar_list[0:60]),copy.deepcopy(CoV_list[0:60]),copy.deepcopy(Rbar_list[0:60]),copy.deepcopy(space_matrix)]
        if t == 7500:
            out7500 = [Ft, rbar, copy.deepcopy(Ft_list[0:75]),copy.deepcopy(rbar_list[0:75]),copy.deepcopy(CoV_list[0:75]),copy.deepcopy(Rbar_list[0:75]),copy.deepcopy(space_matrix)]
        if t == 10000:
            out10000 = [Ft, rbar, copy.deepcopy(Ft_list[0:100]),copy.deepcopy(rbar_list[0:100]),copy.deepcopy(CoV_list[0:100]),copy.deepcopy(Rbar_list[0:100]),copy.deepcopy(space_matrix)]
        if t == 20000:
            out20000 = [Ft, rbar, copy.deepcopy(Ft_list[0:200]),copy.deepcopy(rbar_list[0:200]),copy.deepcopy(CoV_list[0:200]),copy.deepcopy(Rbar_list[0:200]),copy.deepcopy(space_matrix)]
        if t == 50000:
            out50000 = [Ft, rbar, copy.deepcopy(Ft_list[0:500]),copy.deepcopy(rbar_list[0:500]),copy.deepcopy(CoV_list[0:500]),copy.deepcopy(Rbar_list[0:500]),copy.deepcopy(space_matrix)]
                
        if printing == "on":
            print(t)

    #plot matrix
    if showIm == "on":
        plt.imshow(space_matrix, cmap = 'gray');
        plt.colorbar()
    #Final Ft as average Ft across the last 1000 rounds of simulation
    fsteps = 10
    Ft = sum(Ft_list[-fsteps:])/len(Ft_list[-fsteps:])
    rbar = sum(rbar_list[-fsteps:])/len(rbar_list[-fsteps:])
    
    out75000 = [Ft, rbar, Ft_list,rbar_list,CoV_list,Rbar_list,space_matrix]
    
    outlist = [out5000, out6000, out7500, out10000, out20000, out50000, out75000]
    return(outlist)  # output will be added to a list and then saved as a binary file

