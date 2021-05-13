# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 23:50:25 2020

@author: FastRun
"""
import pickle  # use pickle because shelve: 1. splits db in three objects with different extensions in win10 and there are issues calling the file back and 2. is not portbale across OS
import numpy as np
import wall_funcs as w

stones = 3000

filepathFD = "Path/To/File/out_FD_" + str(stones) +"stones"
filepathT = "Path/To/File/out_T_" + str(stones) + "stones"

S = 100

outListFD = list()
outListT = list()        
        
for s in range(S):
    outListFD.append(w.FD(s,0,18, stones=stones))
    outListT.append(w.templateOnly(s,0,18, stones=stones))
    print("FD sim is at:" + str(s))
    
fileFD = open(filepathFD, 'wb')
pickle.dump(outListFD, fileFD)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
fileFD.close()
        
fileT = open(filepathT, 'wb')
pickle.dump(outListT, fileT)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
fileT.close()



# save stats for all 3 stone availability versions in a csv file
# plot data as boxplits using R
# file 1, calculate mean and std
import csv, pickle
import numpy as np
import matplotlib.pyplot as plt

fsteps = 10
## FD files
loadpath = "Path/To/File/out_FD_1000stones"
fileFD1000 = open(loadpath, 'rb')
allResFD1000 = pickle.load(fileFD1000)
fileFD1000.close()
CoVList_FD1000 = [sum(f[4][-fsteps:])/fsteps for f in allResFD1000]
aveCoV_FD1000 = np.mean(CoVList_FD1000)
stdCoV_FD1000 = np.std(CoVList_FD1000)
RbarList_FD1000 = [sum(f[5][-fsteps:])/fsteps for f in allResFD1000]
aveRbar_FD1000 = np.mean(RbarList_FD1000)
stdRbar_FD1000 = np.std(RbarList_FD1000)
aveFt_FD1000 = np.mean([f[0] for f in allResFD1000])
stdFt_FD1000 = np.std([f[0] for f in allResFD1000])
averbar_FD1000 = np.mean([f[1] for f in allResFD1000])
stdrbar_FD1000 = np.std([f[1] for f in allResFD1000])
imFD1000 = allResFD1000[99][6]
    
loadpath = "Path/To/File/out_FD_3000stones"
fileFD3000 = open(loadpath, 'rb')
allResFD3000 = pickle.load(fileFD3000)
fileFD3000.close()
CoVList_FD3000 = [sum(f[4][-fsteps:])/fsteps for f in allResFD3000]
aveCoV_FD3000 = np.mean(CoVList_FD3000)
stdCoV_FD3000 = np.std(CoVList_FD3000)
RbarList_FD3000 = [sum(f[5][-fsteps:])/fsteps for f in allResFD3000]
aveRbar_FD3000 = np.mean(RbarList_FD3000)
stdRbar_FD3000 = np.std(RbarList_FD3000)
aveFt_FD3000 = np.mean([f[0] for f in allResFD3000])
stdFt_FD3000 = np.std([f[0] for f in allResFD3000])
averbar_FD3000 = np.mean([f[1] for f in allResFD3000])
stdrbar_FD3000 = np.std([f[1] for f in allResFD3000])
imFD3000 = allResFD3000[99][6]
    
loadpath = "Path/To/File/out_FD_5000stones"
fileFD5000 = open(loadpath, 'rb')
allResFD5000 = pickle.load(fileFD5000)
fileFD5000.close()
CoVList_FD5000 = [sum(f[4][-fsteps:])/fsteps for f in allResFD5000]
aveCoV_FD5000 = np.mean(CoVList_FD5000)
stdCoV_FD5000 = np.std(CoVList_FD5000)
RbarList_FD5000 = [sum(f[5][-fsteps:])/fsteps for f in allResFD5000]
aveRbar_FD5000 = np.mean(RbarList_FD5000)
stdRbar_FD5000 = np.std(RbarList_FD5000)
aveFt_FD5000 = np.mean([f[0] for f in allResFD5000])
stdFt_FD5000 = np.std([f[0] for f in allResFD5000])
averbar_FD5000 = np.mean([f[1] for f in allResFD5000])
stdrbar_FD5000 = np.std([f[1] for f in allResFD5000])
imFD5000 = allResFD5000[99][6]

## T files
loadpath = "Path/To/File/out_T_1000stones"
fileT1000 = open(loadpath, 'rb')
allResT1000 = pickle.load(fileT1000)
fileT1000.close()
CoVList_T1000 = [sum(f[4][-fsteps:])/fsteps for f in allResT1000]
aveCoV_T1000 = np.mean(CoVList_T1000)
stdCoV_T1000 = np.std(CoVList_T1000)
RbarList_T1000 = [sum(f[5][-fsteps:])/fsteps for f in allResT1000]
aveRbar_T1000 = np.mean(RbarList_T1000)
stdRbar_T1000 = np.std(RbarList_T1000)
aveFt_T1000 = np.mean([f[0] for f in allResT1000])
stdFt_T1000 = np.std([f[0] for f in allResT1000])
averbar_T1000 = np.mean([f[1] for f in allResT1000])
stdrbar_T1000 = np.std([f[1] for f in allResT1000])
imT1000 = allResT1000[99][6]

loadpath = "Path/To/File/out_T_3000stones"
fileT3000 = open(loadpath, 'rb')
allResT3000 = pickle.load(fileT3000)
fileT3000.close()
CoVList_T3000 = [sum(f[4][-fsteps:])/fsteps for f in allResT3000]
aveCoV_T3000 = np.mean(CoVList_T1000)
stdCoV_T3000 = np.std(CoVList_T1000)
RbarList_T3000 = [sum(f[5][-fsteps:])/fsteps for f in allResT3000]
aveRbar_T3000 = np.mean(RbarList_T3000)
stdRbar_T3000 = np.std(RbarList_T3000)
aveFt_T3000 = np.mean([f[0] for f in allResT3000])
stdFt_T3000 = np.std([f[0] for f in allResT3000])
averbar_T3000 = np.mean([f[1] for f in allResT3000])
stdrbar_T3000 = np.std([f[1] for f in allResT3000])
imT3000 = allResT3000[99][6]

loadpath = "Path/To/File/out_T_5000stones"
fileT5000 = open(loadpath, 'rb')
allResT5000 = pickle.load(fileT5000)
fileT5000.close()
CoVList_T5000 = [sum(f[4][-fsteps:])/fsteps for f in allResT5000]
aveCoV_T5000 = np.mean(CoVList_T5000)
stdCoV_T5000 = np.std(CoVList_T5000)
RbarList_T5000 = [sum(f[5][-fsteps:])/fsteps for f in allResT5000]
aveRbar_T5000 = np.mean(RbarList_T5000)
stdRbar_T5000 = np.std(RbarList_T5000)
aveFt_T5000 = np.mean([f[0] for f in allResT5000])
stdFt_T5000 = np.std([f[0] for f in allResT5000])
averbar_T5000 = np.mean([f[1] for f in allResT5000])
stdrbar_T5000 = np.std([f[1] for f in allResT5000])
imT5000 = allResT5000[99][6]

    
savepath = "Path/To/File/stone_n_results.csv"
with open(savepath, mode = "w", newline='') as stat_file:
    stat_writer = csv.writer(stat_file, delimiter = ',', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
    
    stat_writer.writerow(["model","stone n", "mean_CoV", "SD_CoV", "mean_Rbar", "SD_Rbar", "mean_Ft", "SD_Ft", "mean_rbar", "SD_rbar"])
    stat_writer.writerow(["SO", "1000", aveCoV_FD1000, stdCoV_FD1000, aveRbar_FD1000, stdRbar_FD1000, aveFt_FD1000, stdFt_FD1000, averbar_FD1000, stdrbar_FD1000])
    stat_writer.writerow(["SO", "3000", aveCoV_FD3000, stdCoV_FD3000, aveRbar_FD3000, stdRbar_FD3000, aveFt_FD3000, stdFt_FD3000, averbar_FD3000, stdrbar_FD3000])
    stat_writer.writerow(["SO", "5000", aveCoV_FD5000, stdCoV_FD5000, aveRbar_FD5000, stdRbar_FD5000, aveFt_FD5000, stdFt_FD5000, averbar_FD5000, stdrbar_FD5000])
    stat_writer.writerow(["T", "1000", aveCoV_T1000, stdCoV_T1000, aveRbar_T1000, stdRbar_T1000, aveFt_T1000, stdFt_T1000, averbar_T1000, stdrbar_T1000])
    stat_writer.writerow(["T", "3000", aveCoV_T3000, stdCoV_T3000, aveRbar_T3000, stdRbar_T3000, aveFt_T3000, stdFt_T3000, averbar_T3000, stdrbar_T3000])
    stat_writer.writerow(["T", "5000", aveCoV_T5000, stdCoV_T5000, aveRbar_T5000, stdRbar_T5000, aveFt_T5000, stdFt_T5000, averbar_T5000, stdrbar_T5000])

listOfLists = [CoVList_FD1000, CoVList_FD3000, CoVList_FD5000, RbarList_FD1000, RbarList_FD3000, RbarList_FD5000, CoVList_T1000, CoVList_T3000, CoVList_T5000, RbarList_T1000, RbarList_T3000, RbarList_T5000]
savepathL = "Path/To/File/Rbar_and_CoV_values.csv"
with open(savepathL, mode = "w", newline='') as lists_file:
    list_writer = csv.writer(lists_file, delimiter = ',', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
    
    list_writer.writerow(["SO1000_CoV", "SO3000_CoV", "SO5000_CoV", "SO1000_Rbar", "SO3000_Rbar", "SO5000_Rbar", "T1000_CoV", "T3000_CoV", "T5000_CoV", "T1000_Rbar", "T3000_Rbar", "T5000_Rbar"])
    for i in range(len(CoVList_FD1000)):
        temp = []
        for myList in listOfLists:
            temp.append(myList[i])
        list_writer.writerow(temp)
    

# delayed colony size increase
import csv
import numpy as np

stones = 3000
S = 100

outList5000 = list()
outList6000 = list()
outList7500 = list()
outList10000 = list()
outList20000 = list()
outList50000 = list()
outList75000 = list()
for s in range(S):
   out = w.FD_delayed(s,0,18, stones = stones)
   outList5000.append(out[0])
   outList6000.append(out[1])
   outList7500.append(out[2])
   outList10000.append(out[3])
   outList20000.append(out[4])
   outList50000.append(out[5])
   outList75000.append(out[6])
   print("FD_delayed sim is at:" + str(s))
   
filepath5000 = "Path/To/File/out_FD_delayed24_T5000" + str(stones) +"stones"
filepath6000 = "Path/To/File/out_FD_delayed24_T6000" + str(stones) +"stones"
filepath7500 = "Path/To/File/out_FD_delayed24_T7500" + str(stones) +"stones"
filepath10000 = "Path/To/File/out_FD_delayed24_T10000" + str(stones) +"stones"
filepath20000 = "Path/To/File/out_FD_delayed24_T20000" + str(stones) +"stones"
filepath50000 = "Path/To/File/out_FD_delayed24_T50000" + str(stones) +"stones"
filepath75000 = "Path/To/File/out_FD_delayed24_T75000" + str(stones) +"stones"

file5000 = open(filepath5000, 'wb')
pickle.dump(outList5000, file5000)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
file5000.close()

file6000 = open(filepath6000, 'wb')
pickle.dump(outList6000, file6000)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
file6000.close()

file7500 = open(filepath7500, 'wb')
pickle.dump(outList7500, file7500)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
file7500.close()

file10000 = open(filepath10000, 'wb')
pickle.dump(outList10000, file10000)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
file10000.close()

file20000 = open(filepath20000, 'wb')
pickle.dump(outList20000, file20000)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
file20000.close()

file50000 = open(filepath50000, 'wb')
pickle.dump(outList50000, file50000)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
file50000.close()

file75000 = open(filepath75000, 'wb')
pickle.dump(outList75000, file75000)  # each pickled file only contains one object (i.e. if multiple objects need to be saved, they need to be structured as subelements of a bigger objects by the user)
file75000.close()


fsteps = 10
loadpath1 = "Path/To/File/out_FD_delayed24_T50003000stones"
filed24_5000 = open(loadpath1, 'rb')
allResd24_5000 = pickle.load(filed24_5000)
filed24_5000.close()
rbarListd24_5000 = [sum(f[3][-fsteps:])/fsteps for f in allResd24_5000]
imd24_5000 = allResd24_5000[99][6]

fsteps = 10
loadpath2 = "Path/To/File/out_FD_delayed24_T60003000stones"
filed24_6000 = open(loadpath2, 'rb')
allResd24_6000 = pickle.load(filed24_6000)
filed24_6000.close()
rbarListd24_6000 = [sum(f[3][-fsteps:])/fsteps for f in allResd24_6000]
imd24_6000 = allResd24_6000[99][6]

fsteps = 10
loadpath3 = "Path/To/File/out_FD_delayed24_T75003000stones"
filed24_7500 = open(loadpath3, 'rb')
allResd24_7500 = pickle.load(filed24_7500)
filed24_7500.close()
rbarListd24_7500 = [sum(f[3][-fsteps:])/fsteps for f in allResd24_7500]
imd24_7500 = allResd24_7500[99][6]

loadpath4 = "Path/To/File/out_FD_delayed24_T100003000stones"
filed24_10000 = open(loadpath4, 'rb')
allResd24_10000 = pickle.load(filed24_10000)
filed24_10000.close()
rbarListd24_10000 = [sum(f[3][-fsteps:])/fsteps for f in allResd24_10000]
imd24_10000 = allResd24_10000[99][6]

loadpath5 = "Path/To/File/out_FD_delayed24_T200003000stones"
filed24_20000 = open(loadpath5, 'rb')
allResd24_20000 = pickle.load(filed24_20000)
filed24_20000.close()
rbarListd24_20000 = [sum(f[3][-fsteps:])/fsteps for f in allResd24_20000]
imd24_20000 = allResd24_20000[99][6]

loadpath6 = "Path/To/File/out_FD_delayed24_T500003000stones"
filed24_50000 = open(loadpath6, 'rb')
allResd24_50000 = pickle.load(filed24_50000)
filed24_50000.close()
rbarListd24_50000 = [sum(f[3][-fsteps:])/fsteps for f in allResd24_50000]
imd24_50000 = allResd24_50000[99][6]

loadpath7 = "Path/To/File/out_FD_delayed24_T750003000stones"
filed24_75000 = open(loadpath7, 'rb')
allResd24_75000 = pickle.load(filed24_75000)
filed24_75000.close()
rbarListd24_75000 = [sum(f[3][-fsteps:])/fsteps for f in allResd24_75000]
imd24_75000 = allResd24_75000[99][6]

listOfDelayedLists = [rbarListd24_5000, rbarListd24_6000, rbarListd24_7500, rbarListd24_10000, rbarListd24_20000, rbarListd24_50000, rbarListd24_75000]

savepathd = "Path/To/File/Rbar_and_CoV_values_delayed.csv"
with open(savepathd, mode = "w", newline='') as listsd_file:
    list_writer = csv.writer(listsd_file, delimiter = ',', quotechar = '"', quoting = csv.QUOTE_MINIMAL)
    
    list_writer.writerow(["SO_T5000", "SO_T6000", "SO_T7500", "SO_T10000", "SO_T20000", "SO_T50000", "SO_T75000"])
    for i in range(len(rbarListd24_10000)):
        temp = []
        for myList in listOfDelayedLists:
            temp.append(myList[i])
        list_writer.writerow(temp)



# images
f, axarr = plt.subplots(2,3)
axarr[0,0].imshow(imFD1000, cmap = 'gray')
axarr[0,0].xticks([])
axarr[0,1].imshow(imFD3000, cmap = 'gray')
axarr[0,2].imshow(imFD5000, cmap = 'gray')
axarr[1,0].imshow(imT1000, cmap = 'gray')
axarr[1,1].imshow(imT3000, cmap = 'gray')
axarr[1,2].imshow(imT5000, cmap = 'gray')

_, axs = plt.subplots(2, 3)
imgs = [imFD1000, imFD3000, imFD5000, imT1000, imT3000, imT5000]
axs = axs.flatten()
nms = ["stones = 1000", "stones = 3000", "stones = 5000", "", "", ""]
for img, ax, nm in zip(imgs, axs, nms):
    ax.imshow(img, cmap = 'gray')
    ax.axis('off')
plt.savefig('Path/To/File/walls.png', dpi=350)
plt.show()

# for single images
img = imFD1000
im = plt.imshow(img, cmap = 'gray')
plt.colorbar()
plt.axis('off')


# delayed colony size increase - images
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('xtick', labelsize=9) 
matplotlib.rc('ytick', labelsize=9) 

_, axs = plt.subplots(1, 3)
imgs = [imd24_10000, imd24_20000, imd24_75000]
axs = axs.flatten()
nms = ["T = 6000", "T = 10000", "stones = 75000"]
for img, ax, nm in zip(imgs, axs, nms):
    ax.imshow(img, cmap = 'gray')
    ax.xaxis.set_ticks(np.arange(0,81,20))
    ax.yaxis.set_ticks(np.arange(0,81,20))
    if (img==imgs[0]).all() == False:
        ax.get_yaxis().set_visible(False)
plt.show()
plt.savefig('Path/To/File/walls_delayed24.png', dpi=350)

