import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument('--inFile', '-i', type=str, action='store', help='input file')
parser.add_argument('--outFile', '-o', type=str, action='store', help='output file')
args = parser.parse_args()
inFile = args.inFile
outFile = args.outFile
plt.style.use('BME163.mplstyle')

fileName = "BME163_Input_Data_3.txt"

figureWidth = 11
figureHeight = 4
panelWidth = 6
panelHeight = 2.5  #2
relativePanelWidth = panelWidth / figureWidth
relativePanelHeight = panelHeight / figureHeight
plt.figure(figsize=(figureWidth, figureHeight))

x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
y = [75, 80, 85, 90, 95, 100]
Panel = plt.axes([0.12, 0.15, relativePanelWidth, relativePanelHeight])
Panel.set_xticks(x, labels=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ">10"], fontsize = 16)
Panel.set_yticks(y, labels=[75, 80, 85, 90, 95, 100], fontsize = 16)

Panel.set_xlim(0, 12)
Panel.set_ylim(75, 100)

plt.ylabel("Identity (%)",fontsize = 15)
plt.xlabel("Subread coverage",fontsize = 15)

dictreads = {}
keyvals = []

for line in open(fileName):
    a = line.split()
    del a[2:5]
    sub = a[0]
    split = sub.split("_")
    keyvals.append((float(split[3]), float(a[1])))

dictreads[11.0] = []

for each in keyvals:
    key = each[0]
    val = each[1]
    if key not in dictreads and key <=10:
        dictreads[key] = [val]
    elif key not in dictreads and key >10:
        dictreads[11.0].append(val)
    elif key in dictreads.keys():
        dictreads[key].append(val)


def swarmplot(values,key,Panel):  # values obtained from dictionary "dictreads"  #obtained one by one. so an each == an item == (key: value pair), ex 90: [1,2,3,4]
    pointsize = 1
    minDist = pointsize/72
    increment = minDist/5
    panelWidth = 6
    panelHeight = 2
    span = 0.8
    placedpoints = []
    for yval in values:
        placed = False
        xcoord = int(key)
        if len(placedpoints)==0:
            placedpoints.append((xcoord,yval))
            placed = True
        else:
            for shift in np.arange(0,span/2,increment):
                for val in [-1,1]:
                    if placed == False:
                        shift*=val
                        xcoordshift = float(key) + (shift*val)
                        distList = []
                        for point in placedpoints:
                            x = point[0]
                            y = point[1]
                            xdist = (((xcoordshift - x)/12)*panelWidth)
                            ydist = (((yval - y)/25)*panelHeight)
                            dist = (xdist**2 + ydist**2)**.5
                            distList.append(dist)
                        if min(distList)>minDist:
                            placedpoints.append((xcoordshift,yval))
                            placed = True

    for points in placedpoints:
        x=points[0]
        y=points[1]
        Panel.plot(x,y,
           marker='o',
           markersize=pointsize,
           linewidth=0,
           markeredgewidth=0,
           markerfacecolor='black')

for bins in dictreads:
    accuracy=np.random.choice(dictreads[bins],1000)
    swarmplot(accuracy,bins,Panel)
    key = bins
    pos = np.median(accuracy)
    Panel.plot([key - 0.4, key + 0.4], [pos, pos], color='red', lw=1, zorder=100)


plt.savefig(outFile, dpi=600)
