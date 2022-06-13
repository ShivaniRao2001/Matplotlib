import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.ticker as MultipleLocator
import numpy as np
import argparse
import matplotlib.image as mpimg
import math

parser = argparse.ArgumentParser()
parser.add_argument('--inFile', '-i', type=str, action='store', help='input file')
parser.add_argument('--outFile', '-o', type=str, action='store', help='output file')
args = parser.parse_args()
inFile = args.inFile
outFile = args.outFile
plt.style.use('BME163.mplstyle')

input_file = "BME163_Input_Data_4.txt"
fileName = open(input_file)

figureWidth = 5
figureHeight = 3
panelWidth1 = 0.75
panelWidth2 = 2.5
panelHeight = 2.5

plt.figure(figsize=(figureWidth, figureHeight))

relativePanelWidth1 = panelWidth1 / figureWidth
relativePanelWidth2 = panelWidth2 / figureWidth

relativePanelHeight = panelHeight / figureHeight

Panel1 = plt.axes([0.1, 0.1, relativePanelWidth1, relativePanelHeight])
Panel2 = plt.axes([0.44, 0.1, relativePanelWidth2, relativePanelHeight])

Panel1.set_ylabel("Number of genes")
Panel1.set_xlabel("CT")

Panel1.set_yticks([0, 200, 400, 600, 800, 1000, 1200])
Panel1.set_xticks([0, 3, 6, 9, 12, 15, 18, 21])
Panel1.set_xticklabels([0, "", 6, "", 12, "", 18, ""])

Panel1.set_xlim([-1.5, 22.5])
Panel1.set_ylim([0, 1262])

iBlue = (44 / 255, 86 / 255, 134 / 255)
iYellow = (248 / 255, 174 / 255, 51 / 255)
iGreen = (32 / 255, 100 / 255, 113 / 255)

viridis5 = (253/255, 231/255, 37/255)
viridis4 = (94/255, 201/255, 98/255)
viridis3 = (33/255, 145/255, 140/255)
viridis2 = (59/255, 82/255, 139/255)
viridis1 = (68/255, 1/255, 84/255)

color1 = viridis1
color2 = viridis2
color3 = viridis3
color4 = viridis4
color5 = viridis5

R1 = np.linspace(color1[0], color2[0], 25)
G1 = np.linspace(color1[1], color2[1], 25)
B1 = np.linspace(color1[2], color2[2], 25)

R2 = np.linspace(color2[0], color3[0], 25)
G2 = np.linspace(color2[1], color3[1], 25)
B2 = np.linspace(color2[2], color3[2], 25)

R3 = np.linspace(color3[0], color4[0],25)
G3 = np.linspace(color3[1], color4[1],25)
B3 = np.linspace(color3[2], color4[2],25)

R4 = np.linspace(color4[0], color5[0],26)
G4 = np.linspace(color4[1], color5[1],26)
B4 = np.linspace(color4[2], color5[2],26)

R = np.concatenate((R1, R2, R3, R4), axis=None)
G = np.concatenate((G1, G2, G3, G4), axis=None)
B = np.concatenate((B1, B2, B3, B4), axis=None)

allFPKM = []

def normalize(FPKM):
    trueFPKM = np.array(FPKM[0])
    int_array = trueFPKM.astype(int)
    norm = (((int_array - min(int_array))/ (max(int_array) - min(int_array))))*100
    norm_int = norm.astype(int)
    allFPKM.append([norm_int,int(float(FPKM[1]))])

next(fileName)
for line in fileName:
    split = line.split()
    lst = [split[4:12], split[13]]
    ylist = normalize(lst)

allFPKM.sort(reverse=True, key=lambda x: x[1])

for list in allFPKM:
    popped = list.pop()

def moving_average(row, bin):
    smooth_row = []
    for index in range(0, len(row), 1):
        window = row[index:index + bin]
        smoothed = sum(window) / len(window)
        smooth_row.append(int(smoothed))
    return smooth_row

y = 0
for row in allFPKM:
    new = row[0]
    y += 1
    x = -1.5
    new = moving_average(new, 1)
    for value in new:
        # x += 1
        rectangle1 = mplpatches.Rectangle([x, y], 3, 3,
                                          linewidth=0,
                                          facecolor=(R[value], G[value], B[value]))
        Panel1.add_patch(rectangle1)
        x += 3

plt.savefig(outFile, dpi=600)
