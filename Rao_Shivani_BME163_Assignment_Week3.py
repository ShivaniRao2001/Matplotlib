import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse
import matplotlib.patheffects as pe

parser = argparse.ArgumentParser()
parser.add_argument('--inFile', '-i', type=str, action='store', help='input file')
parser.add_argument('--outFile', '-o', type=str, action='store', help='output file')
args = parser.parse_args()
inFile = args.inFile
outFile = args.outFile
plt.style.use('BME163.mplstyle')

FileCellType = "BME163_Input_Data_Week3.celltype.tsv"
FilePosition = "BME163_Input_Data_Week3.position.tsv"

figureWidth = 8
figureHeight = 4
panelWidth = 2
panelHeight = 2
relativePanelWidth = panelWidth / figureWidth
relativePanelHeight = panelHeight / figureHeight
plt.figure(figsize=(figureWidth, figureHeight),linewidth = 3, edgecolor="black")   #figure out o

Panel = plt.axes([0.14, 0.15, relativePanelWidth, relativePanelHeight])
Panel.set_xticks([-20, 0, 20],fontsize = 11.5)
Panel.set_yticks([-40,-30,-20,-10,0,10,20,30])

Panel.set_xlim(-30, 30)
Panel.set_ylim(-40, 30)

Panel.set_xlabel("tSNE 2", fontsize = 11)
Panel.set_ylabel("tSNE 1",fontsize = 11)

seqs = []

for line in open(FileCellType):
    a = line.rstrip().split('\t')
    seq = a[2]
    seqs.append(seq)


xcoords = []  # all coordinates
ycoords = []

for line in open(FilePosition):
    b = line.split()
    x = float(b[1])
    y = float(b[2])
    xcoords.append(x)
    ycoords.append(y)

ys = []  # red
xs = []
indicesred = []
for each in ycoords:
    if each >= -10:
        ys.append(each)
for each in ys:
    indicesred.append(ycoords.index(each))

for each in indicesred:
    xs.append(xcoords[each])

Panel.scatter(xs, ys, s=8.5, c="lightcoral")

restx = []
resty = []
for each in xcoords:
    if each not in xs:
        restx.append(each)
for each in ycoords:
    if each not in ys:
        resty.append(each)
greyx = []
greyy = []
indicesgrey = []

for each in restx:
    if -20 < each < 15:
        greyx.append(each)
for each in greyx:
    indicesgrey.append(restx.index(each))
for each in indicesgrey:
    greyy.append(resty[each])

Panel.scatter(greyx, greyy, c="peru", s=8.5)

Panel.scatter(-11.72618866,-18.19581099, c = "lightcoral", s=8.5)

bluex = []
bluey = []
for each in restx:
    if each not in greyx:
        bluex.append(each)
for each in resty:
    if each not in greyy:
        bluey.append(each)

Panel.scatter(bluex, bluey, c="powderblue", s=8.5)

Panel.text(-0.31392918, 7.88747802, "tCell", va="center", ha = "center",
           path_effects=[pe.withStroke(linewidth=1, foreground="white")])
Panel.text(-4.680359303, -26.60389215, "monocyte", va = "center",ha = "center",
           path_effects=[pe.withStroke(linewidth=1, foreground="white")])
Panel.text(23.02149858, -18.66172484, "bCell", va="center", ha ="center",
           path_effects=[pe.withStroke(linewidth=1, foreground="white")])


plt.savefig(outFile, dpi=600)
