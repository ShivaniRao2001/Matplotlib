import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import math
import argparse

fileName= 'BME163_Input_Data_1.txt'
plt.style.use('BME163.mplstyle')

yvals = [0,5,10,15]
xvals = [0,5,10,15]

parser = argparse.ArgumentParser()
parser.add_argument('--inFile', '-i', type=str, action='store', help='input file')
parser.add_argument('--outFile', '-o', type=str, action='store', help='output file')
args = parser.parse_args()
inFile = args.inFile
outFile = args.outFile

xmin = 0
xmax = 15
ymin = 0
ymax = 20

figureWidth= 5
figureHeight= 2

panelWidth = 1
panelHeight = 1

sidePanelWidth = 0.25
sidePanelHeight = 1

topPanelWidth = 1
topPanelHeight = 0.25

relativePanelWidth = panelWidth/figureWidth
relativePanelHeight = panelHeight/figureHeight
relativeSideWidth = sidePanelWidth/figureWidth
relativeSideHeight = sidePanelHeight/figureHeight
relativeTopWidth = topPanelWidth/figureWidth
relativeTopHeight = topPanelHeight/figureHeight

plt.figure(figsize=(figureWidth, figureHeight))

MainPanel = plt.axes([0.14,0.25,relativePanelWidth,relativePanelHeight])   #(left,bottom,width,height)
SidePanel = plt.axes([0.076,0.25,relativeSideWidth,relativeSideHeight])
TopPanel = plt.axes([0.14,0.785,relativeTopWidth,relativeTopHeight])

MainPanel.set_xlim(0,15)
MainPanel.set_ylim(0,15)
SidePanel.set_xlim(20,0)
SidePanel.set_ylim(0,15)
TopPanel.set_ylim(0,20)
TopPanel.set_xlim(0,15)

TopPanel.set_ylabel("# of\ngenes")
SidePanel.set_ylabel("Replicate 2 (RPM)")
SidePanel.set_xlabel("# of\ngenes")
MainPanel.set_xlabel("Replicate 1 (RPM)")

MainPanel.tick_params(bottom=True, labelbottom=True,
                   left=False, labelleft=False,
                   right=False, labelright=False,
                   top=False, labeltop=False)

MainPanel.set_xticks(xvals)
SidePanel.set_yticks(yvals)
SidePanel.set_xticks([0,20])
SidePanel.set_xticks([0,20])


SidePanel.tick_params(bottom=True, labelbottom=True,
                   left=True, labelleft=True,
                   right=False, labelright=False,
                   top=False, labeltop=False)

TopPanel.tick_params(bottom=False, labelbottom=False,
                   left=True, labelleft=True,
                   right=False, labelright=False,
                   top=False, labeltop=False)

xList = []
yList = []

for line in open(fileName):
    a=line.rstrip().split('\t')
    xs=float(a[1])
    ys = float(a[2])
    xscat = math.log2(xs + 1)    #converted value for scatter plot
    yscat = math.log2(ys+1)   #converted value for scatter plot
    xList.append(xscat)
    yList.append(yscat)   #converted height of the histogram


MainPanel.plot(xList, yList,
               marker='o',
               markersize=2,
               linewidth=0,
               linestyle='--',
               markeredgewidth=0,
               markerfacecolor='black',
               alpha=0.1
               )

bins=np.arange(0,15,0.5)
xHisto,bins=np.histogram(xList,bins)

for index in range(0,len(xHisto),1):
    bottom=0
    left=bins[index]
    width=bins[index+1]-left
    height=np.log2(xHisto[index]+1)
    rectangle1=mplpatches.Rectangle([left,bottom],width,height,
                                 linewidth=0.1,
                                 facecolor="grey",
                                edgecolor = "black"
                                 )
    TopPanel.add_patch(rectangle1)

yHisto,bins = np.histogram(yList,bins)


for index in range(len(yHisto)):
    left =0
    bottom = bins[index]
    height = bins[index+1]-bottom
    width = np.log2(yHisto[index]+1)
    rectangle2 = mplpatches.Rectangle([left,bottom],width,height,
                                      linewidth = 0.1,
                                      facecolor = "grey",
                                      edgecolor = "black")
    SidePanel.add_patch(rectangle2)


plt.savefig(outFile, dpi=600)


















