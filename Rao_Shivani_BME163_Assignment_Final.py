import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse
import matplotlib.image as mpimg
import math

#py Rao_Shivani_BME163_Assignment_Final.py -i1 BME163_Input_Data_5.psl -i2 BME163_Input_Data_6.psl -g gencode.vM12.annotation.gtf -c chr7:45232000-45241000 -o Rao_Shivani_BME163_Assignment_Final.png

parser = argparse.ArgumentParser()
parser.add_argument('--pslFile', '-i1', type=str, action='store', help='input file', default='BME163_Input_Data_5.psl')
parser.add_argument('--psl2File', '-i2', type=str, action='store', help='input file', default='BME163_Input_Data_6.psl')
parser.add_argument('--gtfFile', '-g', type=str, action='store', help='input file', default='gencode.vM12.annotation.gtf')
parser.add_argument('--outFile', '-o', type=str, action='store', help='output file')
parser.add_argument('--chromFile', '-c', type=str, action='store', help='output file')
args = parser.parse_args()
pslFile = args.pslFile
psl2File = args.psl2File
outFile = args.outFile
gtfFile = args.gtfFile
chromFile = args.chromFile

plt.style.use('BME163.mplstyle')

figureWidth = 5
figureHeight = 5
panelWidth = 4.85
panelHeight = 1.5

relativePanelWidth = panelWidth / figureWidth
relativePanelHeight = panelHeight / figureHeight

plt.figure(figsize=(figureWidth, figureHeight))

Panel3 = plt.axes([0.015, 0.01, relativePanelWidth, relativePanelHeight])
Panel2 = plt.axes([0.015, 0.33, relativePanelWidth, relativePanelHeight])
Panel1 = plt.axes([0.015, 0.65, relativePanelWidth, relativePanelHeight])


Panel1.tick_params(bottom=False, labelbottom=False, #gtf
                   left=False, labelleft=False,
                   right=False, labelright=False,
                   top=False, labeltop=False)

Panel2.tick_params(bottom=False, labelbottom=False,
                   left=False, labelleft=False,
                   right=False, labelright=False,
                   top=False, labeltop=False)

Panel3.tick_params(bottom=False, labelbottom=False,
                   left=False, labelleft=False,
                   right=False, labelright=False,
                   top=False, labeltop=False)

iYellow=(248/255,174/255,51/255)
iBlue=(88/255,85/255,120/255)

def readGTF(gtf,chromo,areaS,areaE):

    dict = {}
    transcriptlist = []
    for line in gtf:
        split = line.split()
        # split2  = line.strip().split('\t')
        chr = split[0]
        chrm = chr[3]
        start = int(split[3])
        end = int(split[4])
        typegtf = split[2]
        transcript_id = split[11]

        # transcript = split2[8].split(' transcript_id "')[0].split('"')[1]

        if chrm == chromo:
            if areaS<start<areaE or areaS<end<areaE:
                if typegtf == "CDS" or typegtf == "exon":
                    if transcript_id not in dict.keys():
                        dict[transcript_id] = []
                    dict[transcript_id].append([chrm, start, end, typegtf])

    for key,val in dict.items():
        starts = []
        ends = []
        blockstarts = []
        blockwidths = []
        types = []
        placed = False
        for each in val:
            chr = each[0]
            starts.append(each[1])
            ends.append(each[2])
            blockstarts.append(each[1])
            width = each[2]-each[1]
            blockwidths.append(width)
            types.append(each[3])
        transcriptlist.append([chr,min(starts),max(ends),blockstarts,blockwidths,types,placed])

    return transcriptlist


def readPSL(psl,chromo,areaS,areaE):

    newlist = []
    for line in psl:
        split = line.split()
        start = int(split[15])
        end = int(split[16])
        chr = split[13]
        chrm = chr[3]
        placed = False
        if chrm == chromo:
            if areaE>start >areaS or areaS<end<areaE:
                blockwidths = np.array(split[18].split(',')[:-1], dtype=int)
                blockstarts = np.array(split[20].split(",")[:-1], dtype = int)
                read = [chrm,start,end,blockstarts,blockwidths,placed]
                newlist.append(read)

    return newlist

def plotPSLAlignments(alignment,type):

    ypos = []
    if type == 5:
        color = iYellow
        alignment = sorted(alignment, key=lambda i: i[2], reverse=True)
    elif type == 6:
        color = iBlue
        alignment = sorted(alignment, key=lambda i: i[1], reverse=False)

    for y_pos in np.arange(1,len(alignment)+0.75,1):
        last_plotted = -1
        for each in alignment:
            chrom = each[0]
            start = each[1]
            end = each[2]
            blockstarts = each[3]
            width = each[4]
            placed = each[5]
            ypos.append(y_pos)
            if not placed and start > last_plotted:
                last_plotted = end
                rectangle = mplpatches.Rectangle((start,float(y_pos + 0.2)), end-start, 0.1, facecolor =color)   #blockwidth = end-start
                if type == 5:
                    Panel3.add_patch(rectangle)
                elif type == 6:
                    Panel2.add_patch(rectangle)
                for i in np.arange(0,len(blockstarts),1):
                    blockstart = int(blockstarts[i])
                    blockwidth = int(width[i])
                    rectangle2 = mplpatches.Rectangle((blockstart, int(y_pos)), blockwidth, 0.5, facecolor=color)
                    if type == 5:
                        Panel3.add_patch(rectangle2)
                    elif type == 6:
                        Panel2.add_patch(rectangle2)

                each[5] = True

    return y_pos


def plotGTF(tlist,chromo,areastart,areaend):

    plotted = []
    finished = False

    for each in tlist:
        placed = each[6]
        plotted.append(placed)

    starting = 0.5
    while not finished:
        finished = True
        last_plotted = 1
        for each in range(0, len(tlist)):
            if plotted[each] == False:
                chromosome1 = tlist[each]
                chromosome = chromosome1[0]
                start1 = tlist[each]
                start = start1[1]
                end1 = tlist[each]
                end = end1[2]
                blockstarts1 = tlist[each]
                blockstarts = blockstarts1[3]
                blockwidths1 = tlist[each]
                blockwidths = blockwidths1[4]
                type1 = tlist[each]
                type = type1[5]
                if start > last_plotted:
                    plotted[each] = True
                    rectangle1 = mplpatches.Rectangle((start, starting + 0.2), end - start, 0.1, facecolor='grey',edgecolor='black', linewidth=0.2)
                    Panel1.add_patch(rectangle1)
                    last_plotted = end
                    for y_pos in np.arange(0, len(blockstarts),1):
                        blockstart = blockstarts[y_pos]
                        blockwidth = blockwidths[y_pos]
                        if type[y_pos] == 'exon':
                            rectangle2 = mplpatches.Rectangle((blockstart, starting + 0.04), blockwidth, 0.25, edgecolor='black', linewidth=0.2, facecolor='grey',)
                            Panel1.add_patch(rectangle2)
                        elif type[y_pos] == 'CDS':
                            rectangle3 = mplpatches.Rectangle((blockstart, starting), blockwidth, 0.5, facecolor='grey',edgecolor='black', linewidth=0.2)
                            Panel1.add_patch(rectangle3)
                elif start <= last_plotted:
                    finished = False

        limits = []
        limits.append(starting)
        starting += 1

    return min(limits), max(limits)


openGTF = open(gtfFile)
openPSL1 = open(pslFile)
openPSL2 = open(psl2File)

next(openGTF)
next(openGTF)
next(openGTF)
next(openGTF)
next(openGTF)

chromosome = chromFile[3]
areas = chromFile[5:].split("-")
areaStart = int(areas[0])
areaEnd = int(areas[1])

psl1 = 5
psl2 = 6

PSL1align = readPSL(openPSL1,chromosome,areaStart,areaEnd)  #5 third, yellow
GTFalign = readGTF(openGTF, chromosome, areaStart, areaEnd)
PSL2align = readPSL(openPSL2,chromosome,areaStart,areaEnd)  #6 middle, blue

Panel1.set_xlim(areaStart, areaEnd)
Panel2.set_xlim(areaStart, areaEnd)
Panel3.set_xlim(areaStart, areaEnd)

PSL1plot = plotPSLAlignments(PSL1align, psl1)
PSL2plot = plotPSLAlignments(PSL2align, psl2)
gtfplot = plotGTF(GTFalign,chromosome,areaStart,areaEnd)

min = gtfplot[0]
max = gtfplot[1]

Panel2.set_ylim(0,71)
Panel3.set_ylim(0,450)
Panel1.set_ylim(-0.4,max+2)

plt.savefig(outFile, dpi=2400)