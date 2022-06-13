import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse
import matplotlib.image as mpimg
import math

# note: my computer was struggling with the gz file. It kept crashing everytime it was close to being downloaded so I used the chr15 file for my code instead.

parser = argparse.ArgumentParser()
parser.add_argument('--inFile', '-g', type=str, action='store', help='input file')
parser.add_argument('--bedFile', '-b', type=str, action='store', help='bed file', default='Splice_Locations_chr15.bed')
parser.add_argument('--aFile', '-A', type=str, action='store', help='input file', default='A.png')
parser.add_argument('--cFile', '-C', type=str, action='store', help='input file', default='C.png')
parser.add_argument('--tFile', '-T', type=str, action='store', help='input file', default='T.png')
parser.add_argument('--gFile', '-G', type=str, action='store', help='input file', default='G.png')
parser.add_argument('--outFile', '-o', type=str, action='store', help='output file')
args = parser.parse_args()
inFile = args.inFile
A = mpimg.imread(args.aFile)
C = mpimg.imread(args.cFile)
T = mpimg.imread(args.tFile)
G = mpimg.imread(args.gFile)
outFile = args.outFile
bedFile = args.bedFile
plt.style.use('BME163.mplstyle')


def readFasta(genomeFile):
    fastaDict = {}
    first = False
    for line in open(genomeFile):
        line = line.rstrip()
        if line[0] == '>':
            if first:
                fastaDict[name]=('').join(fastaDict[name])
            name = "chr" + str(line[1:].split(" ",1)[0])
            fastaDict[name] = []
            first = True
        else:
            sequence = line
            fastaDict[name].append(sequence)
    fastaDict[name] = ('').join(fastaDict[name])
    return fastaDict

fastaDict = readFasta(inFile)

junctions = {3: [], 5: []}
sequences = {3: [], 5: []}

def readBed(bedFile):
    dict3, dict5 = {}, {}
    with open(bedFile) as bed:
        for line in bed:
            elements = line.split("\t")
            if elements[3][0] == "3":
                if elements[0] not in dict3:
                    dict3[elements[0]] = []
                dict3[elements[0]].append(elements[1])
            else:
                if elements[0] not in dict5:
                    dict5[elements[0]] = []
                dict5[elements[0]].append(elements[1])

    return dict3, dict5

dict3, dict5 = readBed(bedFile)

sequences5 = []
sequences3 = []
for key, value in dict5.items():
    for point in value:
        sequences5.append(fastaDict[key][int(point)-10:int(point)+10])

for key, value in dict3.items():
    for point in value:
        sequences3.append(fastaDict[key][int(point)-10:int(point)+10])

reverses = []  # 3 prime reverse compliments.
for thisseq in sequences3:  # reverse compliment all 3 prime sequences
    reverse = thisseq[::-1]  # reverse the strand
    reversecompliment = ""
    for letter in reverse:
        if letter == "A":
            reversecompliment += "T"
        elif letter == "T":
            reversecompliment += "A"
        elif letter == "G":
            reversecompliment += "C"
        elif letter == "C":
            reversecompliment += "G"
    reverses.append(reversecompliment)

def frequency(lists):
    dict = {}    #divide by the total
    baselist = ['A', 'G', 'T', 'C']
    for index in range(20):
        dict[index] = {'A': 0, 'G': 0, 'T': 0, 'C': 0}
    for string in lists:
        for n, i in enumerate(string):
            # n index , i letter
            for base in baselist:
                if base == i:
                    dict[n][i] = int(dict[n][i]) + 1

    tot = dict.get(0)
    sumtot = sum(list(tot.values()))
    for key, value in dict.items():
        for k, v in value.items():
            value.update({k: v / sumtot})
    return dict

def height(dict):
    heightfreq = []
    for key, value in dict.items():
        afreq = value.get("A")
        gfreq = value.get("G")
        cfreq = value.get("C")
        tfreq = value.get("T")
        h_a = afreq * (np.log2(afreq))
        h_g = gfreq * (np.log2(gfreq))
        h_c = cfreq * (np.log2(cfreq))
        h_t = tfreq * (np.log2(tfreq))
        total_h = -(h_a + h_g + h_c + h_t)  #for fives, should be 0.16
        bits = 2 - (total_h)   #-1*4*(0.25*np.log2(0.25))
        aheight = afreq * bits
        gheight = gfreq * bits
        cheight = cfreq * bits
        theight = tfreq * bits
        heightfreq.append([(A,aheight),(T,theight),(G,gheight),(C,cheight)])
        # heightfreq[key] = {"A": [afreq, aheight], "G": [gfreq, gheight], "C": [cfreq, cheight], "T": [tfreq, theight]}

    return heightfreq

frequenciesfive = frequency(sequences5)  # dictionary keys ranging from 1-20, contains dictionary as values. with nuc:freq as items. .
frequenciesthree = frequency(reverses)

heightfive = height(frequenciesfive)
heightthree = height(frequenciesthree)

figureWidth = 5
figureHeight = 2
panelWidth = 1.5   #was 1.5
panelHeight = 0.75   #was 0.5

relativePanelWidth = panelWidth / figureWidth
relativePanelHeight = panelHeight / figureHeight

plt.figure(figsize=(figureWidth, figureHeight))

Panel1 = plt.axes([0.1, 0.3, relativePanelWidth, relativePanelHeight])
Panel2 = plt.axes([0.44, 0.3, relativePanelWidth, relativePanelHeight])

Panel1.set_yticks([0, 1, 2])
Panel1.set_xticks([-10, -5, 0, 5, 10])
Panel1.set_xlabel("Distance to\nSplice Site", multialignment="center")
Panel2.set_xlabel("Distance to\nSplice Site", multialignment="center")
Panel1.set_ylabel("Bits")
Panel1.set_title("5'SS")
Panel2.set_title("3'SS")
Panel2.set_xticks([-10, -5, 0, 5, 10])

Panel1.plot([0,0],[0,2], "k-",lw = 0.6)
Panel2.plot([0,0],[0,2], "k-",lw = 0.6)

Panel2.tick_params(bottom=True, labelbottom=True,
                   left=False, labelleft=False,
                   right=False, labelright=False,
                   top=False, labeltop=False)
Panel1.tick_params(bottom=True, labelbottom=True,
                   left=True, labelleft=True,
                   right=False, labelright=False,
                   top=False, labeltop=False)

x = -10
for tuple in heightfive:
    prevHeight = 0
    for base in sorted(tuple, key = lambda x:x[1]):
        Panel1.imshow(base[0], extent=[x, x + 1, prevHeight, prevHeight + base[1]],aspect='auto')
        prevHeight += base[1]
    x += 1

x2 = -10
for tuple in heightthree:
    prevHeight = 0
    for base in sorted(tuple, key = lambda x2:x2[1]):
        Panel2.imshow(base[0], extent=[x2, x2 + 1, prevHeight, prevHeight + base[1]],aspect='auto')
        prevHeight += base[1]
    x2 += 1

Panel1.set_xlim(-10, 10)
Panel2.set_xlim(-10, 10)
Panel2.set_ylim(0, 2)
Panel1.set_ylim(0, 2)

plt.savefig(outFile, dpi=600)