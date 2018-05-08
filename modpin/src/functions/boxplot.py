
import sys
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import string
import argparse
from scipy import stats

def parse_user_arguments(*args, **kwds):
 parser = argparse.ArgumentParser(
  description = "Create a BOXPLOT with a list of data files",
  epilog      = "@oliva's lab 2014")
 parser.add_argument('-l','--list_input_file',dest='input',action = 'store',default='input_list',
  help = 'Input file (default is input_list)')
 parser.add_argument('-f','--format_image',dest='format',action = 'store',default='png',
  help = 'Format of the image in the output file: png(default), pdf, ps, eps and svg ')
 parser.add_argument('-o','--output_image',dest='out',action = 'store',default='output',
  help = 'Output file name for the imag boxplot')
 parser.add_argument('-top','--Top_Y_value',dest='top',action = 'store',default=0,type=float,
  help = 'Top value of Y-axis')
 parser.add_argument('-bottom','--Botom_Y_value',dest='bottom',action = 'store',default=0,type=float,
  help = 'Bottom value of Y-axis')
 parser.add_argument('-s','--show_canvas',dest='canvas',action = 'store_true',
  help = 'Show canvas figure after execution (not shown by default)')
 parser.add_argument('-g','--show_grid',dest='grid',action = 'store_true',
  help = 'Add horizontal gridlines (default is without)')
 parser.add_argument('-legend','--show_legend',dest='legend',action = 'store_true',
  help = 'Add legends (default is without)')
 parser.add_argument('-t','--title_of_plot',dest='title',action = 'store',default='Boxplot',
  help = 'Title of the plot within "" ')
 parser.add_argument('-nc','--set_colored_groups',dest='color',action = 'store',type=int, default=1,
  help = 'Number of sets grouped with the same color. Up to 14 different colors (default 1 => every boxplot has different color)')
 parser.add_argument('-nh','--set_hashed_groups',dest='hashed',action = 'store',type=int, default=1,
  help = 'Number of sets grouped with the same filling hash. Up to 10 different filling hashes (default 1 => every boxplot has different hash)')
 parser.add_argument('-v','--verbose',dest='show',action = 'store_true',
  help = 'Verbose execution')

 options=parser.parse_args()
 return options

def fileExist (file):
 return os.path.exists(file) and os.path.isfile(file)

def printverbose(f,flag,message):
 """Define verbose to print in file 'f', if 'flag'=True, message given as 'message'"""
 if flag: f.write("%s"%(message))

def main():
 options = parse_user_arguments()
 try:
  boxplotlist(options)
 except IOError as e:
  print("I/O error (%s): %s" %(e.errno, e.strerror))

def boxplotlist(options):

 if fileExist(options.input):
  dataname=[]
  dataplot=[]
  fp=open(options.input,"r")
  for line in fp:
   (name,datafile)=line.strip().split()
   if fileExist(datafile):
     dataname.append(name)
     datavalue=[]
     fi=open(datafile,"r")
     printverbose(sys.stdout,options.show,"Open %s\n"%(datafile))
     for datainfo in fi:
       (info,value)=datainfo.strip().split()
       datavalue.append(float(value))
     dataplot.append(np.array(datavalue))
     fi.close()
  fp.close()
 else:
  raise IOError("Input File not found")

 fig, ax1 = plt.subplots(figsize=(10,6))
 title='Boxplot of '+ options.input
 fig.canvas.set_window_title(title)
 plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

 printverbose(sys.stdout,options.show,"Data %s %s\n"%(repr(dataname),repr(dataplot)))
 bp = plt.boxplot(dataplot, notch=0, sym='+', vert=1, whis=1.5)
 plt.setp(bp['boxes'], color='black')
 plt.setp(bp['whiskers'], color='black')
 plt.setp(bp['fliers'], color='red', marker='+')

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
 if options.grid:
   ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
              alpha=0.5)

# Hide these grid behind plot objects
 ax1.set_axisbelow(True)

# Add Title of the plot
 ax1.set_title(options.title)

 ax1.set_xlabel('Samples')
 ax1.set_ylabel('Value')

# Fill the boxes with desired colors and hashes
 boxColors = ['white','darkkhaki','royalblue','salmon','lime','red','violet','cyan','blue',
             'yellow','green','navy','grey','orange','firebrick']
 boxHatch  = [' ','/','\\','|','-','+','x','o','.','*']

 if options.color==0: options.color=1
 if options.hashed==0: options.hashed=1
 numBoxes = len(dataplot)
 medians = range(numBoxes)
 k=0
 kk=0
# Statistics: compare all vs all data
 ks={}
 mw={} 
 for i in range(numBoxes):
  for j in range(numBoxes):
   try:
    ks.setdefault((i,j),stats.ks_2samp(dataplot[i],dataplot[j]))
   except:
    ks.setdefault((i,j),(0.0,1.0))
    e = sys.exc_info()[0]
    sys.stdout.write("Fail Kolmogorov Smirnov Test %s\n"%repr(e))
   try:
    mw.setdefault((i,j),stats.mannwhitneyu(dataplot[i],dataplot[j]))
   except:
    mw.setdefault((i,j),(0.0,0.5))
    e = sys.exc_info()[0]
    sys.stdout.write("Fail Mann-Whitney Test %s\n"%e)

#Boxplot
 for i in range(numBoxes):
   box = bp['boxes'][i]
   boxX = []
   boxY = []
   for j in range(5):
      boxX.append(box.get_xdata()[j])
      boxY.append(box.get_ydata()[j])
   boxCoords = zip(boxX,boxY)
  # Group by color and hashed lines
   n= k * options.color + i % options.color
   if n!=i:
    k=k+1
    if k>14: k=0
    kk=0
   j= i - k * options.color
   m= kk * options.hashed + j %  options.hashed
   if m!=j:
    kk=kk+1
    if kk>9: kk=0
   boxPolygon = Polygon(boxCoords, facecolor=boxColors[k], hatch= boxHatch[kk])
   ax1.add_patch(boxPolygon)
  # Now draw the median lines back over what we just filled in
   med = bp['medians'][i]
   medianX = []
   medianY = []
   for j in range(2):
       medianX.append(med.get_xdata()[j])
       medianY.append(med.get_ydata()[j])
       plt.plot(medianX, medianY, 'black',linewidth=2.5)
       medians[i] = medianY[0]
  # Finally, overplot the sample averages, with horizontal alignment
  # in the center of each box
   plt.plot([np.average(med.get_xdata())], [np.average(dataplot[i])],
           color='white', marker='*', markeredgecolor='black')

# Set the axes ranges and axes labelsa
 ax1.set_xlim(0.5, numBoxes+0.5)
 if options.top ==0:
  top = max([ max([x for x in dataplot[i]]) for i in range(numBoxes) ])
 else:
  top = options.top
 if options.bottom==0:
  bottom = min([ min([x for x in dataplot[i]]) for i in range(numBoxes) ])
 else:
  bottom=options.bottom
 ax1.set_ylim(bottom, top)
 xtickNames = plt.setp(ax1, xticklabels=dataname)
 plt.setp(xtickNames, rotation=45, fontsize=12)

# Finally, add a basic legend
 if options.legend:
  posy=0.12
  posx=0.7
  k=0
  for i in range(numBoxes):
   # Group by color and hashed lines
    n= k * options.color + i % options.color
    if n!=i:
     k=k+1
    plt.figtext(posx, posy,  dataname[i] ,
            backgroundcolor=boxColors[k], color='black', weight='roman',
            size='x-small')
    posy=posy-0.03
    if posy<0.01:
     posy=0.12
     posx=posx+0.12

  plt.figtext(0.14, 0.1, '*', color='white', backgroundcolor='silver',
           weight='roman', size='medium')
  plt.figtext(0.15, 0.095, ' Average Value', color='black', weight='roman',
           size='x-small')

 output=options.out+"."+options.format
 printverbose(sys.stdout,options.show,"output is %s\n"%output)
 plt.savefig(output)

 sys.stdout.write("Kolmogorov-Smirnov Statistic\n") 
 sys.stdout.write("+--------------------------+\n") 
 for i in range(numBoxes):
  for j in range(numBoxes):
    if j>=i :
     sys.stdout.write("%s versus %s : Statistic %f p-value %f Average(%s)= %f Average(%s)= %f\n"%(dataname[i],dataname[j],ks[(i,j)][0],ks[(i,j)][1],dataname[i],np.average(dataplot[i]),dataname[j],np.average(dataplot[j])))

 sys.stdout.write("\nMann-Whitney U Statistic\n") 
 sys.stdout.write("+----------------------+\n") 
 for i in range(numBoxes):
  for j in range(numBoxes):
    if j>=i :
     sys.stdout.write("%s versus %s : Statistic %f p-value %f Average(%s)= %f Average(%s)= %f\n"%(dataname[i],dataname[j],mw[(i,j)][0],2*mw[(i,j)][1],dataname[i],np.average(dataplot[i]),dataname[j],np.average(dataplot[j])))

     
 if options.canvas:
  mpl.use('Agg')
  plt.show()


if __name__=="__main__":
 main()



