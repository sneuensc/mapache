#!/usr/bin/env python

from __future__ import division

import argparse
import pysam
import sys,math
from collections import defaultdict

def average(values):
    '''Computes the arithmetic mean of a list of numbers.'''
    return sum(values, 0.0) / len(values)


parser = argparse.ArgumentParser(description='''Calculate average depth from bam file.''')

parser.add_argument('i', help='input bamfile')
parser.add_argument("--malelimit", action="store", type=float, dest="malelimit",help="Upper R_y limit for assignment as XY/male",default=0.075)                                                            
parser.add_argument("--femalelimit", action="store", type=float, dest="femalelimit",help="Lower R_y limit for assignment as XX/female",default=0.016)                                                      
parser.add_argument("--digits", action="store", type=int, dest="digits",help="Number of decimal digits in R_y output",default=4)
parser.add_argument("--noheader", action="store_true", dest="noheader",help="Do not print header describing the columns in the output",default=False)
parser.add_argument("--female", action="store", type=str, dest="femaleStr",help="Female chromosome name",default="X")
parser.add_argument("--male", action="store", type=str, dest="maleStr",help="Male chromosome name",default="Y")

args = parser.parse_args()

#args = parser.parse_args('libA_hg19.sort.q30.bam.rmdup.bam'.split())
#args = parser.parse_args('alignment/kleb_10_213361.flt.sort.rmdup.bam'.split())
#args = parser.parse_args('A.flt.sort.bam'.split())


# read in from bam
samfile = pysam.Samfile(args.i, 'rb')

# get length of references from header
refl = defaultdict(lambda:0) 
header = samfile.header['SQ']
lengths = {}
for ref in header:
   refl[ref['SN']] = ref['LN']
   lengths[ref['SN']] = []

# parse bam
total_count = 0
not_aligned = 0

d_counts = defaultdict(lambda:0)
d_sums = defaultdict(lambda:0)
s = 0
with samfile as bam_file:
   for alignedRead in bam_file:
      total_count += 1
      if alignedRead.flag == 4:
         not_aligned += 1
         continue
      refname = samfile.getrname( alignedRead.rname )
      d_counts[refname] += 1
      d_sums[refname] += alignedRead.rlen
      s += alignedRead.rlen


print ('Counted %i, with %i not aligned' % (total_count, not_aligned))

if total_count > 0:
    print ('Average length of mapped reads %g' % (s/total_count))
else:
	print ('There are no mapped reads %g' % (total_count)
)
# output
print ('Average depth of %s:' % args.i)
print ('Chrom\tAvgReadDepth\tSumReadLength\tChromLength\tAverageReadLength')
totalLength = 0
totalReadlength = 0
totalReadCount = 0
chrY_counts = None ; chrX_counts = None;

for key in d_counts.keys():
    if args.maleStr in key:
        chrY_counts = d_counts[key]
    else:
    	chrY_counts = 0
    if args.femaleStr in key:
        chrX_counts = d_counts[key]
    else:
    	chrX_counts = 0

    totalLength += refl[key]
    totalReadlength += d_sums[key]
    totalReadCount += d_counts[key]
    avg = d_sums[key] / refl[key]
    avgReadL = 1.0*d_sums[key]/d_counts[key]
    print('%s: %.3f\t%i\t%i\t%.1f' % (key, avg, d_sums[key], refl[key], avgReadL))


if total_count > 0:
    totalAvgReadL = 1.0*totalReadlength/totalReadCount
    print ('Total: %.3f\t%i\t%i\t%.1f' % (totalReadlength / totalLength, totalReadlength, totalLength, totalAvgReadL))
else:
    print ('Total: %.3f\t%i\t%i\t%.1f' % (0, 0, 0, 0))


if chrX_counts > 0:
    print ("Sex assignment (found chroms containing X or Y)")
    print ('nreads on Y: ',chrY_counts)
    print ('nreads on X: ',chrX_counts)
    
    n = chrY_counts+chrX_counts  
    Ry = 1.*chrY_counts/n

    SE = math.sqrt((Ry*(1.0-Ry))/n)      
    confinterval = 1.96*SE
 
    #use criteria to infer chromosomal sex
    gender='NA'
    if (Ry < args.femalelimit) and (Ry > args.malelimit):
        gender='Not Assigned'
    elif Ry==0.0:
        gender='consistent with XX'
    elif Ry+confinterval < args.femalelimit:
        gender='XX'
    elif Ry-confinterval > args.malelimit:
        gender='XY'
    elif Ry-confinterval > args.femalelimit and Ry+confinterval > args.malelimit:
        gender='consistent with XY but not XX'
    elif Ry-confinterval < args.femalelimit and Ry+confinterval < args.malelimit:
        gender='consistent with XX but not XY'
    else:
        gender='Not Assigned'

    if args.noheader == False:
        print ('Nseqs\tNchrY+NchrX\tNchrY\tR_y\tSE\t95% CI\tAssignment')
        print (total_count,'\t',n,'\t',chrY_counts,'\t',round(Ry,args.digits),'\t',round(SE,args.digits),'\t',str(round(Ry-confinterval,args.digits))+'-'+str(round(Ry+confinterval,args.digits)),'\t',gender)


