#!/usr/bin/python

from scipy import stats as stats
import numpy as np
import re,sys

if len(sys.argv) != 7:
    print("Error with arguments\n"+sys.argv[0]+"\n")
    sys.exit(0)

try:
    DEPTH = open(sys.argv[1])
    LOCATION = open(sys.argv[2])
    OUT = open(sys.argv[3],"a")
    OUTNORMAL = open(sys.argv[4],"a")
    depth = sys.argv[5]
    flankingLength = sys.argv[6]
except IOError:
    print("IOError")

for line1 in LOCATION:
    list_hgt=[]
    list_flanking150bp=[]
    arr1 = re.split("\s+",line1.rstrip())
    ori = arr1[2]
    arr3 = re.split(r"\|",arr1[3].rstrip())
    arr4 = re.split(r"-",arr3[1].rstrip())
    start = arr4[0]
    end = arr4[1]
    chro = arr3[0]
    sign = arr3[2]
    length = int(end) - int(start) + 1
    num = 0
    try:
        DEPTH = open(sys.argv[1])
    except IOError:
        print("IOError")
    for line in DEPTH:
        arr = re.split(r"\||\s+",line.rstrip())
        now = arr[0]+"|"+arr[1]+"|"+arr[2]
        if ori==now:
            arr2 = re.split(r"-",arr[1].rstrip())
            start1 = arr2[0]
            end1 = arr2[1]
            location = int(start1) + int(arr[3]) - 1
            if int(start)<=int(location) and int(location)<int(end):
                list_hgt.append(int(arr[4]))
                if int(arr[4]) > 2:
                    num += 1
            elif int(start) - int(flankingLength) <= int(location) and int(location) < int(end) + int(flankingLength):
                list_flanking150bp.append(int(arr[4]))
    DEPTH.close()

    hgt=np.array(list_hgt)
    flanking150bp=np.array(list_flanking150bp)
    pvalue = str(stats.ttest_ind(list_hgt, list_flanking150bp, equal_var = False))
    p = pvalue.split('pvalue=')[1].split(',')[0]
    print(pvalue) 
   
    if((int(length)-int(num)>2) or (float(np.mean(hgt))<float(depth)) or (float(np.mean(hgt)) < float(np.mean(flanking150bp)) and float(p) < 0.05)):
        OUT.write(chro + "|" + start + "-" + end + "|" + sign + "\t0.00\n")
    else:
        OUT.write(chro + "|" + start + "-" + end + "|" + sign + "\t" + str(float(int(num)/int(length))) + "\n")
  
    if((int(length)-int(num)>2) or (float(np.mean(hgt))<float(depth))):
        OUTNORMAL.write(chro + "|" + start + "-" + end + "|" + sign + "\t0.00\n")
    else:
        OUTNORMAL.write(chro + "|" + start + "-" + end + "|" + sign + "\t" + str(float(int(num)/int(length))) + "\n")

    print(DEPTH)
    print(LOCATION)
    print(np.mean(hgt))
    print(np.mean(flanking150bp))
    print(p)
    print(np.var(hgt))
    print(np.var(flanking150bp))

