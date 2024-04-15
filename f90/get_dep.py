#!/usr/bin/python

import re, sys, os

files = []
# r=root, d=directories, f = files
for r,d,f in os.walk("./"):
    for file in f:
        if file.endswith('.f90'):
            files.append([os.path.join(r,file),file[:-4]])

fout = open('./makedep.inc','w+')
for f in files:
    use_line_re = re.compile("^\s*use\s+(\S.+)\s*$")
    ff = open(f[0])
    mod = []
    for line in ff:
        match = use_line_re.search(line)
        if match:
            string = match.group(1).split('use ') #split when it finds 'use ' --> multiple use in single line
            for i in string:
                istrip=i.replace(';','').strip() #strip --> remove trail/end spaces
                istrip = istrip.split(',')
                istrip = istrip[0]
                if istrip[0:3] == 'mod': #string[0:3] --> start at 0, 3-0 characters
                    mod.append(istrip+'.f90')
    ff.close()
    nodup = [i for n, i in enumerate(mod) if i not in mod[:n]] # keep order
#    nodup=list(set(mod)) #do not keep order
    fout.write(f[1]+'.o: '+' '.join(nodup)+'\n')

fout.close()
            
