#! /usr/bin/env python
import os

print """Process is:
    make filenames.txt with input files
    'mkbatchplanet'   with filenames in filenames.txt
                      writes 'batchplanet'
    set any parameters (use.py, freqs in runp.py, ...)
    'batchplanet'     does runp.py with filenames from filenames.txt"""
print
print

try:
    filenames = open('filenames.txt','r')
except IOError:
    print 'filenames.txt not present'
else:
    batchfile = open('batchplanet','w')
    for line in filenames:
        s = 'runp.py '+line
        batchfile.write(s)
    os.chmod('batchplanet',0777)
