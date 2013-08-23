#! /usr/bin/env python
import os, os.path

constituents = os.listdir('.')

for c in constituents:
    if os.path.isdir(c):
        fil = os.listdir(c)
        if 'use.txt' in fil:
            fp = open(os.path.join(c,'use.txt'),'r')
            s = fp.readline().strip()
            print 'Using '+c+' ('+s+')'
        else:
            print 'Not using '+c

