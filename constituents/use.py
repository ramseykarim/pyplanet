#! /usr/bin/env python
import os, os.path, argparse

constituents = os.listdir('.')

o = argparse.ArgumentParser(prefix_chars='-+')
for c in constituents:
    if os.path.isdir(c):
        acon = '--'+c
        o.add_argument(acon,help='toggle '+c,action='store_true')
args = o.parse_args()
argdict = vars(args)

for c in constituents:
    if os.path.isdir(c):
        fil = os.listdir(c)
        if 'use.txt' in fil:
            fp = open(os.path.join(c,'use.txt'),'r')
            s = fp.readline().strip()
            fp.close()
            if argdict[c]:
                print '--- %-8s Toggle to not use %s' % (c,s)
                os.rename(os.path.join(c,'use.txt'),os.path.join(c,'nouse.txt'))
            else:
                print '+++ %-8s Using %s' % (c,s)
        elif 'nouse.txt' in fil:
            fp = open(os.path.join(c,'nouse.txt'),'r')
            s = fp.readline().strip()
            fp.close()
            if argdict[c]:
                print '+++ %-8s Toggle to use %s' % (c,s)
                os.rename(os.path.join(c,'nouse.txt'),os.path.join(c,'use.txt'))
            else:
                print '--- %-8s Not using %s' % (c,s)
        else:
            print '    %-8s Use/nouse.txt not found' % (c)

