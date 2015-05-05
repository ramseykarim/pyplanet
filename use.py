#! /usr/bin/env python
import os, os.path, argparse

constituents = os.listdir('constituents/')
#print constituents

o = argparse.ArgumentParser(prefix_chars='-+')
for c in constituents:
    cpth = os.path.join('constituents',c)
    if os.path.isdir(cpth):
        acon = '--'+c
        o.add_argument(acon,help='toggle '+c,action='store_true')
args = o.parse_args()
argdict = vars(args)
print argdict

for c in constituents:
    cpth = os.path.join('constituents',c)
    if os.path.isdir(cpth):
        fil = os.listdir(cpth)
        if 'use.txt' in fil:
            fp = open(os.path.join(cpth,'use.txt'),'r')
            s = fp.readline().strip()
            fp.close()
            if argdict[c]:
                print '--- %-8s Toggle to not use %s' % (c,s)
                os.rename(os.path.join(cpth,'use.txt'),os.path.join(cpth,'nouse.txt'))
            else:
                print '+++ %-8s Using %s' % (c,s)
        elif 'nouse.txt' in fil:
            fp = open(os.path.join(cpth,'nouse.txt'),'r')
            s = fp.readline().strip()
            fp.close()
            if argdict[c]:
                print '+++ %-8s Toggle to use %s' % (c,s)
                os.rename(os.path.join(cpth,'nouse.txt'),os.path.join(cpth,'use.txt'))
            else:
                print '--- %-8s Not using %s' % (c,s)
        else:
            print '    %-8s Use/nouse.txt not found' % (c)

