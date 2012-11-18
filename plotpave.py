#!/usr/bin/python

import pylab as pl

pave0=open("pavelog0.dat","r").readlines()
pave1=open("pavelog1.dat","r").readlines()

minE,delE,nBin=map(float,pave0[0].split())
nBin=int(nBin)
pave0.pop(0)
pave1.pop(0)
pave0=map(float,pave0)
pave1=map(float,pave1)

ebins=[i*delE+minE for i in range(nBin)]


pl.plot(ebins,pave0,label="Warm-up")
pl.plot(ebins,pave1,label="Total")
pl.legend(loc=0)
pl.show()
