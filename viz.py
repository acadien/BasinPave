#!/usr/bin/python

import pylab as pl
import sys
#mine
from poscarsPlot import plotsimulation
from struct_tools import neighbors

data=open(sys.argv[1],"r").readlines()
if data[0][0]=="#":
    data.pop(0)
atoms=list()
msdNeighb=list()
msdIdeal=list()
energies=list()
volumes=list()
natoms=0


natoms=int(data[1])
pos=14+natoms
for i,line in enumerate(data):
    j=i%pos
    if(j==3):
        energies.append(float(line))
    if(j==5):
        volumes.append(float(line))
    if(j==9):
        atoms.append(list())
    if(j>=9 and j<9+natoms):
        atoms[-1].append(map(float,line.split()))
    if(j==natoms+10):
        msdNeighb.append(float(line))
    if(j==natoms+12):
        msdIdeal.append(float(line))
optimE={38:-173.252,76:-402.384,104:-582.038}
nState=len(msdNeighb)
basis=[[2,0,0],[0,2,0],[0,0,2]]
#Figure 1: MSD & Energy vs step #
pl.subplot(211)
pl.title("%d Atoms"%natoms)
pl.ylabel("Energy")
pl.plot([0,nState],[optimE[natoms]/natoms,optimE[natoms]/natoms],color="black",ls="--",label="Optimal Energy/atom")
pl.plot(range(nState),[i/natoms for i in energies],label="Energy/atom")
pl.subplot(212)
pl.plot(range(nState),msdIdeal,label="MSD-ideal")
pl.legend(loc=0)
pl.xlabel("Step")


#Figure 2: Visualize Lowest Energy Configuration
fig2=pl.figure()
#a2=fig2.add_subplot(111,project='3d')
mc=energies.index(min(energies))
iatoms=atoms[mc]
ienergy=energies[mc]
print ienergy
aa=plotsimulation(basis,iatoms,[natoms],fig2)
bnd=[[0.0,11.0],[0.0,11.0],[0.0,11.0]]
halfNeigh=neighbors(iatoms,bnd,1.5,"half")
for a in range(natoms):
    for b in halfNeigh[a]:
        x,y,z=zip(*[iatoms[a],iatoms[b]])
        aa.plot(x,y,z,c="black",lw=1.5)
pl.show()
