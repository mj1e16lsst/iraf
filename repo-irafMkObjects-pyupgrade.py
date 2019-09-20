
# coding: utf-8

# In[1]:


import keplerSettings as keplerSettings

import subprocess
import os
import csv

from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.io import fits
from astropy.io import ascii

import matplotlib
import matplotlib.pyplot as plt

import numpy as np

from pyraf import gwm

import subprocess


# In[2]:


from pyraf import iraf
from iraf import obsutil
from iraf import psfmeasure
from iraf import phot
from iraf import daofind
from iraf import daophot
from iraf import psf
from iraf import artdata
from iraf import mkobjects
from iraf import starlist


# In[18]:


imagename = keplerSettings.imageName
diffImage = keplerSettings.differenceImageName
#SImage= '/home/mj1e16/Simages/testsimage.fits'
ext =  79 # keplerSettings.ccdExtensions[0]


# In[19]:


outDir = keplerSettings.irafOutputDir # for iraf output files ONLY
imageDir = keplerSettings.irafImageDir # for images generated during process ONLY
simageDir = keplerSettings.simImageDir # directory for simulated images
irafDir = keplerSettings.irafDir #Directory containing iraf installation
starlistDir = keplerSettings.starlistDir # directory containing starlist files
workflowDir = keplerSettings.workflowDir


# In[20]:


def mkobjectspluspsf(imagename,diffimage,diffPlusSim,ext,section,minMag,maxMag,Xmax=1000,Xmin=0,Ymax=1000,Ymin=1000,thresh=50000,nsegs=keplerSettings.nsegs,goodpsf='',spatialDist='uniform',zeropoint=keplerSettings.imageZeroPoint):
    
    daoout = outDir+'{}_{}daofindone'.format(ext,section)
    photout = outDir+'{}_{}photone'.format(ext,section)
    psfout = outDir+'{}_{}psfone'.format(ext,section)
    groupout = outDir+'{}_{}groupfileone'.format(ext,section)
    opstout = outDir+'{}_{}optsone'.format(ext,section)
    starfieldsect = outDir+'starfield_{}.dat'.format(section)
    imseg = imageDir+'imageseg{}.fits'.format(section)
    
    hdu = fits.open(imagename)
    print(imagename)
    image = hdu[ext].section[Xmin:Xmax,Ymin:Ymax]
    hdu = fits.PrimaryHDU(image)
    hdu.writeto(imseg,clobber=True)
    
    with open(starfieldsect,'w') as f:
        f.write("")
    
    daofind(imseg,output=daoout,sigma=1.0,threshold=thresh)
    phot(imseg,skyfile=daoout,coords=daoout,output=photout,interactive='no')
    psf(imseg,photfile=daoout,pstfile=photout,psfimage=psfout,opstfile=opstout,groupfile=groupout,interactive='no')
    starlist(starfieldsect,nstars=1000/nsegs,xmin=Xmin,xmax=Xmax,ymin=Ymin,ymax=Ymax,spatial=spatialDist,minmag=minMag,maxmag=maxMag)
    subprocess.call(['cp',starfieldsect,keplerSettings.starlistDir+'/starfield_{}_{}.dat'.format(section,minMag)])
    try:
        mkobjects(diffimage,output=diffPlusSim,objects=starfieldsect,gain=110,rdnoise=127,star=psfout+'.fits')
        
    except:
        try:
            #subprocess.call(['rm','/home/mj1e16/Simages/testsimage*'])
            mkobjects(diffimage,output=diffPlusSim,objects=starfieldsect,gain=110,rdnoise=127,star=goodpsf)
            print('Waring: Section {} with bad PSF, {} used inplace'.format(section,goodpsf))
        except:
            print('ERROR: PSF for Section {} should be made manually'.format(section))
    
    return psfout+'.fits'


# In[27]:


def findPSFandMkObjects(imageName,diffImage,extension,minMag,maxMag,imageshape=[keplerSettings.astroImageXlength, keplerSettings.astroImageYlength],pixelsize=3.98):
    
    print(imageName)
    
    os.chdir(irafDir)
    
    xlength = imageshape[1]
    ylength = imageshape[0]
    nsegs = keplerSettings.nsegs
    
    xsegment = xlength/nsegs
    ysegment = ylength/nsegs
    PSF = []
    goodpsf=''
    for xsegs in range(nsegs):
        lowx = int(xsegs*xsegment)+keplerSettings.border
        highx = int(lowx+xsegment)+keplerSettings.border-10
        psf = []
        for ysegs in range(0,nsegs):
            lowy = int(ysegs*ysegment)+keplerSettings.border
            highy = int(lowy+ysegment)+keplerSettings.border-10
            
            
            # stars are being lost from images in xseg=0
            SImage= imageDir+'testsimage_{}_{}_{}.fits'.format(xsegs,ysegs,minMag)
            
            goodpsf = mkobjectspluspsf(imageName,diffImage,SImage,extension,'{}_{}'.format(xsegs,ysegs),minMag,maxMag,Xmax=highx,Xmin=lowx,Ymax=highy,Ymin=lowy,nsegs=nsegs**2,goodpsf=goodpsf)
            
            b = diffImage.find('{}_alt'.format(minMag))
            if b == -1:
                diffImage = diffImage[:-5]+'{}_alt.fits'.format(minMag)
            try:
                hdu = fits.open(SImage) #,ignore_missing_end=True)
                imdata = hdu[0].data
                hdu = fits.PrimaryHDU(imdata)
                hdu.writeto(diffImage,clobber=True)
                print('find me')
                print(SImage)
                print(diffImage)
            except:
                print('bad at section {} {}'.format(xsegs,ysegs))

            

            
    return 'complete!'


# In[22]:


def cleanDirectories(dirlist):
    # clean iraf output and input images
    # make sure to have no needed files in these directories
    for directory in dirlist:
        dirfiles = os.listdir(directory)
        for eachFile in dirfiles:
            #print('rm',directory+eachFile)
            subprocess.call(['rm',directory+eachFile])


# In[23]:


def stitchStarList(nsegs,mag,starDir=starlistDir):
    with open(starDir+'starfield_0_0_{}.dat'.format(mag)) as f:
        lines = f.readlines()
    bigstring = ['# X_POS Y_POS MAG\n']
    bigstring.extend(lines[0:18])
    for x in range(nsegs):
        for y in range(nsegs):
            with open(starDir+'starfield_{}_{}_{}.dat'.format(x,y,mag)) as f:
                lines = f.readlines()
            bigstring.extend(lines[18:])


    with open(workflowDir+'starlistFull_{}.dat'.format(mag),'w') as f:
        f.write(''.join(bigstring))


# In[24]:


cleanDirectories([outDir,imageDir]) # cleans iraf output files, will remove all files from these directories


# In[28]:


run = 'not complete'
ext =  63
medianNum = 1
diffImage = simageDir+'diff_{}_{}.fits'.format(ext,medianNum)
for mag in keplerSettings.magRange:
    cleanDirectories([outDir,imageDir]) # cleans iraf output files, will remove all files from these directories
    for x in range(5):
        try:
            run = findPSFandMkObjects(imagename,diffImage,ext,mag,(mag+1))
            if run == 'complete!':
                break
        except:
            pass
    dirlist = os.listdir(starlistDir)
    dirlist = [x for x in dirlist if '.dat' in x]
    for x in dirlist:
        subprocess.call(['mv',starlistDir+x,starlistDir+'{}/{}/'.format(ext,medianNum)])


# In[62]:


mag = -7
for x in range(5):
        try:
            run = findPSFandMkObjects(imagename,diffImage,SImage,ext,mag,(mag+1))
            if run == 'complete!':
                break
        except:
            pass


# In[40]:


for mag in range(-7,0):
    stitchStarList(4,mag,starDir=starlistDir+'{}/{}/'.format(ext,1))


# In[17]:


# def makeDS9RegFile(sexTabList,fileNameBase,tabType):
#     for tables in range(len(sexTabList)):
#         if tabType == 'dao':
#             xcoords = sexTabList[tables]['xcentroid']
#             ycoords = sexTabList[tables]['ycentroid']
#         elif tabType == 'sex':
#             xcoords = sexTabList[tables]['X_IMAGE']
#             ycoords = sexTabList[tables]['Y_IMAGE']
#         elif tabType == 'iraf':
#             xcoords = sexTabList[tables]['X_POS']
#             ycoords = sexTabList[tables]['Y_POS']        
#         else:
#             print('tab type error')
#             break
#         bigString = 'global color=lightgreen\nimage\n'
#         for x in range(len(xcoords)):
#             bigString += 'circle({},{},5)\n'.format(xcoords[x],ycoords[x])
#         fileName = fileNameBase + str(tables)+'.reg'
#         with open(fileName,'w') as f:
#             f.write(bigString)
#         print(fileName)


# In[15]:


#data = ascii.read('/home/mj1e16/keplerPhotometry/starlistFull.dat')


# In[16]:


#makeDS9RegFile([data],'/home/mj1e16/Simages/teststarlist','iraf')


# In[2]:


# for x in range(-7,0):
#     print(x)


# In[ ]:




