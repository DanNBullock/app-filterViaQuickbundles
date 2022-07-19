#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 17:56:15 2022

@author: dan
"""

import os
import sys

sys.path.append('wma_pyTools')
# startDir=os.getcwd()

# import wmaPyTools.roiTools
# import wmaPyTools.analysisTools
# import wmaPyTools.segmentationTools
import wmaPyTools.streamlineTools
# import wmaPyTools.visTools

#os.chdir(startDir)

import os
import json
import numpy as np
import nibabel as nib
import shutil
import itertools


# load inputs from config.json
with open('config.json') as config_json:
	config = json.load(config_json)
    
streamThresh=config['streamThresh']

#try to get the qbThreshes key
try:
    #but it comes in as a string so
    qbThreshes=[int(iqbThresh) for iqbThresh in config['qbThreshes'].split(',')]
except:
    qbThreshes=[30,20,10,5]
    
tractogramIn=nib.streamlines.load(config['tractogram'])
streamlines=tractogramIn.streamlines

#quickbundlesClusters(streamlines, thresholds = [30,20,10,5], nb_pts=100,verbose=False)
clustersOut=wmaPyTools.streamlineTools.quickbundlesClusters(streamlines,thresholds = qbThreshes)

#find surviving and culled indexes
survivingStreamsIndicies, culledStreamIndicies= wmaPyTools.streamlineTools.cullViaClusters(clustersOut.clusters,streamlines,streamThresh)
#create bool vecs for both outputs
survivingStreamsBoolVec=np.zeros(len(streamlines),dtype=bool)
survivingStreamsBoolVec[survivingStreamsIndicies]=True
culledStreamsBoolVec=np.zeros(len(streamlines),dtype=bool)
culledStreamsBoolVec[culledStreamIndicies]=True

#start figuring out how you're going to output this...

from scipy.io import savemat

try:
    classificationDict=wmaPyTools.streamlineTools.matWMC2dict(config['classification'])
    
    survivorClass={}
    culledClass={}
    survivorClass['names']=classificationDict['names']
    culledClass['names']=classificationDict['names']
    survivorClass['index']=list(itertools.compress(classificationDict['index'],survivingStreamsIndicies))
    culledClass['index']=list(itertools.compress(classificationDict['index'],culledStreamIndicies))
    
    #create some output 
    #wmc
    print('saving surviving streams wmc')
    if not os.path.exists(os.path.join('wmc_survived')):
        os.makedirs(os.path.join('wmc_survived'))
    savemat(os.path.join('wmc_survived','classification.mat'),{ "classification": {"names": np.array(survivorClass['names'], dtype=np.object), "index": survivorClass['index'] }})
    print('saving culled streams wmc')
    if not os.path.exists(os.path.join('wmc_culled')):
        os.makedirs(os.path.join('wmc_culled'))
    savemat(os.path.join('wmc_culled','classification.mat'),{ "classification": {"names": np.array(culledClass['names'], dtype=np.object), "index": culledClass['index'] }})
          
    #tck
    print('saving surviving streams tck')
    if not os.path.exists(os.path.join('tck_survived')):
        os.makedirs(os.path.join('tck_survived'))
    streamSubset=streamlines[survivingStreamsIndicies]
    wmaPyTools.streamlineTools.stubbornSaveTractogram(streamSubset, os.path.join('tck_survived','track.tck'))
    print('saving culled streams tck')
    if not os.path.exists(os.path.join('tck_culled')):
        os.makedirs(os.path.join('tck_culled'))
    streamSubset=streamlines[culledStreamIndicies]
    wmaPyTools.streamlineTools.stubbornSaveTractogram(streamSubset,  os.path.join('tck_culled','track.tck'))
     
    #if there's no input file classification
except:
    #wmc doesn't really matter, just create uniform vec outputs for both
    #have to do this because of brainilfe output conventions
    #surviveBool=np.ones(len(survivingStreamsIndicies),dtype=bool)
    survivorClass=wmaPyTools.streamlineTools.updateClassification(survivingStreamsBoolVec,'survivingStreams',existingClassification=None)
    
    #MAcullBool=np.ones(len(culledStreamIndicies),dtype=bool)
    culledClass=wmaPyTools.streamlineTools.updateClassification(culledStreamsBoolVec,'culledStreams',existingClassification=None)
    
    #create some output 
    #wmc
    print('saving surviving streams wmc')
    if not os.path.exists(os.path.join('wmc_survived')):
        os.makedirs(os.path.join('wmc_survived'))
    savemat(os.path.join('wmc_survived','classification.mat'),{ "classification": {"names": np.array(survivorClass['names'], dtype=np.object), "index": survivorClass['index'] }})
    print('saving culled streams wmc')
    if not os.path.exists(os.path.join('wmc_culled')):
        os.makedirs(os.path.join('wmc_culled'))
    savemat(os.path.join('wmc_culled','classification.mat'),{ "classification": {"names": np.array(culledClass['names'], dtype=np.object), "index": culledClass['index'] }})
          
    #tck
    print('saving surviving streams tck')
    if not os.path.exists(os.path.join('tck_survived')):
        os.makedirs(os.path.join('tck_survived'))
    streamSubset=streamlines[survivingStreamsIndicies]
    wmaPyTools.streamlineTools.stubbornSaveTractogram(streamSubset, os.path.join('tck_survived','track.tck'))
    print('saving culled streams tck')
    if not os.path.exists(os.path.join('tck_culled')):
        os.makedirs(os.path.join('tck_culled'))
    streamSubset=streamlines[culledStreamIndicies]
    wmaPyTools.streamlineTools.stubbornSaveTractogram(streamSubset, os.path.join('tck_culled','track.tck'))
    
