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
    
tractopgramIn=nib.streamlines.load(config['tractogram'])
streamlines=tractopgramIn.streamlines

def quickbundlesClusters(streamlines, **kwargs):
    """
    (quickly?, via qbx_and_merge?) perform a quick-bundling of an input
    collection of streamlines, and then extract centroids of the resultant
    clusters as streamlines.
    Parameters
    ----------
    streamlines : nibabel.streamlines.array_sequence.ArraySequence
        The input streamlines 
    **kwargs : keyword arguments for the qbx_and_merge
        Currently only supports [thresholds] and [nb_pts]
        See dipy.segment.bundles.qbx_and_merge for more details
    Returns
    -------
    centroidsAsStreamlines : nibabel.streamlines.array_sequence.ArraySequence
        A streamline object containing the centroids of the clusters resulting
        quickBundle-ification of the input streamlines
    clusters : dipy.segment.clustering.ClusterMapCentroid
        The clusters resulting from the quickBundle-ification of the input 
        streamlines
    """
    from dipy.segment.bundles import qbx_and_merge
    from dipy.tracking.streamline import Streamlines
    
    #fill in parameters if they are there.
    if not 'thresholds' in kwargs.keys():
        thresholds = [30,20,10,5]
    if not 'nb_points' in kwargs.keys():
        nb_pts=50
    #perform the quick, iterave bundling
    clusters=qbx_and_merge(streamlines,thresholds , nb_pts, select_randomly=None, rng=None, verbose=False)
    
    return clusters
    
def cullViaClusters(clusters,streamlines,streamThresh):
    import itertools
    #get the cluster lengths
    clusterLengths=[len(iCluster) for iCluster in clusters]
    
    #find which have meet the thresh criterion
    clustersSurviveThresh=np.greater(clusterLengths,streamThresh)
    
    #get a list of the clusters
    survivingClusters=list(itertools.compress(clusters,clustersSurviveThresh))
    #get the indexes of the streamlines from each
    survivingClusterLists=[iCluster.indices for iCluster in survivingClusters]
    #cat them all together
    survivingStreamsIndicies=list(itertools.chain(*survivingClusterLists))
    
    #find the obverse of the surviving stream set
    culledStreamIndicies=list(set(list(range(0,len(streamlines))))-set(survivingStreamsIndicies))

    return survivingStreamsIndicies, culledStreamIndicies

clustersOut=quickbundlesClusters(streamlines)

survivingStreamsIndicies, culledStreamIndicies= cullViaClusters(clustersOut.clusters,streamlines,streamThresh)

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
    surviveBool=np.ones(len(survivingStreamsIndicies),dtype=bool)
    survivorClass=wmaPyTools.streamlineTools.updateClassification(surviveBool,'survivingStreams',existingClassification=None)
    
    cullBool=np.ones(len(culledStreamIndicies),dtype=bool)
    culledClass=wmaPyTools.streamlineTools.updateClassification(surviveBool,'culledStreams',existingClassification=None)
    
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
    
