#!/usr/bin/env python3

import json
import numpy as np
import nibabel as nib
from dipy.segment.clustering import QuickBundles
from dipy.io.streamline import load_tractogram
from dipy.io.streamline import save_tractogram
import seaborn as sns
import matplotlib.pyplot as plt

def jsonParseNPass():
    # load config.json structure, modified because of windows filepath issues
    with open('config.json','r') as config_f:
        jsonString = config_f.read()
        #windows escape character nonsense
        config = json.loads(jsonString.replace('\\','/'))
        
        # load data
        reference = nib.load(config['dwi'])
        tractogram = load_tractogram(config['tractogram'],reference)
        streamThresh=config['minStreamsInCluster']
        return reference, tractogram, streamThresh

def cullViaClusters(clusters,tractogram,streamThresh):
    clustersSurviveThresh=np.greater(np.asarray(list(map(len, clusters))),streamThresh)
    survivingStreams=[]
    for iclusters in clusters[clustersSurviveThresh]:
        survivingStreams=survivingStreams + iclusters.indices
    culledStreamIndicies=set(list(range(1,len(tractogram.streamlines))))-set(survivingStreams)
    #cull those streamlines
    #don't know what to do about those warnings
    tractogram.streamlines=tractogram.streamlines[survivingStreams]
    culledTractogram=tractogram.copy()
    culledTractogram.streamlines=culledTractogram.streamlines[culledStreamIndicies]
    return tractogram, culledTractogram

#obtain rquisite files from config.json
[reference, tractogram, streamThresh]=jsonParseNPass()
#threshold heuristic designed to target a threshold of 15 at 500k streamlines
#not actually accurate as relationship isn't linear
inputStreamNumber=len(tractogram.streamlines)
qbThresh=np.divide(inputStreamNumber,33333)
#perform quickBundles
qb = QuickBundles(threshold=qbThresh)
clusters = qb.cluster(tractogram.streamlines)

#perform cull
[outTractogram,culledTractogram]=cullViaClusters(clusters,tractogram,streamThresh)

#save tractograms
save_tractogram(outTractogram, 'track.tck')
save_tractogram(culledTractogram, 'culled.tck')

numberCulled=inputStreamNumber-len(outTractogram.streamlines)
print(str(numberCulled) + ' streamlines culled')
#create output diagnostic plots
sns.histplot(data=np.asarray(list(map(len, clusters))),log_scale=True, binwidth=.1)
#I dont know why I have to enter it in for both values
plt.plot([streamThresh,streamThresh], [  1,100], linewidth=4,color='red')
plt.savefig('clusterCullHist.png')
#not working, fix later
#plt.xlabel='cluster size'
#plt.title='cluster size counts'

# end