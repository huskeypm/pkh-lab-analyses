"""
For processing ODE outputs in support of Satin/Despa collaborations 
"""

from matplotlib.ticker import ScalarFormatter
import math
import cPickle as pickle
import runner 
import matplotlib.pylab as plt
import numpy as np
#import analyzeODE as ao 
import downSamplePickles as dsp
from runShannonTest import *
import json
import taufitting as tf
import scipy.fftpack as fftp
runner.init()
import os.path

import analyzeGotran as anG 

class empty:pass
mM_to_uM = 1e3
ms_to_s = 1e-3

plotlyAuth=None

# For plotly support 
def PlotViaPlotly(casesSubset,state):
  downsample=10
  import plotly
  if plotlyAuth==None:
    with open('/net/share/pmke226/PLOTLY', 'r') as f:
        first_line = f.readline()
    plotlyKey = first_line
    plotly.tools.set_credentials_file(
      username='huskeypm', api_key='3x5rz5d19r') #plotlyKey)
  import plotly.tools as tls
  import plotly.plotly as py

  fig = plt.figure()
  ax = plt.subplot(111)
  case1=casesSubset[0]
  case2=casesSubset[1]
  title ="%s_%s_%s"%(state,case1.label,case2.label)
  ax.set_title("%s: %s,%s"%(state,case1.label,case2.label))
  for i,case in enumerate(casesSubset):
    pkg = GetData(case.data,state)
    ax.plot(pkg.t[::downsample],pkg.valsIdx[::downsample], label=case.label)
  plotly_fig = tls.mpl_to_plotly( fig )

  # Adding custom attributes to legend
  plotly_fig['layout']['showlegend'] = True
#layout = go.Layout(
#    xaxis=dict(
#        range=[2, 5]
#    ),
#    yaxis=dict(
#        range=[2, 5]
#    )
#)
#fig = go.Figure(data=data, layout=layout)

  plot_url = py.iplot(plotly_fig, filename = title)
  print "%s",plot_url.resource

### 
### I/O 
###
def makePackage(p,p_idx,s,s_idx,j,j_idx,t):

  return {'p':p,'s':s,'t':t,'j':np.asarray(j),\
           'p_idx':p_idx,'s_idx':s_idx,'j_idx':j_idx}
  
def writePickle(name,p,p_idx,s,s_idx,j,j_idx,t):
  # store to pickle
  # using 'asarray' since my 'j' was getting stored as its transpose 
  #data1 = {'p':p,'s':s,'t':t,'j':np.asarray(j),\
  #         'p_idx':p_idx,'s_idx':s_idx,'j_idx':j_idx}
  
  data1 = makePackage(p,p_idx,s,s_idx,j,j_idx,t)

  #print "j again: ", len(j) 
  #print "j_idx: ",np.shape(j_idx)
  print "name: ", name
  #if ".pkl" not in name[-4:]: 
  if ".pkl" not in name and ".pickle" not in name:
    name += ".pkl"
  output = open(name, 'wb')
  pickle.dump(data1, output)
  output.close()
  print "SUCCESS! Wrote output to", name


import re
def readPickle(name = "PCa0.75kss0.25.pkl",verbose=True,readSubset=None,readConcat=False):          
  raise RuntimeError("Antiquated - use analyzeGotran") 

  if readConcat:
    print name 
    name = re.sub(r'\d+\.p', 'cat.p', name)  # pickle 
    print "Reading concatenated file %s instead" % name

  if verbose: 
    print "Reading " + name  
  pkl_file = open(name, 'rb')
  data1 = pickle.load(pkl_file)  
  pkl_file.close()

  if readSubset!=None:
    print "HARDCODED TO Cai, Ca_jct, Ca_SR for now"
    subsets = ["Cai","Ca_jct1","Ca_SR"]
    nTs = np.shape( data1["t"] )[0]
    ar = np.zeros([nTs, len(subsets)]) 
    si = data1["s"]
    s_idx = data1["s_idx"]
    for i,subset in enumerate(subsets):
      idx = s_idx.index(subset)
      ar[:,i] =  si[:, idx ] 
    # overwrite 
    data1["s"] = ar                  
    data1["s_idx"] = subsets             
  

  return data1  
# subset: None - load all attributes
#         [state1,state2, ...] 
def LoadPickles(caseDict,noOverwrite=False,
                verbose=True,readSubset=None,readConcat=False):
  raise RuntimeError("Antiquated - use analyzeGotran") 
  for key,case in caseDict.iteritems():
    if verbose:
      print "# ", key
      print "Loading"  , case.name

    if hasattr(case,'data') and noOverwrite==True:
      print "Skipping read, since already populated"
    else: 
      case.data = readPickle(case.name,verbose=verbose,readSubset=readSubset, 
                             readConcat=readConcat) 

### taken from fitter.py, used to process data made to put into panda format at the end.
def ProcessDataArray(dataSub,mode,timeRange=[0,1e3],key=None):

      # PRINT NO, NEED TO PASS IN TIME TOO 
      timeSeries = dataSub.t
      idxMin = (np.abs(timeSeries-timeRange[0])).argmin()  # looks for entry closest to timeRange[i]
      idxMax = (np.abs(timeSeries-timeRange[1])).argmin()
      valueTimeSeries = dataSub.valsIdx[idxMin:idxMax]
      #print "obj.timeRange[0]: ", obj.timeRange[0]
      #print "valueTimeSeries: ", valueTimeSeries

      tRange = timeSeries[idxMin:idxMax] - timeSeries[idxMin]
      waveMax = np.argmax(valueTimeSeries)
      tRangeSub = tRange[waveMax:]
      caiSub = valueTimeSeries[waveMax:]

      if key=="Cai": # for debug
        np.savetxt("test%d"%tag,valueTimeSeries)
      #print "dataSub.valsIdx: ", dataSub.valsIdx 
      if mode == "max":
          result = np.max(valueTimeSeries)
      elif mode == "min":
          result = np.min(valueTimeSeries)
      elif mode == "mean":
          result = np.mean(valueTimeSeries)
      elif mode == "amp":
          result = (np.max(valueTimeSeries) - np.min(valueTimeSeries))
      elif mode == "APD":
          #daMaxPlace = 0
          #daMaxHalfPlace = 0
          waveMin = np.argmin(caiSub)
          APDSub = caiSub[:waveMin]
          APDShift = APDSub + abs(APDSub[-1])
          APDMean = (APDShift[0] + APDShift[-1]) / 2
          APDSearch = APDShift - APDMean
          APDTimeidx = np.argmax(APDSearch <= 0)
          #daMin = abs(np.min(valueTimeSeries))
          #daMax = (np.max(valueTimeSeries) + daMin)
          #daMaxHalf = (daMax / 2)
          #for i, val in enumerate(valueTimeSeries):
              #print "val - daMaxHalf: ", (val-daMaxHalf)
          #    if daMax == (val + daMin):
                 #print "did we get here"
          #       daMaxPlace = i
          #    if (val + daMin) - daMaxHalf >= 0:
                 #print "how about here"
          #       daMaxHalfPlace = i 
          result = APDTimeidx * 0.1 * ms_to_s
          #result = (daMaxHalfPlace - daMaxPlace) * 0.1 * ms_to_s
          #result = valueTimeSeries[daMaxPlace] - valueTimeSeries[daMaxHalfPlace]
      elif mode == "tau":
          #tRange = timeSeries[idxMin:idxMax] - timeSeries[idxMin] 
          #waveMax = np.argmax(valueTimeSeries)
          #tRangeSub = tRange[waveMax:]
          #caiSub = valueTimeSeries[waveMax:] 

          fitted = tf.FitExp(tRangeSub,caiSub)
          result = fitted[1]  # Tau value

      else:
          raise RuntimeError("%s is not yet implemented"%output.mode)

      return result

### BDS made this on 01/11/2017 to streamline reducing files. 
def threeDDataReducer(filePath,temp,percents,loadedData,fullDataLimit,reducedDataLimit):

        case = empty()

        num_percents = len(percents)
        diffCai = np.zeros((num_percents, num_percents))
        diffCaSR = np.zeros((num_percents, num_percents))
        maxCaSR_all = np.zeros((num_percents, num_percents))
        minCaSR_all = np.zeros((num_percents, num_percents))
        maxNai_all = np.zeros((num_percents, num_percents))

        limitsFull = np.array(fullDataLimit,dtype=int) # fullData
        limits = np.array(reducedDataLimit,dtype=int) # reduced Data

        for (i ,percent_i) in enumerate(percents):
                for (j, percent_j) in enumerate(percents):

                        print ""
                        print "#######"
                        print "Trying new file!!!"
                        print ""
                        
                        if loadedData == "NKA_mouse":
                                fileName = "mouse_BASELINE_Temp_%sp00_leak%spct_nka%spct_SERCA100p00pct_freq1p0Hz_cat.pkl"%(temp,percent_i,percent_j)
                        elif loadedData == "SERCA_mouse":
                                fileName = "mouse_BASELINE_Temp_%sp00_leak%spct_nka100p00pct_SERCA%spct_freq1p0Hz_cat.pkl"%(temp,percent_i,percent_j)
                        elif loadedData == "NKA_rat":
                                fileName = "rat_BASELINE_Temp_%sp00_leak%spct_nka%spct_SERCA100p00pct_freq1p0Hz_cat.pkl"%(temp,percent_i,percent_j)
                        elif loadedData == "SERCA_rat":
                                fileName = "rat_BASELINE_Temp_%sp00_leak%spct_nka100p00pct_SERCA%spct_freq1p0Hz_cat.pkl"%(temp,percent_i,percent_j)
                        else:
                                raise RuntimeError("Not supported at this time. Dying here. :(")

                        reducedFile = fileName.replace(".pkl","_red.pkl")

                        fileName = filePath + fileName
                        #print fileName
                        if os.path.isfile(filePath + reducedFile):
                            print "File already downsampled!! :)"
                        else:
                            try:
                                dsp.downsample(fileName,10)
                            except IOError:
                                print "Can't find file. Hopefully you already downsampled it."
	return

### BDS made this on 01/11/2017 to load data for 3D plots.
def threeDDataLoader(filePath,temp,percents_one,percents_two,loadedData,fullDataLimit,reducedDataLimit):

        case = empty()

        num_percents_one = len(percents_one)
	num_percents_two = len(percents_two)
        diffCai = np.zeros((num_percents_one, num_percents_two))
        diffCaSR = np.zeros((num_percents_one, num_percents_two))
        maxCaSR_all = np.zeros((num_percents_one, num_percents_two))
        minCaSR_all = np.zeros((num_percents_one, num_percents_two))
        maxNai_all = np.zeros((num_percents_one, num_percents_two))

        limitsFull = np.array(fullDataLimit,dtype=int) # fullData
        limits = np.array(reducedDataLimit,dtype=int) # reduced Data

        for (i ,percent_i) in enumerate(percents_one):
                for (j, percent_j) in enumerate(percents_two):

                        if loadedData == "NKA_mouse":
                                fileName = "mouse_BASELINE_Temp_%sp00_leak%spct_nka%spct_SERCA100p00pct_freq1p0Hz_cat.pkl"%(temp,percent_i,percent_j)
                        elif loadedData == "SERCA_mouse":
                                fileName = "mouse_BASELINE_Temp_%sp00_leak%spct_nka100p00pct_SERCA%spct_freq1p0Hz_cat.pkl"%(temp,percent_i,percent_j)
                        elif loadedData == "NKA_rat":
                                fileName = "rat_BASELINE_Temp_%sp00_leak%spct_nka%spct_SERCA100p00pct_freq1p0Hz_cat.pkl"%(temp,percent_i,percent_j)
                        elif loadedData == "SERCA_rat":
                                fileName = "rat_BASELINE_Temp_%sp00_leak%spct_nka100p00pct_SERCA%spct_freq1p0Hz_cat.pkl"%(temp,percent_i,percent_j)
                        else:
                                raise RuntimeError("Not supported at this time. Dying here. :(")

                        reducedFile = fileName.replace(".pkl","_red.pkl")

                        case.fileName = filePath + reducedFile
                        try:
                            case.data = anG.readPickle(case.fileName,verbose=False)
                            print case.fileName

                            s = case.data['s']
                            s_idx = case.data['s_idx']

                            ## Ca sha-izzle
                            idx = s_idx.index("Cai")
                            maxCai = np.amax(s[limits[0]:limits[1],idx])
                            #print maxCai
                            minCai = np.amin(s[limits[0]:limits[1],idx])
                            #print minCai 
                            diffCai[i,j] = maxCai - minCai
                            #print diffCai
                            if math.isnan(diffCai[i,j]):
                                diffCai[i,j] = 0.0
                                #print diffCai

                            s = case.data['s']
                            s_idx = case.data['s_idx']
                            idx = s_idx.index("Ca_SR")
                            maxCaSR = np.amax(s[limits[0]:limits[1],idx])
                            #print maxCai
                            minCaSR = np.amin(s[limits[0]:limits[1],idx])
                            #print minCai
                            maxCaSR_all[i,j] = maxCaSR
                            if math.isnan(maxCaSR_all[i,j]):
                                maxCaSR_all[i,j] = 0.0
                            minCaSR_all[i,j] = minCaSR
                            if math.isnan(minCaSR_all[i,j]):
                                minCaSR_all[i,j] = 0.0
                            diffCaSR[i,j] = maxCaSR - minCaSR
                            if math.isnan(diffCaSR[i,j]):
                                diffCaSR[i,j] = 0.0

                            ## Na ma-pizzle
                            idx = s_idx.index("Nai")
                            maxNai = np.amax(s[limits[0]:limits[1],idx])
                            #print maxCai
                            maxNai_all[i,j] = maxNai
                            if math.isnan(maxNai_all[i,j]):
                                maxNai_all[i,j] = 0.0
                        
                        except IOError:
                            print "Can't find file: " + case.fileName + " Generating placeholder data."
                            print ":( :( :( :( :("
                            diffCai[i,j] = 0.0
                            diffCaSR[i,j] = 0.0
                            maxCaSR_all[i,j] = 0.0
                            maxNai_all[i,j] = 0.0
                            
                        print "######"
                        print ""
        #print "diffCai: ", diffCai
        #print "diffCaSR: ", diffCaSR
        #print "maxCaSR_all: ", maxCaSR_all
        #print "maxNai_all: ", maxNai_all
        
        return diffCai, diffCaSR, maxCaSR_all, minCaSR_all, maxNai_all

###
### Mostly plotting 
### 
def analyOut(data1,state="Cai",label=""):
  si = data1['s']
  pi = data1['p']

  np.shape(si)
  temp = si[:,runner.model.state_indices(state)]
  # take last third
  idx = np.int(np.shape(temp)[0]*2/3)
  #print idx
  ts = data1['t']
  ts = ts[idx:]
  temp = temp[(1+idx):]
  plt.plot(ts,temp*mM_to_uM,label=label)
    
  pi = data1['p']
  PCa = pi[runner.model.parameter_indices("PCa")]
  ks = pi[runner.model.parameter_indices("ks")]
  minCai = np.min(temp)
  maxCai = np.max(temp)  
    
  return (PCa,ks,minCai,maxCai)
    
#name =  "/tmp/PCa1.00ks1.00.pickle"      
#d = readPickle(name)
#pca,ks,minCai,maxCai = analyOut(d)

def ProcessOneDOutputs(var1Name,names,allVars,state="Cai",xlim=None,ylim=None,offsetMin=False):
  print "WARNING: does not include time steps" 
  for i,name in enumerate(names):               
      print name
      d = anG.readPickle(name+".pickle")
      print np.shape(d['s'])
      s = d['s']
      si = s[:,runner.model.state_indices(state)]         
     
      # recenter each to minimum
      if offsetMin:
        # assume last third is in steady state
        inds  = np.shape(si)[0]
        inds = np.int(inds/3.)
        si = si - np.min(si[-inds:])
          
      plt.plot(si*mM_to_uM,label="%s=%4.2f"%(var1Name,allVars[i]))  

  if offsetMin:
    plt.title("Minima offset to 0. uM")

  plt.legend(loc=2)

  if xlim!=None:
    plt.xlim(xlim)
  if ylim!=None:
    plt.ylim(ylim)

  plt.ylabel("[%s] [uM]" % state)  
  plt.xlabel("timesteps []") # t [ms]") 
  name = state+"transients%s"%(var1Name)
  plt.gcf().savefig(name+".png")

def ProcessTwoDOutputs(allKeys,allVars,state="Cai",ylims=None,stim_period=1000,nameTag=None):
  vars1 = allVars[0]
  var1Name = allKeys[0]
  vars2 = allVars[1]
  var2Name = allKeys[1]

  outsMin = np.zeros([np.shape(vars1)[0],np.shape(vars2)[0]])
  outsMax = np.zeros(outsMin.shape)

  # assuming there exist two iterated var 
  #for i, PCa in enumerate(PCas):
  for i, var1Val in enumerate(vars1):
    plt.figure()  
    plt.title("%s=%3.2f"%(var1Name,var1Val))
    #for j, ks in enumerate(kss):
    for j, var2Val in enumerate(vars2):
        name =namer(var1Name,var1Val,var2Name,var2Val,stim_period=stim_period,tag=nameTag)+".pickle"
        #print name
        try: 
          d = anG.readPickle(name) 
        except: 
          print name + " was not found. Skipping" 
          continue 
        #print np.shape(d['s'])
        dummy,dummy,minCai,maxCai = analyOut(d,state=state,label="%s=%3.2f"%(var2Name,var2Val))
        outsMin[i,j]= minCai 
        outsMax[i,j]= maxCai 
    if ylims!=None:    
      plt.ylim(ylims)  
    plt.ylabel(state+" [uM]")
    plt.xlabel("t [ms]")  
    plt.legend(loc=0,ncol=3)  
    name = state+"transients%s%3.2f"%(var1Name,var1Val)
    plt.gcf().savefig(name.replace(".","p")+".png")
    
  return outsMin,outsMax
    
# labeli - label corresponding to key[i]
def TwoDPlots(allKeys,allVars,outsMin, outsMax,label0="",label1="",state="Cai"):
    xv, yv = np.meshgrid(allVars[1],allVars[0])

    plt.figure()
    plt.subplot(2,2,1)
    plt.title(state+" min")
    plt.pcolormesh(xv,yv,outsMin*mM_to_uM)#,shading='gouraud')
    plt.xlim([np.min(allVars[1]),np.max(allVars[1])])
    plt.ylim([np.min(allVars[0]),np.max(allVars[0])])
    plt.ylabel(label0)
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.title(state+" max")
    plt.pcolormesh(xv,yv,outsMax*mM_to_uM)
    plt.xlim([np.min(allVars[1]),np.max(allVars[1])])
    plt.ylim([np.min(allVars[0]),np.max(allVars[0])])
    plt.ylabel(label0)
    plt.colorbar()
    plt.xlabel(label1)        
    
    plt.subplot(2,2,3)
    plt.title(state+" diff")
    plt.pcolormesh(xv,yv,(outsMax-outsMin)*mM_to_uM)
    plt.ylabel(label0)
    plt.xlim([np.min(allVars[1]),np.max(allVars[1])])
    plt.ylim([np.min(allVars[0]),np.max(allVars[0])])
    plt.colorbar()
    plt.xlabel(label1)         
    plt.tight_layout()
    
    name = state+"_extrema.png"
    plt.gcf().savefig(name) 
    return       

# loads data stored in ode object
# returns a single array with quantity of interest (valsIdx) 
# which is the time series of the idxName 
def GetData(data,idxName):
    raise RuntimeError("Antiquated - use analyzeGotran") 
    #print "getting data" 
    datac = empty()
    datac.t = data['t'] * ms_to_s   # can't guarantee units are in [ms]
    #datac.s = data['s'] / mM_to_uM  # ??? not all states are in uM....
    datac.s_idx = data['s_idx']
    #datac.j = data['j']
    datac.j_idx = data['j_idx']

    if idxName in datac.j_idx:
      datac.v = data['j']   
      datac.v_idx = datac.j_idx
    # states 
    elif idxName in datac.s_idx:
      datac.v = data['s'] 
      datac.v_idx = datac.s_idx
    else:
      print idxName, " not found"
      datac.v =None

    idx = datac.v_idx.index(idxName)
    datac.valsIdx = datac.v[:,idx] 

    return datac

### Below definition made by BDS on 10/13/2016 ###  
def Plot_Pickle_Data(rootOutput,datas,state=None,colors=None,xlabel=None,ylabel=None,
                   Full_image_xlim=None,Zoomed_image_xlim=None,plot_ylim=None,unit_scaler=None,
                   legends=None,time_range=None,legendMover1=None,legendMover2=None,legendPlacer=None):

    Zoomed_image = plt.subplot(1,2,2)
    Full_image = plt.subplot(1,2,1)
    File_Name_Cases = ""
    
    for i,data in enumerate(datas):
        #print state #data
        extracted_data = GetData(data,state)
        #print "here: ",extracted_data.v_idx.index(state)
        idx = extracted_data.v_idx.index(state)
        #print "value: ", extracted_data.v[:,idx]*unit_scaler
        Full_image.plot(extracted_data.t,(extracted_data.v[:,idx]*unit_scaler),colors[i],label=legends[i])
        Zoomed_image.plot(extracted_data.t,(extracted_data.v[:,idx]*unit_scaler),colors[i],label=legends[i])
        File_Name_Cases += "%s_" %legends[i]

    #plt.locator_params(nbins=6)
    Full_image.locator_params(nbins=8)
    Zoomed_image.locator_params(nbins=8)
    Full_image.set_xlim(Full_image_xlim)
    Zoomed_image.set_xlim(Zoomed_image_xlim)
    Full_image.set_ylim(plot_ylim)
    Zoomed_image.set_ylim(plot_ylim)
    
    plt.xlabel(xlabel, weight="bold",fontsize=14)
    plt.ylabel(ylabel, weight="bold",fontsize=14)
    plt.tight_layout()
   
    Zoomed_image.locator_params(axis='x',nbins=6) 
    Zoomed_image.get_xaxis().get_major_formatter().set_useOffset(False)
    
    #Full_image.plot([-1,300],[0.00,0.00],"k")   
    #Zoomed_image.plot([-1,300],[0.00,0.00],"k")
 
    outFile = rootOutput+"Intracellular_%s_%splots.png"%(state,File_Name_Cases)
    print outFile
    
    if legendPlacer == state or legendPlacer == None:
        
        art = []
        lgd = plt.legend(bbox_to_anchor=(legendMover1,legendMover2))
        art.append(lgd)

        plt.gcf().savefig(outFile,additional_artists=art,bbox_inches='tight',dpi=600)
    else:
        plt.gcf().savefig(outFile,dpi=600)
        
    plt.show()
    plt.close()

def Plot_Pickle_Data_One_Plot(rootOutput,datas,state=None,colors=None,xlabel=None,ylabel=None,
                   Full_image_xlim=None,Zoomed_image_xlim=None,plot_ylim=None,unit_scaler=None,
                   legends=None,time_range=None,legendMover1=None,legendMover2=None,legendPlacer=None):

    File_Name_Cases = ""
    
    for i,data in enumerate(datas):
        extracted_data = GetData(data,state)
        print extracted_data
        idx = extracted_data.v_idx.index(state)
        plt.plot(extracted_data.t,(extracted_data.v[:,idx]*unit_scaler),colors[i],label=legends[i])
        File_Name_Cases += "%s_" %legends[i]

    #plt.locator_params(nbins=6)
    plt.locator_params(nbins=8)
    plt.xlim(Full_image_xlim)
    plt.ylim(plot_ylim)
    
    plt.xlabel(xlabel, weight="bold",fontsize=14)
    plt.ylabel(ylabel, weight="bold",fontsize=14)
    plt.tight_layout()
    
    ax = plt.gca()
    ax.get_xaxis().get_major_formatter().set_useOffset(False)

    #Full_image.legend(loc=3)
    outFile = rootOutput+"Intracellular_%s_%splots.png"%(state,File_Name_Cases)
    print outFile
    
    if legendPlacer == state or legendPlacer == None:
        
        art = []
        lgd = plt.legend(bbox_to_anchor=(legendMover1,legendMover2))
        art.append(lgd)

        plt.gcf().savefig(outFile,additional_artists=art,bbox_inches='tight',dpi=300)
    else:
        plt.gcf().savefig(outFile,dpi=300)
        
    plt.show()
    plt.close()


def PlotFluxes(t,j,idx1=None,idx1Name="i_Ca",label1="flux1",idx2=None,label2=None):      

  raise RuntimeError("No longer using this. Try/revise PlotPickle()")
  if idx1==None:
    idx1=runner.model.monitor_indices(idx1Name)    
  #if idx1==None:
  #  idx1=runner.model.monitor_indices(label1)      

  fig, ax1 = plt.subplots()
  rects = []
  labels= []
  #print "WARNING: should pull Cm, Vol, F from shannon model"
  #i_to_j = 2e-2 # [A/F] --> [uM/ms]
  labels.append(label1)
  #rect1 =ax1.plot(t,j[:,idx1],'k',label=label1)
  rect1 =ax1.plot(t,j[idx1,:],'k',label=label1)
  rects.append(rect1[0])
  ax1.set_xlabel('time [ms]')
  ax1.set_ylabel(label1)
  
  
  if idx2!= None:
    ax2 = ax1.twinx()
    labels.append(label2)
    rect2 = ax2.plot(t,j[:,idx2],'k--',label=label2)
    rects.append(rect2[0])
    ax2.set_ylabel(label2)
    leg = ax2.legend( (rects) , (labels), loc=2, fancybox=True)
  else: 
    leg = ax1.legend( (rects) , (labels), loc=2, fancybox=True)
  
  
  
  leg.get_frame().set_alpha(1.0) 

#idxCai = runner.model.state_indices("Cai")
#idxCaSR = runner.model.state_indices("Ca_SR")
#idxjRyR = runner.model.monitor_indices("j_rel_SR")


def TransientBarPlots(cases,caseNames,resultsA,resultsB=False,tag=""):
    plt.subplot(1,3,1)
    idxs = np.arange(len(cases))
    width = 0.5
    plt.bar(idxs,resultsA.pctChgs,width)
    if resultsB:
      plt.bar(idxs+width,resultsB.pctChgs,width,color="r")
    plt.title("$\Delta Ca_i^{2+}$")
    plt.ylabel('$\Delta$ Ca [uM]')
    plt.xticks(idxs+0.5*width, (caseNames),rotation=70)

    plotAsDelta=False
   
    plt.subplot(1,3,2)
    if plotAsDelta:
      plt.bar(idxs,resultsA.pctChgSRs,width)
      if resultsB:
       plt.bar(idxs+width,resultsB.pctChgSRs,width,color="r")
      plt.title("$\Delta Ca_{SR}^{2+}$")
      plt.ylabel('$\Delta$ Ca [uM]')
      plt.xticks(idxs+0.5*width, (caseNames),rotation=70)
    else:
      plt.bar(idxs,resultsA.minCaSRs,width,color='b')
      plt.bar(idxs,resultsA.pctChgSRs,width,color='g',bottom=resultsA.minCaSRs)
      #print resultsA.maxCaSRs
      #if resultsB:
      # plt.bar(idxs+width,resultsB.pctChgSRs,width,color="r")
      plt.title("$Ca_{SR}^{2+}$")
      plt.ylabel("Min/$\Delta Ca^{2+}$ [uM]")
      plt.xticks(idxs+0.5*width, (caseNames),rotation=70)


    ax = plt.subplot(1,3,3)
    rects1 = ax.bar(idxs,resultsA.taus,width)
    if resultsB:
      rects2 = ax.bar(idxs+width,resultsB.taus,width,color="r")
    plt.title("Decay constant")
    plt.ylabel('tau [ms]')
    plt.xticks(idxs+0.5*width, (caseNames),rotation=70)
    #plt.legend()
    plt.ylim([0,300])
    #label1="%3.1f Hz"%(1000/A.)
    #label2="%3.1f Hz"%(1000/B.)
    #ax.legend((rects1[0], rects2[0]), (label1,label2),loc=2)

    plt.tight_layout()
    plt.gcf().savefig(tag+"miscdata.png",dpi=300)

    plt.figure()
    plt.subplot(1,3,1)
    plt.bar(idxs,resultsA.maxRyRs,width)
    if resultsB:
      plt.bar(idxs+width,resultsB.maxRyRs,width,color="r")
    plt.title("$RyR Max$")
    plt.ylabel('$jRyR$ [uM/ms]')
    plt.xticks(idxs+0.5*width, (caseNames),rotation=70)


    plt.tight_layout()
    plt.gcf().savefig(tag+"jRyRmiscdata.png",dpi=300)

# Collect all data 
def ProcessAllTransients(cases,caseTags,pacingInterval,tag="",\
                         cols=[],name=None,root=""):
    baseline = anG.readPickle(root+caseTags[0]+tag+".pickle")
    caseA = anG.readPickle(root+caseTags[1]+tag+".pickle")
    caseB= anG.readPickle(root+caseTags[2]+tag+".pickle")
    #pca1p25vmax1p25 = readOut("PCa1.25ks1.00vMax1.25"+tag+".pickle")
    if len(cases)>3:
      caseC = anG.readPickle(root+caseTags[3]+tag+".pickle")

    
    taus=[]
    pctChgs=[]
    pctChgSRs=[]
    maxRyRs=[]
    minCaSRs=[]
    maxCaSRs=[]
    for i, case in enumerate(cases):
      var = eval(case)
      tau,pctChg,pctChgSR,minCaSR,maxCaSR,maxRyR  = ProcessTransients(var,pacingInterval)
      print minCaSR, maxCaSR
      taus.append(tau)
      pctChgs.append(pctChg)
      pctChgSRs.append(pctChgSR)  
      maxRyRs.append(maxRyR)    
      minCaSRs.append(minCaSR)    
      maxCaSRs.append(maxCaSR)    

    plt.figure()    
    plt.subplot(1,2,1) 
    for i, case in enumerate(cases): 
      var = eval(case)
      s = var['s']
      t = var['t']
      cai = s[1:,idxCai]
      idx = 8000
      plt.plot(t[idx:]*ms_to_s,cai[idx:]*mM_to_uM,cols[i],label=case)
      plt.xlabel("t [ms]")
      plt.ylabel("Ca [uM]")          
    #plt.legend(loc=3,ncol=2)
    plt.title("Cai")
    plt.ylim([0,1.2]) 
    
    plt.subplot(1,2,2)    
    for i, case in enumerate(cases):
      var = eval(case)  
      s = var['s']
      t = var['t']
      casr = s[1:,idxCaSR]
      idx = 8000
      plt.plot(t[idx:]*ms_to_s,casr[idx:]*mM_to_uM,cols[i],label=case)
      plt.xlabel("t [ms]")
      plt.ylabel("Ca [uM]")
    plt.ylim([250,1e3])    
    plt.legend(loc=2,ncol=1)
    plt.title("CaSR")

    plt.tight_layout()
    plt.gcf().savefig(name,dpi=300)

    plt.figure()
    for i, case in enumerate(cases):
      var = eval(case)  
      j = var['j']
      t = var['t']
      jRyR = j[:,idxjRyR]
      idx = 8000
      plt.plot(t[idx:]*ms_to_s,jRyR[idx:]*mM_to_uM,cols[i],label=case)
      plt.xlabel("t [ms]")
      plt.ylabel("j [uM/ms]")
    plt.legend(loc=0,ncol=1)      
    plt.ylim([0,200])    
    plt.xlim([8.5,8.7])    
    plt.title("jRyR")
    plt.tight_layout()
    plt.gcf().savefig("jRyR"+name,dpi=300)
    
    
    
    results = empty()
    results.taus = taus; 
    results.pctChgs = pctChgs; 
    results.pctChgSRs = pctChgSRs
    results.maxRyRs = maxRyRs
    results.maxCaSRs = maxCaSRs
    results.minCaSRs = minCaSRs
    
    return results


# for despa data 
def DespaPlots(caseDict, stateName=None, monitorName=None,doLegend=True,label="",loadPickle=False,\
               xlim=[7000,10000],ylim=None,
               root=None,normalize=False):
    if loadPickle:
      LoadPickles(caseDicts)
#print "python runShannonTest.py "+" ".join(healthyArgs)+" -name "+healthyName+" &"

    if ylim:
      ylimUpdate=False
    else:
      ylim=[1e9,-1e9]
      ylimUpdate=True
    for key, case in caseDict.iteritems():
      t = case.data['t']
      j = case.data['j']
      s = case.data['s']
       
      if stateName!=None:  
        idxCai = runner.model.state_indices(stateName)
        #sca = s[1:,idxCai] - np.min(s[-ind:,idxCai])
        sca = s[1:,idxCai]  
      elif monitorName!=None:
        idxCai = runner.model.monitor_indices(monitorName)
        #sca = s[1:,idxCai] - np.min(s[-ind:,idxCai])
        sca = j[:,idxCai]  

      print "%s Diast/Systolic %f/%f " % \
            (case.tag,np.min(sca[xlim[0]:xlim[1]]),np.max(sca[xlim[0]:xlim[1]]))

      
      if normalize:
        sca = sca - np.min(sca[xlim[0]:xlim[1]])
        sca = sca/np.max(sca[xlim[0]:xlim[1]])
      plt.plot(t,sca,label=case.label)
            
      # determine bounds   
      if ylimUpdate:
        ylim[0] = np.min([ylim[0],np.min(sca)])
        ylim[1] = np.max([ylim[1],np.max(sca)])

    #plt.title("%s transients (offset by diastolic [%s])" %(stateName,stateName))    
   # plt.ylabel("[%s] - min([%s])" %(stateName,stateName))
    if stateName!=None:  
      plt.title("%s transients" %(stateName))    
      plt.ylabel("[%s] [mM])" %(stateName))
      fileName = stateName  
    elif monitorName!=None:
      plt.title("%s " %(monitorName))    
      plt.ylabel(label)
      fileName = monitorName
        
    if doLegend:
      plt.legend(loc=0)    
    plt.xlabel("t [ms]")
    if xlim!=None:
      plt.xlim(xlim)
    plt.ylim(ylim)
    if root==None:
      root=""
  
    if normalize:
      normalize="Norm"
    else:
      normalize=""        
    plt.gcf().savefig(root+"despa_%s%s.png"%(fileName,normalize),dpi=300)

###
### Data processing 
###

## Compute quantities of interest from transient data
# transient decay rates, 
# changes in SR calcium 
# ryr flux rates 
def ProcessTransients(case,pacingInterval,tstart=8000):
    # compute rate of \catwo transient decline 
    tau = tf.GetTau(case,pacingInterval,tstart=tstart,idxCai=idxCai)
    
    
    # get transient amplitudes 
    # grab suitable spot for statistics 
    tsub, caisub= tf.GetInterval(case,pacingInterval,tstart=tstart,idx=idxCai)
    delCa = (np.max(caisub) - np.min(caisub))*mM_to_uM
    
    # get SR transient amplitudes 
    tsub, casrsub= tf.GetInterval(case,pacingInterval,tstart=tstart,idx=idxCaSR)
    delCaSR = (np.max(casrsub) - np.min(casrsub))*mM_to_uM
    minCaSRsub = (np.min(casrsub))*mM_to_uM
    maxCaSRsub = (np.max(casrsub))*mM_to_uM
    
    # get max RyR
    tsub, jRyRsub= tf.GetInterval(case,pacingInterval,tstart=tstart,getFlux=True, idx=idxjRyR)
    #plt.figure()
    #plt.plot(jRyRsub)
    #jRyR = j[:, runner.model.monitor_indices(idx) ]
    maxRyR = np.max(jRyRsub)*mM_to_uM
    
    print tau, delCa, delCaSR, maxRyR
    #pctChg = GetExtreme(t,cai,subMin=7900,subMax=8500,si=7910)
    #pctChgSR = GetExtreme(t,casr,subMin=7900,subMax=8500,si=7910)
    
    return tau,delCa, delCaSR, minCaSRsub,maxCaSRsub,maxRyR

## Power spectral density stuff 
        
def doPSD(sca):
  dm = sca - np.mean(sca,axis=0)
  M = fftp.fft(dm,axis=0)
  psd = np.abs(M)**2
  return psd  


#psd = doPSD(sca)
#plt.plot(psd)

def PlotFrequencies(s1,statei=0):
    # real domain 
    plt.figure()
    plt.pcolormesh(s1.T)
    plt.figure()
    plt.plot(s1[:,statei],label="%d"%statei)
    plt.legend()
    
    # get PSD
 #   dm = s1 - np.mean(s1)
 #   S = fftp.fft(dm,axis=0)
 #   psd2 = np.abs(S*S)
    
    psd2 = doPSD(s1)
    
    plt.figure()
    plt.pcolormesh(psd2.T)
    plt.figure()

    energy = np.sum(psd2,axis=0)
    psdScaled = psd2/energy
#plt.plot(psd2[:,0])
    plt.plot(psdScaled[:,statei],label="%d"%statei)

    return psdScaled

def PSDAnaly(s1,ranger=[2,200],verbose=True):
    ## raw data 
    #plt.figure()
    #plt.plot(s1)
    
    #plt.figure()
    #pcolormesh(s1.T,cmap=cm.gray)
    #plt.colorbar()
    
    ## Demeaned (looks for zeros) 
    dc = np.mean(s1,axis=0)
    dm = s1 - dc

    # abs. value 
    daMax = np.max(s1,axis=0)
    daMaxAbs = np.abs(np.max(s1,axis=0))
    daMin = np.min(s1,axis=0)
    daMinAbs = np.abs(np.min(s1,axis=0))
    #amp = daMax-daMin
    # compare abs value of min/max and see which corresponds to the largest value 
    # columnwise stack the minabs/maxabs, then get highest of the two numbers
    amp = np.max( np.column_stack([daMinAbs,daMaxAbs]), axis=1)
    #print amp
    eps = 1e-9
    nonZero = np.argwhere(np.abs(daMax)>eps)
    dm[:,nonZero] = dm[:,nonZero]/daMax[nonZero] 
    isZero = np.argwhere(np.abs(daMax)<=eps)
    dm[:,isZero] = eps

    if verbose: 
      plt.figure()
      plt.bar(np.arange(np.shape(dc)[0]),dc)
      plt.legend()
    #plt.figure()
    #pcolormesh(dm.T,cmap=cm.gray)
    #plt.colorbar()
    
    # Power Spectral Density 
    S = fftp.fft(dm,axis=0)
    psd2 = np.abs(S*S)
    # log to make peaks more observable 
    psd2[ np.argwhere( np.abs(psd2) < eps )  ]= eps
    psd2 = np.log(psd2)
    #print np.shape(psd2)
    # grab low freq 
    psd2 = psd2[ranger[0]:ranger[1],:]
    
    if verbose:
      plt.figure()
      from matplotlib import cm 
      plt.pcolormesh(psd2.T,cmap=cm.gray)
      plt.colorbar()
      plt.legend()
    
    
    #plt.figure()
    #plt.plot(psd2[ranger[0]:ranger[1],0])
    return daMin,daMax,amp,dc, psd2

# State Decomposition Analysis (SDA)
# indSS = 2e3 # collect statistics after this time point [ms] (looking for steady state)
def donorm(subj,ref):
    eps = 1e-9
    # means 
    nonZero = np.argwhere(np.abs(ref.dc)>eps)
    subj.dcn = np.zeros( np.shape(ref.dc) )
    subj.dcn[ nonZero ] = subj.dc[ nonZero ] / ref.dc[ nonZero ]

    #maxima
    nonZero = np.argwhere(np.abs(ref.max)>eps)
    subj.maxn = np.zeros( np.shape(ref.max) )
    subj.maxn[ nonZero ] = subj.max[ nonZero ] / ref.max[ nonZero ]

    #maxima
    nonZero = np.argwhere(np.abs(ref.min)>eps)
    subj.minn = np.zeros( np.shape(ref.min) )
    subj.minn[ nonZero ] = subj.min[ nonZero ] / ref.min[ nonZero ]

    #maxima
    nonZero = np.argwhere(np.abs(ref.amp)>eps)
    subj.ampn = np.zeros( np.shape(ref.amp) )
    subj.ampn[ nonZero ] = subj.amp[ nonZero ] / ref.amp[ nonZero ]

    #maxdiff
    nonZero = np.argwhere(np.abs(ref.amp)>eps)
    subj.maxdiff= np.zeros( np.shape(ref.amp) )
    subj.maxdiff[ nonZero ] = subj.maxdiff[ nonZero ] - ref.maxdiff[ nonZero ]

def Autolabel(rects,labels,ax):
# attach some text labels
	for rect,label in zip(rects,labels):
        	height = rect.get_height()
        	ax.text(rect.get_x() + rect.get_width()/2., 1.15*height,
                	"%4.1f [pA/pF]"%label,
                	rotation='vertical',
                	ha='center', va='bottom')


def StateDecompositionAnalysisBetter(rootOutput,caseDict,
                  wanted,
                  indSS=[2e3,-1], # Time aftter which data is used for comparative analysis. -1 signifies getting last time step  [ms]
                  xlim=None,
                  root="./",
                  ranked=20,
                  ignoreList = ["dNa_SL_buf"],
                  Title=None,
                  mode="states", # fluxes
                  cols = ["r","b"],
                  doPlot=False,
                  sortby="mean", # max, min, amp
                  legendMover1=None,
                  legendMover2=None,
                  j_to_i = 26.92, # convert j into i, done for mouse model
                  tag="",  # optional tag
                  autolabel = False
		 ):

    # decide on which data to pull
    if mode=="states":
        print "Pulling states"
        v_key = 's'
    elif mode=="fluxes":
        print "Pulling fluxes"
        v_key = 'j'
    else:
        raise RuntimeError(mode +" not understood")

    # comparing cases  
    subCaseDict = dict()

    for idx, wantedi in enumerate(wanted):
        for key, case in caseDict.iteritems():
            if key == wantedi:
                print "Selecting: ", key
                case.idx = idx
                subCaseDict[case.tag] = case

                t = case.data['t']
                ref= case.data[v_key]
                vcp = np.copy(ref)
                # Rescale j by j_to_i conversion factor (uM/ms --> pA/pF)
                # do NOT do this (it's wrong) 
                # INEFFICIENT, BUT NEED TO WRAP THIS UP 
                #fluxIdxs = case.data[v_key+"_idx"]
                #for i in np.arange(np.shape(vcp)[1]):
                #   if "j_" in fluxIdxs[i]:
                #     vcp[:,i] = vcp[:,i] # * j_to_i
                #     print "Converting ",fluxIdxs[i]," [uM/ms] --> [pA/pF] FOR MICE!!"

                sub = vcp[int(indSS[0]):int(indSS[1])]
                case.min,case.max, case.amp,case.dc, case.psd2 = PSDAnaly(sub,verbose=False)
                wanted[case.idx] = case
                donorm(wanted[idx],wanted[0])

    ## displaycomparison
    # get shortest traj
    lastT = 1e30
    for key, case in subCaseDict.iteritems():
        ti = case.data['t']
        # pick whichever T is smallest: max length of either trajector or the xlim bound
        if xlim != None:
            lastT = np.min([ti[-1],lastT,xlim[1]])
        else:
            lastT = np.min([ti[-1],lastT])
    if xlim==None:
        xlim = [0,lastT]

    # sort from largest to smallest change 
    if sortby == "mean" or sortby == "dc":
        valueChg = wanted[1].dcn-wanted[0].dcn
    elif sortby == "max":
        valueChg = wanted[1].maxn-wanted[0].maxn
    elif sortby == "min":
        valueChg = wanted[1].minn-wanted[0].minn
    elif sortby == "amp":
        valueChg = wanted[1].amp-wanted[0].amp
    sort_index = (np.argsort(np.abs(valueChg)))[::-1]
    # verify we're getting the right numbers 
    #print valueChg[sort_index] #[0:5]]
    daIdx=97

    #Assuming that both pickle files have same states/ode model. 
    v_idx = case.data['%s_idx'%v_key]
    value_inds = {key: idx for (idx, key) in enumerate(v_idx)}
    # create reverse lookup
    value_inds_rev = {idx: key for (idx, key) in enumerate(v_idx)}
    #print value_inds_rev[ sort_index[daIdx] ]

    # grabbing top-twenty modulated states
    bestidx    = []
    bestvalues = []
    lines = []
    lines.append("name caseno val diff_val\n")

    stored = 0
    for i,idx in enumerate(sort_index):
        if idx not in value_inds_rev:
            raise ValueError("Unknown state/flux: '{0}'".format(idx))

        if value_inds_rev[idx] in ignoreList:
            #print "Skipping %d/%d:"%(i,idx), value_inds_rev[idx]
            continue

        # print raw values
        if sortby == "mean":
            line=value_inds_rev[idx]+" pct %4.2f "%valueChg[idx]+ \
                                      " 0 %4.1e %4.1e/%4.1e"% (wanted[0].max[idx], wanted[0].dc[idx],wanted[0].dcn[idx])+\
                                      " 1 %4.1e %4.1e/%4.1e"% (wanted[1].max[idx], wanted[1].dc[idx],wanted[1].dcn[idx])
        elif sortby == "max":
            line = "MAX %-20s"%value_inds_rev[idx]+ " 0 %4.1e 1 %4.1e [%3.1f] "%( wanted[0].max[idx],
                                                                                 wanted[1].max[idx],100*wanted[1].max[idx]/wanted[0].max[idx]-100)
        elif sortby == "min":
            line = "MEAN %-20s"%value_inds_rev[idx]+ " 0 %4.1e 1 %4.1e [%3.1f] "%(wanted[0].min[idx],
                                                                                  wanted[1].min[idx],100*wanted[1].min[idx]/wanted[0].min[idx]-100)
        elif sortby == "amp":
            line = "AMP  %-20s"%value_inds_rev[idx]+ " 0 %4.1e 1 %4.1e [%4.2f,pctdiff %3.1f] "%(wanted[0].amp[idx],
                                                                                  wanted[1].amp[idx],
                                                                                  wanted[1].amp[idx]-wanted[0].amp[idx],
                                                                                  100*wanted[1].amp[idx]/wanted[0].amp[idx]-100)
        print "idx(%d/%d): "%(i,idx),line

        # store 
        lines.append(line)
        bestidx.append(idx)
        bestvalues.append(value_inds_rev[idx])
        #print beststates[i]
        #indices.append(state_inds[state])

        stored += 1
        #print "Stored: ", stored
        #print "Ranked: ", ranked
        if stored >= ranked and ranked > 0:
            print "Pulling the rip cord!!!!!"
            break

    ## Plot State Data 
    if mode == "states":
        plotValues = ["V","Cai","Ca_SR"] + bestvalues
    elif mode =="fluxes":
        plotValues = ["i_Cab"] + bestvalues
    if doPlot:
        for label in plotValues:
            PlotValueComparison(subCaseDict,v_key,label,  cols = ["k","b","r","g"],xmin=indSS[0])

    ## Plot comparative data 
    width = 0.2
    plt.figure()
    fig, ax = plt.subplots()
    ind = np.arange(stored)

    # store info for later analysis
    wanted[0].bestidx =  bestidx
    wanted[0].bestvalues =  bestvalues
    wanted[1].bestidx = bestidx
    wanted[1].bestvalues = bestvalues

    vals0s = wanted[0].dcn
    vals1s = wanted[1].dcn
    vals2s = wanted[2].dcn
    vals3s = wanted[3].dcn

    # bar plot 
    rects1 = ax.bar(ind, 100*vals0s[bestidx], width,color='k')
    rects2 = ax.bar(ind+width, 100*vals1s[bestidx], width,color='b')
    rects3 = ax.bar(ind+(2*width), 100*vals2s[bestidx], width,color='r')
    rects4 = ax.bar(ind+(3*width), 100*vals3s[bestidx], width,color='g')

    # autolabel
    if autolabel == True:
    	labels = wanted[0].amp[bestidx]
    	Autolabel(rects1,labels,ax)

    ax.set_xticks(ind+width)
    ax.set_xticklabels(bestvalues,rotation=90)

    lb1 = wanted[0].label
    lb2 = wanted[1].label
    lb3 = wanted[2].label
    lb4 = wanted[3].label

    plt.title(Title)
    ax.set_ylabel("% of Control")
    plt.tight_layout()

    art = []
    lgd = plt.legend( (rects1[0], rects2[0],rects3[0],rects4[0]), (lb1,lb2,lb3,lb4),bbox_to_anchor=(legendMover1,legendMover2))
    art.append(lgd)

    #plt.gcf().savefig(root+versionPrefix+"comparative.png",dpi=300)
    outFile = rootOutput+"comparative_%s_%s_%s_%s_%s_%s"%(mode,wanted[0].label,wanted[1].label,wanted[2].label,wanted[3].label,sortby)
    outFile+=tag
    plt.gcf().savefig(outFile+".png",additional_artists=art, bbox_inches="tight",dpi=1200)

    f = open(outFile+'.txt', 'w')
    json.dump(lines, f)
    f.close()

    return wanted

def PlotMorotti(cases,
                case1Name='rabbit_5',case2Name='mouse_5',
                trange=[2.0e3,2.2e3],
                root ="./" # path for prionting figures 
                ):
    case1=cases[case1Name]
    case2=cases[case2Name]    
    states = ["V", "Cai", "Nai"]
    ctr=0
    for i,state in enumerate(states):
        PlotPickleData(case1.data,data2=case2.data,idxName=state, 
                          ylabel=state,trange=trange,
                          case1legend=case1.caseName,
                          case2legend=case2.caseName)
        plt.tight_layout()
        #title = case1.caseName+"_mouserabbit_comp"+state
        title = "mouse_rabbit_compare_%.2d_"%ctr+state
        plt.gcf().savefig(root+"/"+title+".png")
        ctr+=1

    fluxes = ["i_Na", "i_CaL", 
              "i_kur",# I think this is Morotti's I_K,slow
              "i_ss",
              "i_tof", "i_tos","i_Ks","i_Kr","i_K1",
              "i_NaCa","i_NaK",
              "i_Kp"]

    #I_Ca, IK,Slow, Iss, Ito, IKs, IKr, IK1, INaCa, INaK
    for i,flux in enumerate(fluxes):
        PlotPickleData(case1.data,data2=case2.data,idxName=flux,
                          ylabel="%s [A/F]"%flux,trange=trange,
                            case1legend=case1.caseName,
                          case2legend=case2.caseName)                        
        #title = case1.caseName+"_mouserabbit_comp"+flux
        title = "mouse_rabbit_compare_%.2d_"%ctr+flux
        ctr+=1
        plt.gcf().savefig(root+"/"+title+".png")
        
def shiftNnorm(ctl_case):
    minCtl_Cai_case = np.amin(ctl_case)
    Ctl_Cai_case =  ctl_case -  minCtl_Cai_case
    maxCtl_Cai_case = np.amax(Ctl_Cai_case)
    Ctl_Cai_case /= maxCtl_Cai_case
    
    return Ctl_Cai_case

def normShizzle(case0p5Hz,tRange = [290e3,291e3]):
    s = case0p5Hz.data['s']
    s_idx = case0p5Hz.data['s_idx']
    idx = s_idx.index("Cai")
    t = case0p5Hz.data['t']

    cais = s[tRange[0]:tRange[1],idx]
    min0p5Hz_Cai = np.amin(cais)
    #print minWT0p5Hz_Cai
    max0p5Hz_Cai = np.amax(cais)
    max0p5Hz_Caix = np.argmax(cais)
    print max0p5Hz_Caix

    Comp_Cai_shift0p5Hz = cais - min0p5Hz_Cai
    #print Comp_Cai_shift
    maxComp_Cai_shift0p5Hz = np.amax(Comp_Cai_shift0p5Hz)
    Comp_Cai_shift0p5Hz/=maxComp_Cai_shift0p5Hz
    #print Comp_Cai_shift

    ms_to_s = 1e-3
    time=t[tRange[0]:tRange[1]]*ms_to_s
    
    return time, Comp_Cai_shift0p5Hz

        
