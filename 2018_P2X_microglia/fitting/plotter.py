import matplotlib.pylab as plt
import analyzeGotran as aG
import numpy as np

data1=aG.readPickle("1000_cat.pickle")
#data2=aG.readPickle("100_cat.pickle")
#data3=aG.readPickle("1000_cat.pickle")

CaData1 = aG.GetData(data1,"I_ptxf")
#CaData2 = aG.GetData(data2,"Cai")
#CaData3 = aG.GetData(data3,"Cai")

Ca1 = CaData1.valsIdx - min(CaData1.valsIdx)
#Ca2 = CaData2.valsIdx - min(CaData2.valsIdx)
#Ca3 = CaData3.valsIdx - min(CaData3.valsIdx)

#I_ptxf1 = aG.GetData(data1,"I_ptxs")
#I_ptxf2 = aG.GetData(data2,"I_ptxf")
#I_ptxf3 = aG.GetData(data3,"I_ptxf")

#Ip2x41 = I_ptxf1.valsIdx #- min(CaData1.valsIdx)
#Ip2x42 = I_ptxf2.valsIdx #- min(CaData2.valsIdx)
#Ip2x43 = I_ptxf3.valsIdx #- min(CaData3.valsIdx)

I_ptxs1 = aG.GetData(data1,"I_ptxs")
#I_ptxs2 = aG.GetData(data2,"I_ptxs")
#I_ptxs3 = aG.GetData(data3,"I_ptxs")

Ip2x71 = I_ptxs1.valsIdx #-min(CaData1.valsIdx)
#Ip2x72 = I_ptxs2.valsIdx #-min(CaData2.valsIdx)
#Ip2x73 = I_ptxs3.valsIdx #-min(CaData3.valsIdx)

#plt.figure(figsize=(21,7))
#plt.subplot(1,3,1)
#plt.title('Ca transients')
#plt.plot(CaData1.t,Ca1*1000,'r-',label='10 uM')
#plt.plot(CaData2.t,Ca2*1000,'b-',label='100 uM')
#plt.plot(CaData3.t,Ca3*1000,'g-',label='1 mM')
#plt.xlim(1000,1750)
#plt.legend(loc=0)

#plt.subplot(1,3,2)
#plt.title('I-p2x4')
#plt.plot(CaData1.t,Ip2x41,'r-',label='10 uM')
#plt.plot(CaData2.t,Ip2x42,'b-',label='100 uM')
#plt.plot(CaData3.t,Ip2x43,'g-',label='1 mM')
#plt.xlim(1000,1750)
#plt.legend(loc=0)

#plt.subplot(1,3,3)
#plt.title('I-p2x7')
plt.plot(CaData1.t,Ip2x71,'r-',label='1000 uM')
#plt.plot(CaData2.t,Ip2x72,'b-',label='100 uM')
#plt.plot(CaData3.t,Ip2x73,'g-',label='1 mM')
#plt.xlim(1000,1750)
#plt.legend(loc=0)

fileName = "test.png"
plt.gcf().savefig(fileName)

print "Printed ", fileName
