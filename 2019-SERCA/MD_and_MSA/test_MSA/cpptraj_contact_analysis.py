def cpptraj_analysis(h1,h2,h3,a1,a2,a3,figname,xlim,ylim,color_title,title,x_axis,y_axis):
    import matplotlib as plt
    import numpy as np
    
    holo1 = np.zeros([xlim,xlim])
    fholo1 = np.loadtxt(h1)
    rows = np.arange( np.shape(fholo1)[0] )
    for i in rows:
        holo1[ int(fholo1[i,0]), int(fholo1[i,1]) ] = np.abs(fholo1[i,2])
    
    holo2 = np.zeros([xlim,xlim])
    fholo2 = np.loadtxt(h2)
    rows = np.arange( np.shape(fholo2)[0] )
    for i in rows:
        holo2[ int(fholo2[i,0]), int(fholo2[i,1]) ] = np.abs(fholo2[i,2])
    
    holo3 = np.zeros([xlim,xlim])
    fholo3 = np.loadtxt(h3)
    rows = np.arange( np.shape(fholo3)[0] )
    for i in rows:
        holo3[ int(fholo3[i,0]), int(fholo3[i,1]) ] = np.abs(fholo3[i,2])
    
    apo1 = np.zeros([xlim,xlim])
    fapo1 = np.loadtxt(a1)
    rows = np.arange( np.shape(fapo1)[0] )
    for i in rows:
        apo1[ int(fapo1[i,0]), int(fapo1[i,1]) ] = np.abs(fapo1[i,2])
    
    apo2 = np.zeros([xlim,xlim])
    fapo2 = np.loadtxt(a2)
    rows = np.arange( np.shape(fapo2)[0] )
    for i in rows:
        apo2[ int(fapo2[i,0]), int(fapo2[i,1]) ] = np.abs(fapo2[i,2])
    
    apo3 = np.zeros([xlim,xlim])
    fapo3 = np.loadtxt(a3)
    rows = np.arange( np.shape(fapo3)[0] )
    for i in rows:
        apo3[ int(fapo3[i,0]), int(fapo3[i,1]) ] = np.abs(fapo3[i,2])
    
    HOLOS = (holo1 + holo2 + holo3) / 3.0
    APOS = (apo1+apo2+apo3) / 3.0
    diff = HOLOS-APOS
    (nRows,nCols) = np.shape(diff)
    diffp = np.zeros([nRows,nCols])
    y = np.zeros([nRows,nCols])
    for i in range(nRows):
      diff[i,i]=0
    for i in range(nRows):
        cols = np.arange(nRows-i)+i
        for j in cols:    
          diffp[i,j] = diff[i,j]
          diffp[j,i] = diff[i,j]
    for i in range(nRows):
        cols = np.arange(nRows-i)+i
        for j in cols:
            if diffp[i,j]>0:
                y[i,j]=diffp[i,j]
            elif diff[i,j]<0:
                y[j,i]=diffp[i,j]
        
    font = {'weight' : 'bold','size'   : 10}
    plt.rc('font', **font)
    plt.pyplot.pcolormesh(y.T, cmap='seismic')
    plt.pyplot.xlim(0,xlim)
    plt.pyplot.ylim(0,ylim)
    plt.pyplot.colorbar(label = color_title)
    plt.pyplot.title(title,**font)
    plt.pyplot.xlabel(x_axis, **font)
    plt.pyplot.ylabel(y_axis, **font)
    plt.pyplot.gcf().savefig(figname, dpi=1200.)
    plt.pyplot.show()
