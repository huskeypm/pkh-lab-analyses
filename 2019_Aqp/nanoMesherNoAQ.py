import numpy as np
import math as m
nm = 1.0e-9
aqpH = 3.0e-9

def Cylinder(centers,radius=1e-9,boxSize=15e-9, length = 1e-8, revH = 1e-8):

    text = "lc = 1e-9;\n"
    text += "lc2 = 10e-10;\n"
    lineNames = ""
    text+="radius = %s;\n"%(radius)
    text+="boxSize = %s;\n"%(boxSize)
    text+="length = %s;\n"%length
    text+="revH =%s;\n"%(revH)
    text+="aqpH =%s;\n"%(aqpH)
    text+="""
Point(1) = { -boxSize, -boxSize, 0, lc};
Point(2) = { -boxSize, boxSize, 0, lc};
Point(3) = { boxSize, boxSize, 0, lc};
Point(4) = { boxSize, -boxSize, 0, lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(1) = {1,2,3,4};

Point(5) = { -boxSize, -boxSize, length, lc};
Point(6) = { -boxSize, boxSize, length, lc};
Point(7) = { boxSize, boxSize, length, lc};
Point(8) = { boxSize, -boxSize, length, lc};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
Line Loop(2) = {5,6,7,8};

Point(9) = { -boxSize, -boxSize, -revH, lc};
Point(10) = { -boxSize, boxSize, -revH, lc};
Point(11) = { boxSize, boxSize, -revH, lc};
Point(12) = { boxSize, -boxSize, -revH, lc};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,9};
Line Loop(6) = {9,10,11,12};
Plane Surface(3) = {6};

Line(17) = {9,1};
Line(18) = {4,12};
Line Loop(7) = {17,-4,18,12};
Plane Surface(4) = {7};

Line(19) = {3,11};
Line Loop(8) = {-18,-3,19,11};
Plane Surface(5) = {8};

Line(20) = {2,10};
Line Loop(9) = {-19,-2,20,10};
Plane Surface(6) = {9};

Line Loop(10) = {-20,-1,-17,9};
Plane Surface(7)={10};



Point(13) = { -boxSize, -boxSize, length+revH, lc};
Point(14) = { -boxSize, boxSize, length+revH, lc};
Point(15) = { boxSize, boxSize, length+revH, lc};
Point(16) = { boxSize, -boxSize, length+revH, lc};
Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};
Line(16) = {16,13};
Line Loop(11) = {13,14,15,16};
Plane Surface(8) = {11};

Line(21) = {13,5};
Line(22) = {8,16};
Line Loop(12) = {-21,-16,-22,8};
Plane Surface(9) = {12};

Line(23) = {15,7};
Line Loop(13) = {22,-15,23,7};
Plane Surface(10) = {13};

Line(24) = {14,6};
Line Loop(14) = {-23,-14,24,6};
Plane Surface(11) = {14};

Line Loop(15) = {-24,-13,21,5};
Plane Surface(12)={15};



"""

    #text+="Plane Surface(1) = {1};\n"
    poreTops = []
    ptID=50
    cirNum = 50
    for i,center in enumerate(centers):

        text+="Point(%s) = { %s, %s, 0, lc2};\n"%(ptID,center[0],center[1])
        ptID+=1

        text+="Point(%s) = { %s, %s, 0, lc2};\n"%(ptID,center[0]+radius,center[1])
        ptID+=1

        text+="Point(%s) = { %s, %s, 0, lc2};\n"%(ptID,center[0],center[1]-radius)
        ptID+=1

        text+="Point(%s) = { %s, %s, 0, lc2};\n"%(ptID,center[0]-radius,center[1])
        ptID+=1

        text+="Point(%s) = { %s, %s, 0, lc2};\n"%(ptID,center[0],center[1]+radius)
        ptID+=1

        text+="Circle(%s) = { %s, %s, %s};\n"%(cirNum,ptID-4,ptID-5,ptID-3)
        cirNum+=1

        text+="Circle(%s) = { %s, %s, %s};\n"%(cirNum,ptID-3,ptID-5,ptID-2)
        cirNum+=1

        text+="Circle(%s) = { %s, %s, %s};\n"%(cirNum,ptID-2,ptID-5,ptID-1)
        cirNum+=1

        text+="Circle(%s) = { %s, %s, %s};\n"%(cirNum,ptID-1,ptID-5,ptID-4)




        text+="Line Loop(%s) = {%s,%s,%s,%s};\n"%(i+20,cirNum-3,cirNum-2,cirNum-1,cirNum)
        text+="Plane Surface(%s) = {%s};\n\n"%(i+20,i+20)
        cirNum+=1

        if i==0:
            lineNames+= str(i+20)
            poreTops.append(str(i+20))
        else:
            lineNames+= "," + str(i+20)
            poreTops.append(str(i+20))
    cirNum+=1 #THIS IS A HACK TO GET extruded circle numbers
    #text+= "Plane Surface(1) = {1,%s};\n"%(lineNames)
    #text+="NewSurf1[]=Translate{0,0,%s}{ Duplicata{Surface{1};}};"%(length)
    #text+="Plane Surface(2) = {2,%s};\n"%(lineNames)


    longString=""
    shorterString=""
    listExtrudedCircleNames = ""
    volumeString=""
    poreBot = []
    for i,center in enumerate(centers):
        text+="""out%s[]=Extrude{0,0,length}{Surface{%s};};\n"""%(i+1,i+20)
        shorterString +="out%s[0],"%(i+1)
        poreBot.append("out%s[0]"%(i+1))
        volumeString+="out%s[1],"%(i+1)
        longString +="out%s[2],out%s[3],out%s[4],out%s[5],"%(i+1,i+1,i+1,i+1)
        text+="Line Loop(%s) = {%s,%s,%s,%s};\n"%(cirNum,cirNum,cirNum+1,cirNum+2,cirNum+3)
        listExtrudedCircleNames+= "%s,"%(cirNum)
        cirNum+=22
    longString=longString[:-1]
    shorterString=shorterString[:-1]
    volumeString=volumeString[:-1]
    listExtrudedCircleNames= listExtrudedCircleNames[:-1]
    ########################### TESTING TO ADD AQP


    text+= "Plane Surface(1) = {1,%s};\n"%(lineNames)
    #text+="NewSurf1[]=Translate{0,0,%s}{ Duplicata{Surface{1};}};"%(length)
    #texts= text[:-2] #this is stupid ugly
    text+= "Plane Surface(2) = {2,%s};\n"%(listExtrudedCircleNames)
    text+= "Surface Loop(1) = {1,%s,6,5,3,7,4};\n"%(lineNames)
    text+= "Volume(%s) = {1};\n"%(cirNum)
    text+= "Surface Loop(2) = {2,%s,9,10,11,8,12};\n"%(shorterString)
    text+= "Volume(%s) = {2};\n"%(cirNum+1)



    text+= "Physical Surface(1) = {1,2};\n"  ### This is the z=0 wall not pores
    #text+= "Physical Surface(2) = {%s};\n"%(lineNames) ### These are the pore bots
    text+="Physical Surface(2) = {%s};\n"%(longString)  ### These are the pore walls

    text+="Physical Surface(3) = {3};\n"  ### This is the top of the top reservoir
    text+="Physical Surface(4) = {8};\n"  ### This is the bottom of the bottom res
    text+="Physical Surface(5) = {%s};\n"%(shorterString)   ### These are the pore bots
    text+="Physical Surface(6) = {%s};\n"%(lineNames) ### These are the pore tops
    text+="Physical Volume(1) = {%s};\n"%(cirNum)
    text+="Physical Volume(2) = {%s};\n"%(cirNum+1)
    text+="Physical Volume(3) = {%s};\n"%(volumeString)





    text+="\n\n\n\n\n\n\n\n\n\n\n"
    return text#, counter-1
def getCenters(num,radius,bR=0):
    numS = num**2
    i=0
    centers = [];

    if num%2!=0:
        while i <num:
            j=0
            x = radius*(-(3.*num-4)/2+3*i-0.5)
            while j < num:

                y = radius*(-(3.*num-4)/2+3*j-0.5)
                if (abs(x)-radius/2)**2+(abs(y)-radius/2)**2<bR**2:
                    this = False
                else:
                    centers.append((x,y))
                j+=1
            i+=1

    if num%2==0:
        while i <num:
            j=0
            x = radius*(-(3.*num-4)/2+3*i-0.5)
            while j < num:
                y = radius*(-(3.*num-4)/2+3*j-0.5)
                if (abs(x)-radius/2)**2+(abs(y)-radius/2)**2<bR**2:
                    this = False
                else:
                    centers.append((x,y))
                j+=1
            i+=1
    return centers






def Build(num, radius,boxSize, length, revH,name="test"):
    bR=0
    centers = getCenters(num,radius,bR)
    text = Cylinder(centers,radius,boxSize,length,revH)
    print "post cylinder"
    #fileName = "noAQ_{}_{}".format(num,boxSize)+".geo"

    fileName = name+".geo"
    fileLines = []
    fileLines.append(text)
    with open(fileName, "w") as text_file:
                for fileLine in fileLines:
                        #print "fileLine: ", fileLine
                        text_file.write(fileLine)


    import os
    mshName = fileName.replace(".geo",".msh")
    xmlName = fileName.replace(".geo",".xml")
    print "pre OS"
    line = "gmsh -3 %s"%fileName
    os.system(line)
    line="dolfin-convert %s %s"%(mshName,xmlName)
    os.system(line)
    print "post dolfin-convert"
    myPath = os.path.abspath(fileName)
    myPath = myPath.replace(".geo",".hdf5")
    print "myPath",myPath
    import fenics as fem
    mesh =  fem.Mesh(xmlName)
    facets =fem.MeshFunction("size_t", mesh, name+"_facet_region.xml")
    cells = fem.MeshFunction("size_t", mesh, name+"_physical_region.xml")
    hdf5=fem.HDF5File(mesh.mpi_comm(),myPath,'w')
    hdf5.write(mesh,'mesh')
    hdf5.write(facets,'facets')
    hdf5.write(cells,'cells')
    hdf5.close()
    print "post OS"
    print xmlName
    #print "The last physical surface number is: ", physNum
    return xmlName







#!/usr/bin/env python
import sys
#
# Message printed when program run without arguments
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose:
  Ryan's nanoporous w/ occlusion mesh generator

Usage:
"""
  msg+="to generate relevant xml provide a string with the number of nanopores desired per dimension (entering 9 gives 9x9 nanopores), radius of the nanopore, radius of the occlusion, dimension of reservoir (ex 15 makes a 15x15), length, and reservoir depth all lengths in nm then string for fileName:python %s 9,1e-9,9e-9,15e-9,10e-9,10e-9,fileName" %(scriptName)
  msg+="""


Notes:
  Pay attention to the units in the input arguments!
  Example: python nanoMesherNoAQ.py -makem "8,1e-9,15e-9,10e-9,10e-9,nameNoAQ"


"""
  return msg



if __name__ == "__main__":
  import sys
  import numpy as np
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  # Loops over each argument in the command line
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-makem"):
      mystr = sys.argv[i+1]
      myspl=mystr.split(',')
      print myspl[0]
      print myspl[1]
      Build(num=float(myspl[0]),radius = float(myspl[1]),boxSize=float(myspl[2]),length=float(myspl[3]),revH = float(myspl[4]),name=myspl[5])
      quit()
  raise RuntimeError("Arguments not understood")
