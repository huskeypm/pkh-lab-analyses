//This would be really nice if we could get it to work.

//Use solid geometry kernel
SetFactory("OpenCASCADE");

//Geometric parameters
radius = 3e0; // radius of pore
hsp = 3e0; // horizontal spacing
vsp = 3e0; // vertical spacing
RevW = 2*(hsp + radius); // width of reservoir
RevH = 2*(vsp + radius); // height of reservoir
RevD = 40e0; // depth of reservoir
length = 90e0; // length of pore
incRadius = 1.5e0; // radius of spherical inclusion
incDist = 1.5e0; //distance from pore surface to inclusion
vShift = radius - incRadius - incDist; //Distance from porce center to sphere center
incDepth =  2.0e0; // how far inside the pore the inclusion is

//Characteristic lengths for mesh refinement
lc1 = 1e0; //5e0;
lc2 = 1e0;
lc3 = 2e-1;

//Calculate key coordinates
//These coordinates are defined for positve values only, so they can be negated about any axis.
//Note that there may not be a geometric feature at all negative coordinates, though.
X0 = 0.0; //Center of model
X1 = incRadius; //Edge of inclusion
X2 = radius; //Edge of pore
X3 = radius + hsp;  //Model boundary
Y0 = 0.0; //Center of model
Y1 = vShift - incRadius; //Edge of inclusion farthest from pore wall
Y2 = vShift; //Center of inclusion
Y3 = vShift + incRadius; //Edge of inclusion closest to pore wall
Y4 = radius; //Edge of pore
Y5 = radius + vsp; //Model boundary
Z0 = 0.0; //Midline of the channel
Z1 = length/2.0 - incDepth - incRadius; //Edge of inclusion closest to midline
Z2 = length/2.0 - incDepth; //Center of inclusion
Z3 = length/2.0 - incDepth + incRadius; //Edge of inclusion farthest from midline
Z4 = length/2.0; //End of pore
Z5 = length/2.0 + RevD; //Model boundary

//Cylinder
cyl=newv;
Cylinder(cyl)={X0, Y0, -Z4, 0,0,length, radius};
//Cylinder(cyl)={X0, Y0, -Z5, 0,0,2*Z5, radius};

//Spherical inclusion
sph=newv;
Sphere(sph)={X0, Y2, Z2, incRadius};
//Get its surface points for mesh density specification
/*sphere_surf()=Unique(Abs(Boundary{Volume{sph}; }));
sphere_lines()=Unique(Abs(Boundary{Surface{sphere_surf()}; }));
sphere_points()=Unique(Abs(Boundary{Line{sphere_lines()}; }));*/

//Boolean subtract sphere from cylinder
cyl_sph=newv;
BooleanDifference(cyl_sph) = { Volume{cyl}; Delete; }{ Volume{sph}; Delete; };
/*cyl_surf() = Unique(Abs(Boundary{Volume{cyl_sph}; }));
cyl_lines() = Unique(Abs(Boundary{Surface{cyl_surf()}; }));
cyl_points() = Unique(Abs(Boundary{Line{cyl_lines()}; }));*/

//Top reservoir
toprev=newv;
Box(toprev)={-X3, -Y5, Z4, RevW, RevH, RevD};

//Bottom reservoir
botrev=newv;
Box(botrev)={-X3, -Y5, -Z5, RevW, RevH, RevD};

//Take union of everything
stage1=newv;
BooleanUnion(stage1)= {Volume{cyl_sph}; Delete; }{Volume{toprev}; Delete; };
allv=newv;
BooleanUnion(allv)= {Volume{stage1}; Delete; }{Volume{botrev}; Delete; };
//BooleanUnion(allv)= {Volume{cyl_sph}; Volume{toprev}; Volume{botrev}; Delete; }{};

//Get boundary surfaces
/*all_surf() = Unique(Abs(Boundary{ Volume{allv}; }));
all_lines() = Unique(Abs(Boundary{Surface{all_surf()}; }));
all_points() = Unique(Abs(Boundary{Line{all_lines()}; }));*/

//Characteristic length(s)
/*Characteristic Length{all_points()}= lc1;*/
//Characteristic Length{cyl_points()}= lc2;
//Characteristic Length{sphere_points()}= lc3;


//Mesh
//Mesh 3;