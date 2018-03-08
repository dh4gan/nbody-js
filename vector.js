/* vector.js
 *  Written 05-March-2018 by dh4gan
 *
 *  Part of a trio that defines Javascript objects to be used for N Body simulation
 *  in the browser
 *
 
 
 *  WARNING: Scripts must be included in following order:
 
 *  vector.js - Defines a 3D cartesian vector
 *  body.js - Defines a body with mass, position, velocity and orbital data
 *  nbodysystem.js - collects bodies into a system and integrates them
 
 
 
 
 ************************************
   Vector:
   
   Attributes:  (x,y,z, mag)
   Methods:     print, add, subtract, scale, unit, setZero,
                dot, cross,
                rotateX, rotateY, rotateZ
 *************************************
 
 
 * *************************************
 *  Body: 
 *  Attributes: (size,colour,pos,vel,acc, a, e, i, argPer, longAscend, trueAnomaly)
 *  Methods: calcOrbitalAngularMomentum, calcEccentricity, calcTrueAnomaly, 
    calcOrbitFromVector, calcVectorFromOrbit, calcPeriod
    draw
    ***************************************
 *  **********************************************
 *  nBodySystem:
 *  Attributes: (bodies, timestep, totalmass, G , softeningLength, rCOM, vCOM, aCOM)
 *  Methods: calcCOM, setupOrbits, calcTotalEnergy, calcTotalAngularMomentum, calcTimestep,

    drawSystem
    ***********************************************
 *
 *  RK4 or Hermite?
 */


/*
 **************************
 Vector object constructor
 **************************
 
 x,y,z  - cartesian coordinates
 mag    - magnitude of vector

 */

function Vector(xi,yi,zi){
    this.x = xi;
    this.y = yi;
    this.z = zi;
    this.mag = Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
    
};

/*

  Vector object methods

 */


// Print the vector to the console
Vector.prototype.print = function(){
    console.log("("+this.x+", "+this.y+", "+this.z+")");
};

// Compute vector magnitude
Vector.prototype.getMag = function(){
    return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
};


// Add vector this + v
Vector.prototype.add = function(v){
    return new Vector(this.x+v.x, this.y+v.y, this.z+v.z);
};

// Subtract v from vector this
Vector.prototype.subtract = function(v){
    return new Vector(this.x-v.x, this.y-v.y, this.z-v.z);
};


// Multiply vector by scalar
Vector.prototype.scale = function(fac){
    return new Vector(this.x*fac, this.y*fac, this.z*fac);
};

// Return the unit vector
Vector.prototype.unit = function(){   
    var unitVector = this.scale(1.0/this.mag);
    return unitVector;
};

// Sets vector to zero
Vector.prototype.setZero = function(){
    this.x = 0.0;
    this.y = 0.0;
    this.z = 0.0;
};

// Compute the dot product with vector v
Vector.prototype.dot = function(v){
    return this.x*v.x+ this.y*v.y + this.z*v.z;
};

// Compute the cross product with vector v
Vector.prototype.cross = function(v){

    xcross = this.y*v.z - this.z*v.y;
    ycross = -this.x*v.z + this.z*v.x;
    zcross = this.x*v.y - this.y*v.x;

    return new Vector(xcross, ycross,zcross);
};

// Rotate vector in x axis
Vector.prototype.rotateX = function(angle){

    var oldVector = new Vector(this.x,this.y,this.z);
    var sinang = Math.sin(angle)
    var cosang = Math.cos(angle);

    this.y = oldVector.y*cosang - oldVector.z*sinang;
    this.z = oldVector.y*sinang + oldVector.z*cosang;
};

// Rotate vector in y axis
Vector.prototype.rotateY = function(angle){

    var oldVector = new Vector(this.x,this.y,this.z);
    var sinang = Math.sin(angle)
    var cosang = Math.cos(angle);

    this.x = oldVector.x*cosang + oldVector.z*sinang;
    this.z = -oldVector.x*sinang + oldVector.z*cosang;
};

// Rotate vector in z axis
Vector.prototype.rotateZ = function(angle){

    var oldVector = new Vector(this.x,this.y,this.z);
    var sinang = Math.sin(angle)
    var cosang = Math.cos(angle);

    this.x = oldVector.x*cosang - oldVector.y*sinang;
    this.y = oldVector.x*sinang + oldVector.y*cosang;
};


function testVector(){

    var vector1 = new Vector(1,0,0);
    var vector2 = new Vector(0,1,0);

    var dot = vector1.dot(vector2);

    vector1.print();
    vector2.print();
    console.log("dot product is "+dot);

    vector3 = vector2.add(vector1);
    vector3.print();
    var vector4 = vector1.cross(vector2);

    vector4.print();

    vector4.rotateX(1.5707);
    vector4.print();

}



