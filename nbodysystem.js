/* eslint linebreak-style: ["error", "windows"] */
/* vector.js
 *  Written 05-March-2018 by dh4gan
 *
 *  Part of a trio that defines Javascript objects to
 *  be used for N Body simulation in the browser
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
 *  Attributes: (size,colour,pos,vel,acc,
 * a, e, i, argPer, longAscend, trueAnomaly)
 *  Methods: calcOrbitalAngularMomentum, calcEccentricity, calcTrueAnomaly,
    calcOrbitFromVector, calcVectorFromOrbit, calcPeriod
    draw
    ***************************************

 *  **********************************************
 *  nBodySystem:
 *  Attributes: (bodies, timestep, totalmass, G ,
 *  softeningLength, rCOM, vCOM, aCOM)
 *  Methods: calcCOM, setupOrbits, calcTotalEnergy,
 *  calcTotalAngularMomentum, calcTimestep,

    drawSystem
    ***********************************************
 *
 *  RK4 or Hermite?
 */

/**
 * nBodySystem constructor
 * @param {double} timestep
 * @param {double} G
 * @param {double} softeningLength
 * @param {array} bodies
 *
 */
function NBodySystem(timestep = 0.001, G = 1.0,
    softeningLength = 1.0e-5, bodies = []) {
  this.timestep = timestep;
  this.G = G;
  this.softeningLength = softeningLength;
  this.bodies = bodies;
  this.N = bodies.length;
}

NBodySystem.prototype.calcTotalMass =
function calcTotalMass() {
  this.totalmass = 0.0;
  for (i = 0; i< this.N; i++) {
    this.totalmass += this.bodies[i].mass;
  }
};

NBodySystem.prototype.calcCOM =
function calcCOM() {
};


NBodySystem.prototype.addBody = function addBody(body) {
  this.bodies.push(body);
  this.N++;
  this.calcTotalMass();
  this.calcCOM();
};

NBodySystem.prototype.setupOrbits =
function setupOrbits() {
};

NBodySystem.prototype.calcTotalEnergy =
function calcTotalEnergy() {
};

NBodySystem.prototype.calcTotalAngularMomentum =
function calcTotalAngularMomentum() {
};

NBodySystem.prototype.calcTimestep =
function calcTimestep() {
};

NBodySystem.prototype.drawSystem =
function drawSystem(npoints, canvasID, pixscale) {
  pixscale = 100.0;
  for (ibody = 0; ibody < this.N; ibody++) {
    this.bodies[ibody].drawOrbit(this.G, this.totalmass,
        npoints, canvasID, pixscale);
    this.bodies[ibody].draw2D(canvasID, pixscale);
  }
};

/**
 * Quick method to test NBodySystem
 */
function testSystem() {
  const pixscale = 100.0;
  const system = new NBodySystem();
  system.calcTotalMass();

  system.addBody(new Body(1.0, 10.0, 'yellow',
      new Vector(0.0, 0.0, 0.0), new Vector(0.0, 0.0, 0.0)));

  system.addBody(createBodyFromOrbit(0.001, 10.0, 'green',
      system.G, system.totalmass+0.001, 1.0, 0.1, 0.0, 0.0, 0.0, 0.0));

  system.addBody(createBodyFromOrbit(0.01, 10.0, 'blue', system.G,
      system.totalmass+0.01, 3.0, 0.4, 0.0, 0.0, 1.3, 0.0, 4.0));
  system.drawSystem(100, 'myCanvas', pixscale);
};

testSystem();
