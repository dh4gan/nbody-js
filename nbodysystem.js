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
 * @param {array} bodies
 * @param {double} timestep
 * @param {double} G
 * @param {double} softeningLength
 *
 */
function NBodySystem(bodies, timestep, G, softeningLength) {
  this.bodies = bodies;
  this.timestep = timestep;
  this.G = G;
  this.softeningLength = softeningLength;

  calcTotalMass();
}

NBodySystem.prototype.addBody = function addBody(body) {
  bodies.push(body);
  calcTotalMass();
  calcCOM();
};


NBodySystem.prototype.calcTotalMass =
function calcTotalMass() {
};

NBodySystem.prototype.calcCOM =
function calcCOM() {
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
function drawSystem() {
};

