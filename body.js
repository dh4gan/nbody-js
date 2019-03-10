/* eslint linebreak-style: ["error", "windows"] */
/* vector.js
 *  Written 05-March-2018 by dh4gan
 *
 *  Part of a trio that defines Javascript objects to be used
 *  for N Body simulation in the browser
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
 *  Attributes: (size, colour, pos, vel, acc, a, e, i,
 *  argPer, longAscend, trueAnomaly)
 *
 *  Methods: calcOrbitalAngularMomentum, calcEccentricity, calcTrueAnomaly,
    calcOrbitFromVector, calcVectorFromOrbit, calcPeriod
    draw
    ***************************************

 *  **********************************************
 *  nBodySystem:
 *  Attributes: (bodies, timestep, totalmass, G , softeningLength,
 *  rCOM, vCOM, aCOM)
 *  Methods: calcCOM, setupOrbits, calcTotalEnergy,
 *  calcTotalAngularMomentum, calcTimestep,

    drawSystem
    ***********************************************
 *
 *  RK4 or Hermite?
 */


/*
 ******************
 Body object constructor
 ******************
 size       - radius of body in pixels
 colour     - colour of body (string)
 pos        - position (vector)
 vel        - velocity (vector)

 a          - semimajor axis
 e          - eccentricity
 i          - inclination
 longascend - longitude of the ascending node
 argper     - argument of periapsis
 trueanom   - true anomaly
 period     - orbital period

 angmom     - orbital angular momentum (vector)
 eccvec     - eccentricity vector (Laplace-Runge-Lenz vector)

 */

/**
 * Body constructor
 * @param {double} mass - mass of the body
 * @param {double} size - radius of the body in pixels
 * @param {string} colour - colour of body
 * @param {vector} position - position
 * @param {vector} velocity - velocity
 *
*/
function Body(mass, size, colour, position, velocity) {
  this.mass = mass;
  this.size = size;
  this.colour = colour;
  this.position = position;
  this.velocity = velocity;

  this.a = 0.0; // semimajor axis
  this.e = 0.0; // eccentricity
  this.i = 0.0; // inclination
  this.longascend = 0.0; // longitude of the ascending node
  this.argper = 0.0; // argument of periapsis
  this.eccanom = 0.0; // eccentric anomaly
  this.trueanom = 0.0; // true anomaly
  this.meananom = 0.0; // mean anomaly
  this.period = 0.0; // orbital period

  this.angmom = new Vector(0.0, 0.0, 0.0);
  this.eccvec = new Vector(0.0, 0.0, 0.0);
}

// Body methods


Body.prototype.clone = function clone() {
  return new Body(this.mass, this.size, this.colour, this.pos, this.vel);
};

// Calculate orbital angular momentum of body
Body.prototype.calcOrbitalAngularMomentum =
function calcOrbitalAngularMomentum() {
  this.angmom = this.pos.cross(this.vel);
};


// Calculate orbital eccentricity vector
// (assumes orbital angular momentum up to date)

Body.prototype.calcEccentricity = function calcEccentricity(G, totalmass) {
  const magpos = this.pos.getMag();
  const magvel = this.vel.getMag();
  const vdotr = this.vel.dot(this.pos);
  const gravparam = G * totalmass;

  if (magpos === 0.0) {
    this.eccvec.setZero();
  } else {
    this.eccvec = this.pos.scale(magvel *
      magvel).subtract(this.vel.scale(vdotr / gravparam));

    console.log(vdotr, magvel);
    this.eccvec = this.eccvec.subtract(this.pos.scale(1.0 / magpos));
  }

  this.e = this.eccvec.getMag();
};


// Calculates orbital parameters from input position and velocity

Body.prototype.calcOrbitFromVector =
function calcOrbitFromVector(G, totalmass) {
  this.calcOrbitalAngularMomentum();
  this.calcEccentricity(G, totalmass);

  const angmag = this.angmom.getMag();
  const nplane = new Vector(0.0, 0.0, 0.0);
  let nscalar = 0.0;
  let magpos = 0.0;

  // Compute semimajor axis

  this.a = angmag * angmag / (G * totalmass * (1.0 - this.e * this.e));

  // Compute inclination

  if (angmag > 0.0) {
    this.i = Math.acos(this.angmom.z / angmag);
  } else {
    this.i = 0.0;
  }

  // Calculate Longitude of the Ascending Node

  if (this.i === 0.0) {
    this.longascend = 0.0;

    nplane.x = angmag;
    nplane.y = 0.0;
    nplane.z = 0.0;
    nscalar = nplane.getMag();
  } else {
    nplane.x = -this.angmom.y;
    nplane.y = this.angmom.x;
    nplane.z = 0.0;

    nscalar = nplane.getMag();
    this.longascend = Math.acos(nplane.x / nscalar);

    if (nplane.y < 0.0) {
      this.longascend = 2.0 * Math.pi - this.longascend;
    }
  }

  // Calculate true anomaly

  magpos = this.pos.getMag();

  // If orbit circular, no inclination, then use the position vector itself

  if (this.e === 0.0 && this.i === 0.0) {
    this.trueanom = Math.acos(this.pos.x / magpos);
    if (this.vel.x < 0.0) {
      this.trueanom = 2.0 * Math.pi - this.trueanom;
    }
  } else if (this.e === 0.0) {
    // If orbit circular and inclination non-zero,
    // then use the orbital plane vector

    let ndotR = nplane.dot(this.pos);
    ndotR /= (magpos * nscalar);

    const ndotV = nplane.dot(this.vel);

    this.trueanom = Math.acos(ndotR);

    if (ndotV > 0.0) {
      this.trueanom = 2.0 * Math.pi - this.trueanom;
    } else {
    // For non-circular orbits use the eccentricity vector
      edotR = this.eccvec.dot(this.pos);
      edotR /= (magpos * this.e);

      rdotV = this.vel.dot(this.pos);

      this.trueanom = Math.acos(edotR);

      if (rdotV < 0.0) {
        this.trueanom = 2.0 * Math.pi - this.trueanom;
      }
    }
  }

  // Finally, calculate the longitude of periapsis
  // First calculate the argument of periapsis

  if (this.e > 0.0) {
    edotn = this.eccvec.dot(nplane);
    edotn /= (nscalar * this.e);

    this.argper = Math.acos(edotn);
    if (this.eccvec.z < 0.0) {
      this.argper = 2.0 * Math.pi - this.argper;
    }

    this.longper = this.argper + this.longascend;
  } else {
    this.argper = 0.0;
    this.longper = 0.0;
  }
};

// Given orbital parameters, computes position and velocity of body

Body.prototype.calcVectorFromOrbit =
function calcVectorFromOrbit(G, totalmass) {
  console.log('Orbit', this.a, this.e, this.i,
      this.longascend, this.argper, this.trueanom);
  /* 1. calculate distance from CoM using semimajor axis,
eccentricity and true anomaly */

  const magpos = this.a * (1.0 - this.e * this.e) /
(1.0 + this.e * Math.cos(this.trueanom));

  /* 2. Calculate position vector in orbital plane */

  this.pos.x = magpos * Math.cos(this.trueanom);
  this.pos.y = magpos * Math.sin(this.trueanom);
  this.pos.z = 0.0;

  /* 3. Calculate velocity vector in orbital plane */
  const semiLatusRectum = Math.abs(this.a * (1.0 - this.e * this.e));
  const gravparam = G * totalmass;

  let magvel = 0.0;
  if (semiLatusRectum > 0.0) {
    magvel = Math.sqrt(gravparam / semiLatusRectum);
  } else {
    magvel = 0.0;
  }

  this.vel.x = -magvel * Math.sin(this.trueanom);
  this.vel.y = magvel * (Math.cos(this.trueanom) + this.e);
  this.vel.z = 0.0;

  /* 4. Begin rotations:
     * Firstly, Rotation around z axis by argument of Periapsis */

  if (this.argper != 0.0) {
    this.pos.rotateZ(-1 * this.argper);
    this.vel.rotateZ(-1 * this.argper);
  }

  /* Secondly, Rotate around x by inclination */

  if (this.i != 0.0) {
    this.pos.rotateX(-1 * this.i);
    this.vel.rotateX(-1 * this.i);
  }

  /* Lastly, Rotate around z by longitude of ascending node */

  if (this.longascend != 0.0) {
    this.pos.rotateZ(-1 * this.longascend);
    this.vel.rotateZ(-1 * this.longascend);
  }
};


// Returns the orbital period of the body
Body.prototype.calcPeriod = function calcPeriod(G, totalmass) {
  const period = Math.sqrt(4.0 * Math.pi * Math.pi *
    this.a * this.a * this.a / (G * totalmass));
  return period;
};

createBodyFromOrbit = function(mass, size, colour, G,
    totalmass, a, e, i, longascend, argper, trueanom) {
  const zeroVector1 = new Vector(0.0, 0.0, 0.0);
  const zeroVector2 = new Vector(0.0, 0.0, 0.0);
  newBody = new Body(mass, size, colour, zeroVector1, zeroVector2);
  newBody.a = a;
  newBody.e = e;
  newBody.i = i;
  newBody.longascend = longascend;
  newBody.argper = argper;
  newBody.trueanom = trueanom;

  newBody.vel.print();
  newBody.calcVectorFromOrbit(G, totalmass);


  console.log('New body parameters');
  newBody.pos.print();
  newBody.vel.print();
  return newBody;
};

// Body drawing methods
Body.prototype.draw2D = function(canvasID, pixscale) {
  const canvas = document.getElementById(canvasID);
  const context = canvas.getContext('2d');
  const centerX = canvas.width / 2;
  const centerY = canvas.height / 2;

  const offSetX = centerX + pixscale*this.pos.x;
  const offSetY = centerY + pixscale*this.pos.y;

  context.beginPath();
  context.arc(offSetX, offSetY, this.size,
      0, 2 * Math.PI, false);

  const grad = context.createRadialGradient(
      offSetX, offSetY, 0.0, offSetX, offSetY, this.size);

  grad.addColorStop(0, this.colour);
  grad.addColorStop(1, 'white');

  context.fillStyle = grad;
  context.fill();
};


// Draws an ellipse corresponding to a given orbit
Body.prototype.drawOrbit = function(G, totalmass, npoints, canvasID, pixscale) {
  const dummyBody = this.clone();
  dummyBody.calcOrbitFromVector(G, totalmass);

  const dphi = 2.0 * Math.PI / npoints;

  const canvas = document.getElementById(canvasID);
  const context = canvas.getContext('2d');
  const centerX = canvas.width / 2;
  const centerY = canvas.height / 2;
  let offSetX;
  let offSetY;

  context.moveTo(centerX + dummyBody.pos.x, centerY + dummyBody.pos.y);
  context.beginPath();
  context.strokeStyle = this.colour;
  context.linewidth = 0.5;
  context.setLineDash([2, 3]);

  for (i = 0; i < npoints; i++) {
    dummyBody.trueanom += dphi;
    if (dummyBody.trueanom > 2.0 * Math.Pi) {
      dummyBody.trueanom -= 2.0 * Math.Pi;
    }

    dummyBody.calcVectorFromOrbit(G, totalmass);

    offSetX = centerX + pixscale * (dummyBody.pos.x);
    offSetY = centerY + pixscale * (dummyBody.pos.y);
    context.lineTo(offSetX, offSetY);
    context.stroke();
    context.moveTo(offSetX, offSetY);
  }
  context.globalAlpha = 1.0;
};

/**
 * simple test function to check the body stuff works
 */
function testBody() {
  const pos1 = new Vector(0.707106, 0.707106, 0.0);
  const vel1 = new Vector(-0.707106, 0.707106, 0.0);
  const body1 = new Body(0.001, 10.0, 'blue', pos1, vel1);


  const npoints = 100;
  const totalmass = 1.0;
  const G = 1.0;

  const pixscale = 100.0;

  const body2 = createBodyFromOrbit(0.001, 10.0, 'green',
      G, totalmass, 1.0, 0.1, 0.0, 0.0, 0.0, 0.0);

  body1.calcOrbitFromVector(G, totalmass);

  body1.drawOrbit(G, totalmass, npoints, 'myCanvas', pixscale);
  body1.draw2D('myCanvas', pixscale);

  body2.drawOrbit(G, totalmass, npoints, 'myCanvas', pixscale);
  body2.draw2D('myCanvas', pixscale);
}


testBody();
