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
 *  Attributes: (size, colour, position, velocity, acc, a, e, i,
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
 position        - position (vector)
 velocity        - velocity (vector)

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

  this.acceleration = new Vector();
  this.jerk = new Vector();
  this.snap = new Vector();
  this.crackle = new Vector();

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

  this.parentBody = null;
}

// Body methods

/**
* Converts a number to a string
* @param {number} x - number to convert
* @param {number} sigfig - number of significant figures
* @return {string} - string of x, to sigfig significant figures
*/
function precise(x, sigfig=3) {
  return Number.parseFloat(x).toPrecision(sigfig);
}

Body.prototype.print = function print() {
  console.log('Body: '+this.mass);
};

Body.prototype.printOrbit = function printOrbit() {
  const orbitString = '(a,e,i,argper,longascend,trueanom) = ('+
   precise(this.a) + '  '+ precise(this.e) + '  '+ precise(this.i)+ '  ' +
   precise(this.argper) + '  '+ precise(this.longascend) +
   '  ' + precise(this.trueanom)+')';
  return orbitString;
};

Body.prototype.printVectors = function printVectors() {
  const vectorString = '<br>Position: '+this.position.toString()+
    '<br>Velocity:'+this.velocity.toString()+
    '<br>Acceleration:'+this.acceleration.toString()+
    '<br>Jerk:'+this.jerk.toString();
  return vectorString;
};


// Calculate orbital angular momentum of body
Body.prototype.calcOrbitalAngularMomentum =
function calcOrbitalAngularMomentum() {
  this.angmom = this.position.cross(this.velocity);
};

Body.prototype.changeFrame =
function changeFrame(position, velocity) {
  this.position = this.position.subtract(position);
  this.velocity = this.velocity.subtract(velocity);
};

// Calculate orbital eccentricity vector
// (assumes orbital angular momentum up to date)

Body.prototype.calcEccentricity = function calcEccentricity(G, totalmass) {
  const magpos = this.position.getMag();
  const magvel = this.velocity.getMag();
  const vdotr = this.velocity.dot(this.position);
  const gravparam = G * totalmass;

  if (magpos < 1.0e-5) {
    this.eccvec.setZero();
  } else {
    this.eccvec = this.position.scale(magvel *
      magvel/gravparam - 1.0/magpos)
        .subtract(this.velocity.scale(vdotr / gravparam));
  }

  this.e = this.eccvec.getMag();
};


// Calculates orbital parameters from input position and velocity

Body.prototype.calcOrbitFromVector =
function calcOrbitFromVector(G, totalmass, parentBody=null) {
  let hostMass = totalmass;
  if (parentBody) {
    this.changeFrame(parentBody.position, parentBody.velocity);
    hostMass = parentBody.mass;
  }

  this.calcOrbitalAngularMomentum();
  this.calcEccentricity(G, hostMass);

  const angmag = this.angmom.getMag();
  const nplane = new Vector(0.0, 0.0, 0.0);
  let nscalar = 0.0;
  let magpos = 0.0;

  // Compute semimajor axis

  this.a = angmag * angmag / (G * hostMass * (1.0 - this.e * this.e));

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
      this.longascend = 2.0 * Math.PI - this.longascend;
    }
  }

  // Calculate true anomaly

  magpos = this.position.getMag();

  // If orbit circular, no inclination, then use the position vector itself

  if (this.e === 0.0 && this.i === 0.0) {
    if (magpos > 0.0) {
      this.trueanom = Math.acos(this.position.x / magpos);
    } else {
      this.trueanom = 0.0;
    }
    if (this.velocity.x < 0.0) {
      this.trueanom = 2.0 * Math.PI - this.trueanom;
    }
  } else if (this.e === 0.0) {
    // If orbit circular and inclination non-zero,
    // then use the orbital plane vector

    let ndotR = nplane.dot(this.position);
    ndotR /= (magpos * nscalar);

    const ndotV = nplane.dot(this.velocity);

    this.trueanom = Math.acos(ndotR);

    if (ndotV > 0.0) {
      this.trueanom = 2.0 * Math.PI - this.trueanom;
    }
  } else {
    // For non-circular orbits use the eccentricity vector
    edotR = this.eccvec.dot(this.position);
    edotR /= (magpos * this.e);

    rdotV = this.velocity.dot(this.position);

    edotR = edotR>1.0 ? 1.0 : edotR < -1.0 ? -1.0 : edotR;
    this.trueanom = Math.acos(edotR);


    if (rdotV < 0.0) {
      this.trueanom = 2.0 * Math.PI - this.trueanom;
    }
  }

  // Finally, calculate the longitude of periapsis
  // First calculate the argument of periapsis

  if (this.e > 0.0) {
    edotn = this.eccvec.dot(nplane);
    edotn /= (nscalar * this.e);

    this.argper = Math.acos(edotn);
    if (this.eccvec.z < 0.0) {
      this.argper = 2.0 * Math.PI - this.argper;
    }

    this.longper = this.argper + this.longascend;
  } else {
    this.argper = 0.0;
    this.longper = 0.0;
  }

  if (parentBody) {
    this.changeFrame(parentBody.position.scale(-1), parentBody.velocity.scale(-1));
  }
};

// Given orbital parameters, computes position and velocity of body
// We need a static version of this method for drawing orbits

calcPositionVelocityFromOrbit =
function calcPositionVelocityFromOrbit(G, totalmass,
    a, e, i, argper, longascend, trueanom) {
  /* 1. calculate distance from CoM using semimajor axis,
eccentricity and true anomaly */

  const magpos = a * (1.0 - e * e) /
(1.0 + e * Math.cos(trueanom));

  const position = new Vector();
  const velocity = new Vector();
  /* 2. Calculate position vector in orbital plane */

  position.x = magpos * Math.cos(trueanom);
  position.y = magpos * Math.sin(trueanom);
  position.z = 0.0;

  /* 3. Calculate velocity vector in orbital plane */
  const semiLatusRectum = Math.abs(a * (1.0 - e * e));
  const gravparam = G * totalmass;

  let magvel = 0.0;
  if (semiLatusRectum > 0.0) {
    magvel = Math.sqrt(gravparam / semiLatusRectum);
  } else {
    magvel = 0.0;
  }

  velocity.x = -magvel * Math.sin(trueanom);
  velocity.y = magvel * (Math.cos(trueanom) + e);
  velocity.z = 0.0;

  /* 4. Begin rotations:
     * Firstly, Rotation around z axis by argument of Periapsis */

  if (argper != 0.0) {
    position.rotateZ(-1 * argper);
    velocity.rotateZ(-1 * argper);
  }

  /* Secondly, Rotate around x by inclination */

  if (i != 0.0) {
    position.rotateX(-1 * i);
    velocity.rotateX(-1 * i);
  }

  /* Lastly, Rotate around z by longitude of ascending node */

  if (longascend != 0.0) {
    position.rotateZ(-1 * longascend);
    velocity.rotateZ(-1 * longascend);
  }

  return [position, velocity];
};


Body.prototype.calcVectorFromOrbit =
function calcVectorFromOrbit(G, totalmass) {
  const self=this;

  let hostMass = totalmass;
  if (this.parentBody) {
    hostMass = this.parentBody.mass;
  }

  const positionVelocity = calcPositionVelocityFromOrbit(G, hostMass,
      self.a, self.e, self.i,
      self.argper, self.longascend, self.trueanom);
  this.position = positionVelocity[0];
  this.velocity = positionVelocity[1];

  if (this.parentBody) {
    this.changeFrame(this.parentBody.position.scale(-1), this.parentBody.velocity.scale(-1));
  }
};


// Returns the orbital period of the body
Body.prototype.calcPeriod = function calcPeriod(G, totalmass) {
  const period = Math.sqrt(4.0 * Math.pi * Math.pi *
    this.a * this.a * this.a / (G * totalmass));
  return period;
};

createBodyFromOrbit = function createBodyFromOrbit(mass, size, colour, G,
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

  newBody.velocity.print();
  newBody.calcVectorFromOrbit(G, totalmass);


  console.log('New body parameters');
  newBody.position.print();
  newBody.velocity.print();
  return newBody;
};

// Body drawing methods
Body.prototype.draw2D = function draw2D(canvasID, pixscale) {
  const canvas = document.getElementById(canvasID);
  const context = canvas.getContext('2d');

  const centerX = canvas.width / 2;
  const centerY = canvas.height / 2;

  const offSetX = centerX + pixscale*this.position.x;
  const offSetY = centerY + pixscale*this.position.y;

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
Body.prototype.drawOrbit = function drawOrbit(G, totalmass,
    npoints, canvasID, pixscale) {
  self=this;

  const canvas = document.getElementById(canvasID);
  const context = canvas.getContext('2d');

  if (self.a >0.0 || self.position.mag > 1.0e-5) {
    const minoraxis = this.a*Math.sqrt((1.0-this.e*this.e));
    const focusDistance = Math.sqrt(this.a*this.a - minoraxis*minoraxis);

    let focusX = canvas.width / 2 -
    focusDistance*pixscale*Math.cos(2.0*Math.PI - this.argper);

    let focusY = canvas.height / 2 -
    focusDistance*pixscale*Math.sin(2.0*Math.PI - this.argper);

    if (this.parentBody) {
      focusX +=this.parentBody.position.x * pixscale;
      focusY +=this.parentBody.position.y * pixscale;
    }

    context.beginPath();
    context.strokeStyle = this.colour;
    context.linewidth = 0.01;

    context.ellipse(focusX, focusY,
        this.a*pixscale, minoraxis*pixscale, this.argper, 0.0, 2.0 * Math.PI);
    context.stroke();
  } else {
    // If orbit open, use the multiple point draw
    const dphi = 2.0 * Math.PI / npoints;
    let centerX = canvas.width / 2;
    let centerY = canvas.height / 2;

    if (this.parentBody) {
      centerX +=this.parentBody.position.x * pixscale;
      centerY +=this.parentBody.position.y * pixscale;
    }

    let offSetX = centerX + self.position.x * pixscale;
    let offSetY = centerY + self.position.y * pixscale;

    context.moveTo(offSetX, offSetY);
    context.beginPath();
    context.strokeStyle = self.colour;
    context.linewidth = 0.5;
    context.setLineDash([2, 3]);

    let nu = 0.0;
    for (i = 0; i < npoints; i++) {
      nu +=dphi;

      const positionVelocity = calcPositionVelocityFromOrbit(G, totalmass,
          self.a, self.e, self.i,
          self.argper, self.longascend, nu);

      offSetX = centerX + pixscale * (positionVelocity[0].x);
      offSetY = centerY + pixscale * (positionVelocity[0].y);
      context.lineTo(offSetX, offSetY);
      context.stroke();
      context.moveTo(offSetX, offSetY);
    }
  }

  context.globalAlpha=1.0;
};


Body.prototype.calcAccelJerk =
function calcAccelJerk(G, bodyarray, softeningLength) {
  /* Author: dh4gan
     * Calculate the gravitational acceleration on a Body
     * given an array of to bodies
     * Method :
     * Loop over each body in the system and calculate the vector r
     * from the body in question, then work out the force vector
     * and therefore the acceleration vector
     *
     * a_i = -G * M_i
     *    ------------- * r^hat_i
     *         |r_i|**2
     *
     * j_i = -G * M_i
     *  ------------- * v^hat_i - 3 alpha_i a_i
     *        |r_i|**2
     *
     * alpha_i = (r_i.v_i)/|r_i|**2
     *
     *
     * r^hat is the unit vector in the direction of the vector connecting
     * M1 is the mass of the body
     * M2 is the mass of the other body
     * |r| is the distance between the two bodies
     *
     *
     * Argument :
     *   bodyarray : a vector containing all the Body objects
     *    in the system that the Body in question belongs to
     *
     */
  let b;
  const N = bodyarray.length;

  let alpha;
  let factor;
  let rmag;
  let r2;
  let r3;
  let r21;

  relativeVelocity = new Vector();
  relativePosition = new Vector();

  accelterm = new Vector();
  jerkterm = new Vector();
  jerkterm1 = new Vector();
  jerkterm2 = new Vector();

  for (b = 0; b < N; b++) {
    // Get relative position and velocity

    relativePosition = this.position.relativeVector(bodyarray[b].position);
    relativeVelocity = this.velocity.relativeVector(bodyarray[b].velocity);

    rmag = relativePosition.getMag();

    // Since the body in question is inside the body array
    // I need to make sure the body I am looping through
    // isn't itself so skip the loop if distance effectively 0

    if (rmag < 1.0e-2 * softeningLength) {
      continue;
    }

    // Add gravitational softening
    r2 = rmag * rmag + softeningLength * softeningLength;
    r21 = 1.0 / r2;
    rmag = Math.sqrt(r2);
    r3 = rmag * rmag * rmag;

    // Define this factor, as it's useful
    factor = -G * bodyarray[b].mass / r3;

    // Calculate acceleration term
    accelterm = relativePosition.scale(factor);

    // acceleration = acceleration - accelterm
    this.acceleration = accelterm.relativeVector(this.acceleration);

    // now jerk - calculate alpha term
    alpha = relativeVelocity.dot(relativePosition);

    alpha = alpha * r21;
    jerkterm1 = relativeVelocity.scale(factor);
    jerkterm2 = accelterm.scale(-3.0 * alpha);

    // jerkterm = jerkterm2 - jerkterm1
    jerkterm = jerkterm1.add1(jerkterm2);

    // jerk = jerk - jerkterm
    this.jerk = jerkterm.relativeVector(this.jerk);
  } // End of loop
  // End of method
};

Body.prototype.calcSnapCrackle = function calcSnapCrackle(G, bodyarray,
    softeningLength) {
  /* Author: dh4gan 11/3/13
     * Calculate the snap and crackle on a Body given an array of bodies
     * WARNING: This method does not work unless calcAccelJerk is called first
     * This method will calculate a lot of the same terms as calcAccelJerk, but
     * this is the only way these calculations work on a Body by Body basis
     */

  let b;
  const N = bodyarray.length;

  let alpha;
  let beta;
  let gamma;
  let factor;
  let rmag;
  let r2;
  let r3;
  let r21;
  let vmag;
  let v2;

  zerovector = new Vector();
  relativeVelocity = new Vector();
  relativePosition = new Vector();
  relativeAcceleration = new Vector();
  relativeJerk = new Vector();

  accelterm = new Vector();

  jerkterm = new Vector();
  jerkterm1 = new Vector();
  jerkterm2 = new Vector();

  snapterm = new Vector();
  snapterm1 = new Vector();
  snapterm2 = new Vector();
  snapterm3 = new Vector();

  crackleterm = new Vector();
  crackleterm1 = new Vector();
  crackleterm2 = new Vector();
  crackleterm3 = new Vector();
  crackleterm4 = new Vector();

  for (b = 0; b < N; b++) {
    // Get relative position, velocity, acceleration and jerk

    relativePosition = this.position.relativeVector(bodyarray[b].position);
    relativeVelocity = this.velocity.relativeVector(bodyarray[b].velocity);
    relativeAcceleration = this.acceleration.relativeVector(
        bodyarray[b].acceleration);
    relativeJerk = this.jerk.relativeVector(bodyarray[b].jerk);

    rmag = relativePosition.getMag();
    vmag = relativeVelocity.getMag();

    v2 = vmag * vmag;

    // Since the body in question is inside the body array
    // I need to make sure the body I am looping through
    // isn't itself so skip the loop if distance ==0

    if (rmag < 1.0e-2 * softeningLength) {
      continue;
    }

    // Add gravitational softening
    r2 = rmag * rmag + softeningLength * softeningLength;
    r21 = 1.0 / r2;
    rmag = Math.sqrt(r2);
    r3 = rmag * rmag * rmag;

    // Define this factor, as it's useful
    factor = G * bodyarray[b].mass / r3;

    // Calculate acceleration term
    accelterm = relativePosition.scale(factor);

    // now jerk - calculate alpha term
    alpha = relativeVelocity.dot(relativePosition) * r21;

    jerkterm1 = relativeVelocity.scale(factor);
    jerkterm2 = accelterm.scale(3 * alpha);

    // jerkterm = jerkterm2 - jerkterm1
    jerkterm = jerkterm2.relativeVector(jerkterm1);

    // calculate snap terms
    beta = (v2 + relativePosition.dot(relativeAcceleration))
    * r21 + alpha * alpha;

    snapterm1 = relativeAcceleration.scale(factor);
    snapterm2 = jerkterm.scale(-6 * alpha);
    snapterm3 = accelterm.scale(-3 * beta);

    snapterm = snapterm1.add2(snapterm2, snapterm3);

    this.snap = snapterm.relativeVector(this.snap);

    // Finally crackle terms

    gamma = (3.0 * relativeVelocity.dot(relativeAcceleration)
    + relativePosition.dot(relativeJerk)) * r21;
    gamma = gamma + alpha * (3.0 * beta - 4.0 * alpha * alpha);

    crackleterm1 = relativeJerk.scale(factor);
    crackleterm2 = snapterm.scale(-9.0 * alpha);
    crackleterm3 = jerkterm.scale(-9.0 * beta);
    crackleterm4 = accelterm.scale(-3.0 * gamma);

    crackleterm = crackleterm1.add3(crackleterm2, crackleterm3,
        crackleterm4);

    this.crackle = crackleterm.relativeVector(this.crackle);
  } // End of loop
  // End of method
};

Body.prototype.calcTimestep = function calcTimestep(greekEta) {
  /* Author: dh4gan
    * Calculate the preferred timestep for the Body
    * given its acceleration, jerk, snap and crackle
    */

  const tolerance = 1e-20;

  const normJ = this.jerk.getMag();
  const normA = this.acceleration.getMag();
  const normC = this.crackle.getMag();
  const normS = this.snap.getMag();

  // If numerator zero, give a warning

  if (normA * normS + normJ * normJ < tolerance) {
    console.log('warning in calcTimestep: numerator zero for '+this.name);
    console.log(normJ, normA, normS, tolerance);
    this.position.print();
    console.log(this.position.getMag() + '  ' + this.velocity.getMag());
    this.timestep = 0.0;
  } else if (normC * normJ + normS * normS < tolerance) {
    // If denominator zero, give a warning - set timestep very large
    console.log('warning in calcTimestep: denominator zero for '+this.name);
    this.timestep = 1.0e30;
  } else {
    // Otherwise calculate timestep as normal
    this.timestep = Math.sqrt(greekEta * (normA * normS + normJ * normJ) /
    (normC * normJ + normS * normS));
  }
};
