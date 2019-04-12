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

  this.frameRate = 0.00000000001;

  this.angmom = new Vector();
  this.positionCOM = new Vector();
  this.velocityCOM = new Vector();
  this.totalEnergy = undefined;
  this.initialEnergy = undefined;
  this.totalMass = 0.0;

  this.nOrbitPoints = 100; // Number of points for drawing orbits
  this.canvasID = 'myCanvas';
  this.pixscale = 100.0;
  this.timestepTolerance = 0.00001;
}

NBodySystem.prototype.calcTotalMass =
function calcTotalMass() {
  this.totalMass = 0.0;
  for (let i = 0; i< this.N; i++) {
    this.totalMass += this.bodies[i].mass;
  }
};

NBodySystem.prototype.calcCOM =
    function calcCOM() {
      this.positionCOM.setZero();
      this.velocityCOM.setZero();

      if (this.totalmass > 0.0) {
	    for (let ibody=0; ibody<this.N; ibody++) {
          this.positionCOM = this.positionCOM.add(this.bodies[ibody].position.scale(this.bodies[ibody].mass));
          this.velocityCOM = this.velocityCOM.add(this.bodies[ibody].velocity.scale(this.bodies[ibody].mass));
	    }

	    this.positionCOM = this.positionCOM.scale(1.0/this.totalMass);
	    this.velocityCOM = this.velocityCOM.scale(1.0/this.totalMass);
      }
    };

NBodySystem.prototype.calcTotalEnergy =
    function calcTotalEnergy() {
      let kineticEnergy = 0.0;
      let potentialEnergy = 0.0;

      for (let ibody=0; ibody< this.N; ibody++) {
	    const speed = this.bodies[ibody].velocity.getMag();
	    kineticEnergy += 0.5*this.bodies[ibody].mass*speed*speed;

	    for (let jbody=0; jbody< this.N; jbody++) {
          if (ibody === jbody) {
		    continue;
          }
          separation = this.bodies[ibody].position.subtract(this.bodies[jbody].position).getMag();
          potentialEnergy += separation >0.0 ? this.G*this.bodies[ibody].mass*this.bodies[jbody].mass/separation : 0.0;
	    }
      }

      if (this.totalEnergy == undefined)
      {
	    this.initialEnergy = kineticEnergy-potentialEnergy;
      }
      this.totalEnergy = kineticEnergy - potentialEnergy;

      if (this.initialEnergy !== undefined)
      {
	    this.dE = 100*(this.totalEnergy-this.initialEnergy)/this.initialEnergy;
      }
    };

NBodySystem.prototype.addBody = function addBody(body) {
  this.bodies.push(body);
  this.N++;
  this.calcTotalMass();
  this.calcCOM();
};


NBodySystem.prototype.changeFrame = function changeFrame(framePosition, frameVelocity) {
  for (let ibody=0; ibody< this.N; ibody++)
  {
    this.bodies[ibody].position = this.bodies[ibody].position.subtract(framePosition);
    this.bodies[ibody].velocity = this.bodies[ibody].velocity.subtract(frameVelocity);
  }
};

NBodySystem.prototype.changeToCOMFrame = function changeToCOMFrame() {
  this.changeFrame(this.positionCOM, this.velocityCOM);
};


NBodySystem.prototype.setupOrbits =
function setupOrbits() {
};


NBodySystem.prototype.calcTotalAngularMomentum =
    function calcTotalAngularMomentum() {
      this.angmom.setZero();
      for (let ibody=0; ibody < this.N; ibody++) {
        this.bodies[ibody].calcOrbitalAngularMomentum();
        this.angmom.add1(this.bodies[ibody].angmom);
      }
    };

NBodySystem.prototype.calcTimestep =
    function calcTimestep(dtmax = 1.0e30) {
      dtmin = 1.0e30;
      for (let ibody=0; ibody < this.N; ibody++) {
	    this.bodies[ibody].calcTimestep(this.timestepTolerance);
	    if (this.bodies[ibody].timestep < dtmin) {
          dtmin = this.bodies[ibody].timestep;
	    }
      }
      if (dtmin > dtmax) dtmin = dtmax;
      this.timestep = dtmin;
    };

NBodySystem.prototype.calcForces =
    function calcForces(bodyList) {
      const NBodies = bodyList.length;

      for (let i=0; i< NBodies; i++) {
	    bodyList[i].acceleration.setZero();
	    bodyList[i].jerk.setZero();
	    bodyList[i].snap.setZero();
	    bodyList[i].crackle.setZero();


	    bodyList[i].calcAccelJerk(this.G, bodyList, this.softeningLength);
	    bodyList[i].calcSnapCrackle(this.G, bodyList, this.softeningLength);
      }
    };


NBodySystem.prototype.evolveSystem =
    function evolveSystem() {
      /* Author: dh4gan
	 * This method updates the positions of the bodies via
	 * a 4th order Hermite integration (via a predictor-corrector algorithm)
	 */

      const dtmax = 0.5 * this.frameRate;
      /* i.  Calculate initial accelerations, jerks, snaps and crackles */
      this.calcForces(this.bodies);

      /* ii. Calculate initial global timestep */
      this.calcTimestep(dtmax);
      this.calcTotalEnergy();
      this.calcTotalAngularMomentum();


      const tend = this.time + this.frameRate;

      while (this.time < tend)
      {
	    const t2 = this.timestep * this.timestep;
	    const t3 = this.timestep * t2;

        /* iii. Perform the prediction step */

        // Arrays to store old co-ordinates
        const pos_old = [];
        const vel_old = [];
        const acc_old = [];
        const jerk_old = [];
        /* Calculate predicted positions and velocities */
        for (let i = 0; i < this.N; i++)
	    {
          // Store old velocity, acceleration and jerk
          pos_old.push(cloneVector(this.bodies[i].position));
          vel_old.push(cloneVector(this.bodies[i].velocity));
          acc_old.push(cloneVector(this.bodies[i].acceleration));
          jerk_old.push(cloneVector(this.bodies[i].jerk));

	    // Pull the body object's data //
	    const pos = this.bodies[i].position;
	    const vel = this.bodies[i].velocity;
	    const acc = this.bodies[i].acceleration;
	    const jerk = this.bodies[i].jerk;


	    // Calculate predicted position and velocity //
	    this.bodies[i].position = pos.add3(vel.scale(this.timestep), acc.scale(
		    0.5 * t2), jerk.scale(t3 / 6.0));

	    this.bodies[i].velocity = vel.add2(acc.scale(this.timestep), jerk.scale(
		    0.5 * t2));
	    }


	    /* iv. Use predicted positions and velocities to calculate
	     * predicted accelerations, jerks, snaps and crackles */

	    this.calcForces(this.bodies);

        /* v. Perform the correction step */
	    for (let i = 0; i < this.N; i++)
	    {
          const pos = this.bodies[i].position;
          const vel = this.bodies[i].velocity;
          const acc = this.bodies[i].acceleration;
          const jerk = this.bodies[i].jerk;

          let accterm = this.bodies[i].acceleration.add1(acc_old[i]).scale(0.5 * this.timestep);

          const jerkterm = this.bodies[i].jerk.relativeVector(jerk_old[i]).scale(t2/ 12.0);

          this.bodies[i].velocity = vel_old[i].add2(accterm, jerkterm);

          accterm = this.bodies[i].acceleration.relativeVector(acc_old[i]).scale(t2 / 12.0);
          const velterm = vel_old[i].add1(vel).scale(0.5 * this.timestep);
          this.bodies[i].position = pos_old[i].add2(velterm, accterm);
	    }

        this.calcForces(this.bodies);


	    this.time = this.time + this.timestep;
	    this.calcTimestep(dtmax);
	    this.calcTotalEnergy();
	    this.calcTotalAngularMomentum();
      }
    };

NBodySystem.prototype.drawSystem =
function drawSystem() {
  const canvas = document.getElementById(this.canvasID);
  const context = canvas.getContext('2d');
  context.clearRect(0, 0, canvas.width, canvas.height);

  this.evolveSystem();

  for (let ibody = 0; ibody < this.N; ibody++) {
    if (this.bodies[ibody].a > 1.0e-5)
    {
	  this.bodies[ibody].calcOrbitFromVector(this.G, this.totalMass);

      this.bodies[ibody].drawOrbit(this.G, this.totalMass,
          this.nOrbitPoints, this.canvasID, this.pixscale);
    }
    this.bodies[ibody].draw2D(this.canvasID, this.pixscale);

    const bodyElementId = 'Body'+(ibody+1);
    document.getElementById(bodyElementId).innerHTML = bodyElementId+this.bodies[ibody].printOrbit()+ '<br>'+this.bodies[ibody].printVectors();
  }

  document.getElementById('time').innerHTML = 'Time = '+this.time.toPrecision(4).toString() + '<br>Energy= '+this.totalEnergy.toPrecision(4).toString() + '<br> dE = '+this.dE.toPrecision(4).toString()+'%';
};

NBodySystem.prototype.Run =
    function Run(milliseconds=1000)
    {
      this.calcCOM();
      this.changeToCOMFrame();

      const interval = setInterval(this.drawSystem.bind(this), milliseconds);
    };

/**
 * Quick method to setup a test NBodySystem
 */
function testSystem() {
  const system = new NBodySystem();
  system.calcTotalMass();

  system.time = 0.0;

  system.addBody(new Body(1.0, 10.0, 'yellow',
      new Vector(0.0, 0.0, 0.0), new Vector(0.0, 0.0, 0.0)));

  system.addBody(createBodyFromOrbit(0.001, 10.0, 'green',
      system.G, system.totalMass+0.001, 1.0, 0.3, 0.0, 0.0, 1.7, 0.0));

  system.addBody(createBodyFromOrbit(0.001, 10.0, 'blue', system.G,
      system.totalMass+0.002, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

  system.frameRate = 0.01;

  return system;
};


const system = testSystem();

system.Run(30);

