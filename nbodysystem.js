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

    this.frameRate = 0.01

    this.angmom = new Vector();
    this.positionCOM = new Vector();
    this.velocityCOM = new Vector();
    this.totalEnergy = 0.0;
    this.totalMass = 0.0;

    this.nOrbitPoints = 100 // Number of points for drawing orbits
    this.canvasID = "myCanvas"
    this.pixscale = 100.0;
    this.timestepTolerance = 0.001;
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

	if(this.totalmass > 0.0) {
	
	    for(let ibody=0; ibody<this.N; ibody++) {
		this.positionCOM = this.positionCOM.add(this.bodies[ibody].position.scale(this.bodies[ibody].mass));
		this.velocityCOM = this.velocityCOM.add(this.bodies[ibody].velocity.scale(this.bodies[ibody].mass));
		
	    }

	    this.positionCOM = this.positionCOM.scale(1.0/this.totalMass)
	    this.velocityCOM = this.velocityCOM.scale(1.0/this.totalMass);
	}
    };

NBodySystem.prototype.calcTotalEnergy =
    function calcTotalEnergy() {

	let kineticEnergy = 0.0;
	let potentialEnergy = 0.0;
	
	for (let ibody=0; ibody< this.N; ibody++) {
	    let speed = this.bodies[ibody].velocity.getMag();
	    kineticEnergy += 0.5*this.bodies[ibody].mass*speed*speed;

	    for (let jbody=0; jbody< this.N; jbody++) {
		if(ibody === jbody){
		    continue;
		}
		separation = this.bodies[ibody].position.subtract(this.bodies[jbody].position).getMag();
		potentialEnergy += separation >0.0 ? this.G*this.bodies[ibody].mass*this.bodies[jbody].mass/separation : 0.0;
		
	    }
	}

	this.totalEnergy = kineticEnergy - potentialEnergy;
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



NBodySystem.prototype.calcTotalAngularMomentum =
    function calcTotalAngularMomentum() {

	this.angmom.setZero();
	for(let ibody=0; ibody < this.N; ibody++) {
	    this.bodies[ibody].calcOrbitalAngularMomentum();
	    this.angmom.add1(this.bodies[ibody].angmom);
	}
};

NBodySystem.prototype.calcTimestep =
    function calcTimestep(dtmax = 1.0e30) {

	dtmin = 1.0e30;
	for (let ibody=0; ibody < this.N; ibody++) {
	    this.bodies[ibody].calcTimestep(this.timestepTolerance);
	    if(this.bodies[ibody].timestep < dtmin) {
		dtmin = this.bodies[ibody].timestep;
	    }
	}
	if(dtmin > dtmax) dtmin = dtmax;
	this.timestep = dtmin;
};

NBodySystem.prototype.calcForces =
    function calcForces(bodyList) {

	let NBodies = bodyList.length;

	for (let i=0; i< NBodies; i++) {
	    
	    bodyList[i].acceleration.setZero();
	    bodyList[i].jerk.setZero();
	    bodyList[i].snap.setZero();
	    bodyList[i].crackle.setZero();
        

	    
	    bodyList[i].calcAccelJerk(this.G, bodyList, this.softeningLength);
	    bodyList[i].calcSnapCrackle(this.G, bodyList, this.softeningLength);
        
	}
    }


NBodySystem.prototype.evolveSystem =
    function evolveSystem() {
	/* Author: dh4gan
	 * This method updates the positions of the bodies via 
	 * a 4th order Hermite integration (via a predictor-corrector algorithm)
	 */

	let dtmax = 0.5 * this.frameRate;
	/* i.  Calculate initial accelerations, jerks, snaps and crackles */
	this.calcForces(this.bodies);

	/* ii. Calculate initial global timestep */
	this.calcTimestep(dtmax);
	this.calcTotalEnergy();
	this.calcTotalAngularMomentum();
    
	/* iii. Clone the current body array */

	let predicted = [];

	for (let i=0; i<this.N; i++) {
	    predicted.push(this.bodies[i].clone());
	}
        
        let counter = 0;

	let tend = this.time + this.frameRate;
	
	while (this.time < tend)
	{
	    let t2 = this.timestep * this.timestep;
	    let t3 = this.timestep * t2;
	    
	/* Calculate predicted positions and velocities */
	for (let i = 0; i < this.N; i++)
	    {
        
	    // Pull the body object's data //
	    let pos = this.bodies[i].position;
	    let vel = this.bodies[i].velocity;
	    let acc = this.bodies[i].acceleration;
	    let jerk = this.bodies[i].jerk;


	    // 1. Calculate predicted position and velocity //
	    predicted[i].position = pos.add3(vel.scale(this.timestep), acc.scale(
		    0.5 * t2), jerk.scale(t3 / 6.0));

	    predicted[i].velocity = vel.add2(acc.scale(this.timestep), jerk.scale(
		    0.5 * t2));

	    }
	

	    /* 2. Use predicted positions and velocities to calculate
	     * predicted accelerations, jerks, snaps and crackles */

	    this.calcForces(predicted);


	    for (let i = 0; i < this.N; i++)
	    {

		let pos_p = predicted[i].position;
		let vel_p = predicted[i].velocity;
		let acc_p = predicted[i].acceleration;
		let jerk_p = predicted[i].jerk;

		let pos = this.bodies[i].position;
		let vel = this.bodies[i].velocity;
		let acc = this.bodies[i].acceleration;
		let jerk = this.bodies[i].jerk;

		let accterm = acc_p.add1(acc).scale(0.5 * this.timestep);

		let jerkterm = jerk_p.relativeVector(jerk).scale(t2/ 12.0);

		this.bodies[i].velocity = vel.add2(accterm, jerkterm);

		accterm = acc_p.relativeVector(acc).scale(t2 / 12.0);
		let velterm = this.bodies[i].velocity.add1(vel).scale(0.5 * this.timestep);
		this.bodies[i].position = pos.add2(velterm, accterm);
		
	    }
	
        this.calcForces(this.bodies);
        

	    this.time = this.time + this.timestep;
	    this.calcTimestep(dtmax);
	    this.calcTotalEnergy();
	    this.calcTotalAngularMomentum();
        
        counter++;
       
	}
    };

NBodySystem.prototype.drawSystem =
function drawSystem() {
    
  
    var canvas = document.getElementById(this.canvasID);
    var context = canvas.getContext('2d');
    context.clearRect(0,0,canvas.width, canvas.height);

     this.evolveSystem();
    
  for (let ibody = 0; ibody < this.N; ibody++) {
      
      if(this.bodies[ibody].a > 1.0e-5)
      {
	  this.bodies[ibody].calcOrbitFromVector(this.G, this.totalMass);
	  
    this.bodies[ibody].drawOrbit(this.G, this.totalMass,
        this.nOrbitPoints, this.canvasID, this.pixscale);
      }
    this.bodies[ibody].draw2D(this.canvasID, this.pixscale);

      let bodyElementId = "Body"+(ibody+1);
      document.getElementById(bodyElementId).innerHTML = bodyElementId+this.bodies[ibody].printOrbit()+ "<br>"+this.bodies[ibody].printVectors();

      
  }

    document.getElementById("time").innerHTML = "Time = "+this.time.toPrecision(4).toString();
};

NBodySystem.prototype.Run =
    function Run(milliseconds=1000)
{

    var _this = this;
 
    var interval = setInterval(this.drawSystem.bind(this),milliseconds);

}

/**
 * Quick method to setup a test NBodySystem
 */
function testSystem() {
  const pixscale = 100.0;
  const system = new NBodySystem();
  system.calcTotalMass();
    
    system.time = 0.0;
 
  system.addBody(new Body(1.0, 10.0, 'yellow',
      new Vector(0.0, 0.0, 0.0), new Vector(0.0, 0.0, 0.0)));

  //system.addBody(createBodyFromOrbit(0.001, 10.0, 'green',
  //    system.G, system.totalMass+0.001, 1.0, 0.1, 0.0, 0.0, 0.0, 0.0));

  system.addBody(createBodyFromOrbit(0.01, 10.0, 'blue', system.G,
      system.totalMass+0.01, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    
    system.frameRate = 0.01

    return system;
    
};




var system = testSystem();

system.Run(10);

