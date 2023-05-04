"use strict";

// Particle constructor
function Particle(x, y, z, mass) {
  this.position = new THREE.Vector3(); // position
  this.previous = new THREE.Vector3(); // previous
  this.original = new THREE.Vector3(); // original
  initParameterizedPosition(x, y, this.position);
  initParameterizedPosition(x, y, this.previous);
  initParameterizedPosition(x, y, this.original);

  this.netForce = new THREE.Vector3(); // net force acting on particle
  this.mass = mass; // mass of the particle
  this.correction = new THREE.Vector3(); // offset to apply to enforce constraints
}

// Snap a particle back to its original position
Particle.prototype.lockToOriginal = function() {
  this.position.copy(this.original);
  this.previous.copy(this.original);
};

// Snap a particle back to its previous position
Particle.prototype.lock = function() {
  this.position.copy(this.previous);
  this.previous.copy(this.previous);
};

// Add the given force to a particle's total netForce.
// Params:
// * force: THREE.Vector3 - the force to add
Particle.prototype.addForce = function(force) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 1 lines of code.
  this.netForce.add(force);
  // ----------- STUDENT CODE END ------------
};

// Perform Verlet integration on this particle with the provided
// timestep deltaT.
// Params:
// * deltaT: Number - the length of time dt over which to integrate
Particle.prototype.integrate = function(deltaT) {
  const DAMPING = SceneParams.DAMPING;

  // ----------- STUDENT CODE BEGIN ------------
  // You need to:
  // (1) Save the old (i.e. current) position into this.previous.
  // (2) Compute the new position of this particle using Verlet integration,
  //     and store it into this.position.
  // (3) Reset the net force acting on the particle (i.e. make it (0, 0, 0) again).
  // ----------- Our reference solution uses 13 lines of code.
  let oldPos = this.position.clone();
  this.position = this.position.multiplyScalar(2 - DAMPING).sub(
    this.previous.multiplyScalar(1 - DAMPING)).add(
      this.netForce.multiplyScalar(deltaT*deltaT / this.mass));
  this.previous = oldPos;
  this.netForce = new THREE.Vector3();
  // ----------- STUDENT CODE END ------------
};

// Handle collisions between this Particle and the provided floor.
// Note: the fields of floor are documented for completeness, but you
//       *WILL NOT* need to use all of them.
// Params:
// * floor: An object representing the floor of the scene, with properties:
//    - mesh: THREE.Mesh - the physical representation in the scene
//    - geometry: THREE.PlaneBufferGeometry - the abstract geometric representation
//    - material: THREE.MeshPhongMaterial - material information for lighting
Particle.prototype.handleFloorCollision = function(floor) {
  let floorMesh = floor.mesh;
  let floorPosition = floorMesh.position;
  const EPS = 3;
  // ----------- STUDENT CODE BEGIN ------------
  // Handle collision of this particle with the floor.
  // ----------- Our reference solution uses 4 lines of code.
  if (this.position.y <= floorPosition.y + EPS){
    this.position.y = floorPosition.y + EPS;
  }
  // ----------- STUDENT CODE END ------------
};

// Handle collisions between this Particle and the provided sphere.
// Note: the fields of sphere are documented for completeness, but you
//       *WILL NOT* need to use all of them.
// Params:
// * sphere: An object representing a sphere in the scene, with properties:
//    - mesh: THREE.Mesh - the physical representation in the scene
//    - geometry: THREE.SphereGeometry - the abstract geometric representation
//    - material: THREE.MeshPhongMaterial - material information for lighting
//    - radius: number - the radius of the sphere
//    - position: THREE.Vector3 - the sphere's position in this frame
//    - prevPosition: THREE.Vector3 - the sphere's position in the previous frame
Particle.prototype.handleSphereCollision = function(sphere) {
  if (sphere.mesh.visible) {
    const friction = SceneParams.friction;
    let spherePosition = sphere.position.clone();
    let prevSpherePosition = sphere.prevPosition.clone();
    let EPS = 5; // empirically determined
    // ----------- STUDENT CODE BEGIN ------------
    // Handle collision of this particle with the sphere.
    // As with the floor, use EPS to prevent clipping.
    let posFriction = new THREE.Vector3();
    let posNoFriction = new THREE.Vector3();

    let dist = spherePosition.distanceTo(this.position);

    if (dist <= sphere.radius + EPS){
      let scalar = (sphere.radius + EPS) / dist;
      posNoFriction = this.position.sub(spherePosition).multiplyScalar(scalar);
      posNoFriction.add(spherePosition);

      this.position = posNoFriction;

      if (sphere.position.distanceTo(this.previous) > sphere.radius + EPS){
        posFriction = this.previous.clone().add(sphere.position).sub(prevSpherePosition);
        this.position.multiplyScalar(1 - friction).add(posFriction.multiplyScalar(friction));
        
      }
    }
    // ----------- Our reference solution uses 28 lines of code.
    // ----------- STUDENT CODE END ------------
  }
};

// Handle collisions between this Particle and the provided axis-aligned box.
// Note: the fields of box are documented for completeness, but you
//       *WILL NOT* need to use all of them.
// Params:
// * box: An object representing an axis-aligned box in the scene, with properties:
//    - mesh: THREE.Mesh - the physical representation in the scene
//    - geometry: THREE.BoxGeometry - the abstract geometric representation
//    - material: THREE.MeshPhongMaterial - material information for lighting
//    - boundingBox: THREE.Box3 - the bounding box of the box in the scene
Particle.prototype.handleBoxCollision = function(box) {
  if (box.mesh.visible) {
    const friction = SceneParams.friction;
    let boundingBox = box.boundingBox.clone();
    const EPS = 10; // empirically determined
    // ----------- STUDENT CODE BEGIN ------------
    // Handle collision of this particle with the axis-aligned box.
    // As before, use EPS to prevent clipping
    let posFriction = new THREE.Vector3();
    let posNoFriction = new THREE.Vector3();

    boundingBox.expandByScalar(EPS);

    let xMaxDist = boundingBox.max.x - this.position.x; 
    if(xMaxDist < 0) {
      return;
    }
    
    let xMinDist = this.position.x - boundingBox.min.x; 
    if(xMinDist < 0) {
      return;
    }

    let yMaxDist = boundingBox.max.y - this.position.y; 
    if(yMaxDist < 0) {
      return;
    }

    let yMinDist= this.position.y - boundingBox.min.y; 
    if(yMinDist < 0){
      return;
    }

    let zMaxDist = boundingBox.max.z - this.position.z; 
    if(zMaxDist < 0) {
      return;
    }

    let zMinDist = this.position.z - boundingBox.min.z; 
    if(zMinDist < 0) {
      return;
    }

    let minDist = Math.min(xMaxDist, xMinDist, yMaxDist, yMinDist, zMaxDist, zMinDist);

    if(xMaxDist == minDist){
      this.position.x = boundingBox.max.x;
    }
    else if(xMinDist == minDist){
      this.position.x = boundingBox.min.x;
    }
    else if(yMaxDist == minDist){
      this.position.y = boundingBox.max.y;
    }
    else if(yMinDist == minDist){
      this.position.y = boundingBox.min.y;
    }
    else if(zMaxDist == minDist){
      this.position.z = boundingBox.max.z;
    }
    else if(zMinDist == minDist){
      this.position.z =  boundingBox.min.z;
    }

    posNoFriction = this.position.clone();

    if (!boundingBox.containsPoint(this.previous)){
      posFriction = this.previous.clone();
      this.position = posNoFriction.multiplyScalar(1-friction).add(posFriction.multiplyScalar(friction));
    }
    
    // ----------- Our reference solution uses 66 lines of code.
    // ----------- STUDENT CODE END ------------
  }
};

// ------------------------ Don't worry about this ---------------------------
// Apply the cached correction vector to this particle's position, and
// then zero out the correction vector.
// Particle.prototype.applyCorrection = function() {
//   this.position.add(this.correction);
//   this.correction.set(0,0,0);
// }
