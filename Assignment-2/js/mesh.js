// In this file you will implement traversal and analysis for your assignment.
// Make sure to familiarize yourself with the utility functions in meshUtils.js
// they might be useful for the second part of your assignment!

////////////////////////////////////////////////////////////////////////////////
// Traversal
////////////////////////////////////////////////////////////////////////////////

// Return all vertices on face f
Mesh.prototype.verticesOnFace = function(f) {
  const vertices = [];
  let he = f.halfedge;
  const first = he;
  while (true) {
    vertices.push(he.vertex);
    he = he.next;
    if (he === first) {
      break;
    }
  }
  return vertices;
};

// Return all halfedges on face f
Mesh.prototype.edgesOnFace = function(f) {
  const halfedges = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  let he = f.halfedge;
  const first = he;
  do {
    he = he.next;
    halfedges.push(he);
  } while(he.id != first.id);
  // ----------- STUDENT CODE END ------------

  return halfedges;
};

// Return all faces adjacent to input face, not
// including input face.
Mesh.prototype.facesOnFace = function(f) {
  const faces = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  let he = f.halfedge;
  const first = he;
  do {
    he = he.next;
    faces.push(he.opposite.face);
  } while(he.id != first.id);
  // ----------- STUDENT CODE END ------------

  return faces;
};

// Return one-ring neighbors of input vertex, not
// including the input vertex itself
Mesh.prototype.verticesOnVertex = function(v) {
  const vertices = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  let he = v.halfedge;
  const first = he;
  do {
    he = he.opposite.next;
    vertices.push(he.vertex);
  } while (he.id != first.id);
  // ----------- STUDENT CODE END ------------

  return vertices;
};

// Return all halfedges that point away from v
Mesh.prototype.edgesOnVertex = function(v) {
  const halfedges = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  let he = v.halfedge;
  const first = he;
  do {
    he = he.opposite.next;
    halfedges.push(he);
  } while(he.id != first.id);
  // ----------- STUDENT CODE END ------------

  return halfedges;
};

// Return all faces that include v as a vertex.
Mesh.prototype.facesOnVertex = function(v) {
  const faces = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  let he = v.halfedge;
  const first = he;
  do {
    he = he.opposite.next;
    faces.push(he.face);
  } while(he.id != first.id);
  // ----------- STUDENT CODE END ------------

  return faces;
};

// Return the vertices that form the endpoints of a given edge
Mesh.prototype.verticesOnEdge = function(e) {
  const vertices = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 2 lines of code.
  vertices.push(e.vertex);
  vertices.push(e.opposite.vertex);
  // ----------- STUDENT CODE END ------------

  return vertices;
};

// Return the faces that include a given edge
Mesh.prototype.facesOnEdge = function(e) {
  const faces = [];
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 2 lines of code.
  faces.push(e.face);
  faces.push(e.opposite.face);
  // ----------- STUDENT CODE END ------------
  return faces;
};

// Return the edge pointing from v1 to v2
Mesh.prototype.edgeBetweenVertices = function(v1, v2) {
  let out_he = undefined;
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 11 lines of code.
  let he = this.edgesOnVertex(v1);
  for(let i = 0; i < he.length; i++) {
    if(he[i].vertex.position.equals(v2.position)){
      out_he = he[i];
    }
  }
  // ----------- STUDENT CODE END ------------
  return out_he;
};

////////////////////////////////////////////////////////////////////////////////
// Analysis
////////////////////////////////////////////////////////////////////////////////
function triangleArea(v1, v2, v3){
  const triArea = new THREE.Triangle(v1.position, v2.position, v3.position).getArea();
  return (triArea);
}
// Return the surface area of a provided face f.
Mesh.prototype.calculateFaceArea = function(f) {
  let area = 0.0;
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 21 lines of code.
  const vertices = this.verticesOnFace(f);

  for(let i = 0; i < vertices.length-2; i++){
    area += triangleArea(vertices[i], vertices[i+1], vertices[i+2]);
  }
  // ----------- STUDENT CODE END ------------
  f.area = area;
  return area;
};

// Update the area attributes of all faces in the mesh
Mesh.prototype.calculateFacesArea = function() {
  for (let i = 0; i < this.faces.length; ++i) {
    this.calculateFaceArea(this.faces[i]);
  }
};

// Calculate the vertex normal at a given vertex,
// using the face normals of bordering faces, weighted by face area
Mesh.prototype.calculateVertexNormal = function(v) {
  const v_normal = new THREE.Vector3(0, 0, 0);
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 11 lines of code.
  const fs = this.facesOnVertex(v);
  for (let f of fs){
    let norm = f.normal.clone();
    let t = norm.multiplyScalar(f.area);
    v_normal.add(t);
  }
  v_normal.normalize();
  // ----------- STUDENT CODE END ------------
  v.normal = v_normal;
  return v_normal;
};

// update the vertex normals of every vertex in the mesh
Mesh.prototype.updateVertexNormals = function() {
  for (let i = 0; i < this.vertices.length; ++i) {
    this.calculateVertexNormal(this.vertices[i]);
  }
};

// compute the average length of edges touching v
Mesh.prototype.averageEdgeLength = function(v) {
  let avg = 0.0;

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  const vertices = this.verticesOnVertex(v);
  for (let vert of vertices){
    avg += vert.position.distanceTo(v.position);
  }
  avg /= vertices.length;
  // ----------- STUDENT CODE END ------------

  return avg;
};

////////////////////////////////////////////////////////////////////////////////
// Topology
////////////////////////////////////////////////////////////////////////////////

// Given a face in the shape of an arbitrary polygon,
// split that face so it consists only of several triangular faces. 
Mesh.prototype.triangulateFace = function(f) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 8 lines of code.
  // ----------- STUDENT CODE END ------------
};
