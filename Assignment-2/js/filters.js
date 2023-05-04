var Filters = Filters || {};

// Space for your helper functions
// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 105 lines of code.
// ----------- STUDENT CODE END ------------

// Translate all selected vertices in the mesh by the given x,y,z offsets.
Filters.translation = function(mesh, x, y, z) {
  // console.log("===== DEBUG START =====");
  // console.log(mesh);
  // const v1Index = 1;
  // const v1 = mesh.vertices[v1Index];
  // const v2Index = 2;
  // const v2 = mesh.vertices[v2Index];
  // const face3Index = 3;
  // const face3 = mesh.faces[face3Index];
  // const he1Index = 1;
  // const he1 = mesh.halfedges[he1Index];

  // console.log("== Traversal ==");
  // console.log(`edgesOnFace(${face3Index})`)
  // console.log(mesh.edgesOnFace(face3));
  // console.log(`facesOnFace(${face3Index})`);
  // console.log(mesh.facesOnFace(face3));
  // console.log(`verticesOnVertex(${v1Index})`);
  // console.log(mesh.verticesOnVertex(v1));
  // console.log(`edgesOnVertex(${v1Index})`);
  // console.log(mesh.edgesOnVertex(v1));
  // console.log(`facesOnVertex(${v1Index})`);
  // console.log(mesh.facesOnVertex(v1));
  // console.log(`verticesOnEdge(${he1Index})`);
  // console.log(mesh.verticesOnEdge(he1));
  // console.log(`facesOnEdge(${he1Index})`);
  // console.log(mesh.facesOnEdge(he1));
  // console.log(`edgeBetweenVertices(${v1Index},${v2Index})`);
  // console.log(mesh.edgeBetweenVertices(v1, v2));

  // console.log("== Analysis ==");
  // console.log(`calculateFacesArea(${face3Index})`);
  // console.log(mesh.calculateFaceArea(face3));
  // console.log(`averageEdgeLength(${v1Index})`);
  // console.log(mesh.averageEdgeLength(v1));

  // console.log("===== DEBUG END =====");

  const t = new THREE.Vector3(x, y, z);

  const verts = mesh.getModifiableVertices();

  const n_vertices = verts.length;
  for (let i = 0; i < n_vertices; ++i) {
    verts[i].position.add(t);
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Given x,y,z, the desired rotation around each axis, in radians,
// apply this rotation to all selected vertices in the mesh.
Filters.rotation = function(mesh, x, y, z) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 5 lines of code.
  const t = new THREE.Euler(x,y,z);

  for (let i = 0; i < verts.length; ++i) {
    verts[i].position.applyEuler(t);
  }
  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce("Rotation is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Uniformly scale the position of all selected vertices in the mesh
// by the provided scale factor s
Filters.scale = function(mesh, s) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 4 lines of code.
  for (let i = 0; i < verts.length; ++i) {
    verts[i].position.multiplyScalar(s);
  }
  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce("Scaling is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// estimate the per-vertex gaussian vurvature of the mesh at each vertex.
// set that vertex's color to some value based on its curvature value.
// (the precise mapping of curvature to color is left to you)
Filters.curvature = function(mesh) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 102 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("Curvature is not implemented yet");
};

// Apply a random offset to each selected vertex in the direction of its normal
// scale the random offset by the provided factor and by
// the average length of edges at that vertex
Filters.noise = function(mesh, factor) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 13 lines of code.
  for (let v of verts){
    const change = v.normal.clone();
    change.multiplyScalar(Math.random() * factor * mesh.averageEdgeLength(v));
    v.position.add(change);
  }
  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce("Noise is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Smooth the mesh using the specified weighting scheme.
// In the standard case, this is done using uniform Laplacian smoothing,
// by moving each vertex towards the average position of its neighbors.
//
// Arguments:
//  - mesh: the mesh to smooth
//  - iter: the number of iterations of smoothing to apply
//  - delta: a scaling factor for the amount of smoothing
//  - curvFlow: a bool. if true, use cotangent weights instead of uniform (requires triangular mesh)
//  - scaleDep: a bool. if true, scale offsets differently for each vertex (see spec.)
//  - implicit: a bool. if true, perform implicit smoothing (see spec.)
//
// Note that the reference solution calls a giant utility function so the line
// count is not terribly representative of the true solution
//
// For implicit, you will want to compute I - M*L*delta, where I is the identity
// matrix, M is a diagonal "mass" matrix, and L is a Laplacian matrix. Then
// you will want to call math.lup() on your result in order to decompose the
// matrix. Finally, call math.lusolve() to compute the X,Y, and Z positions of
// vertices. Note that the decomposition step allows for fast solving multiple
// times. It would be possible to replace a few of these steps with simple matrix
// inversion; however, matrix inversion is far slower and less numerically stable
//
Filters.smooth = function(mesh, iter, delta, curvFlow, scaleDep, implicit) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 16 lines of code.
  for (let i = 0; i < iter; i++){
    // matrix so you dont modify original values yet
    let offsets = [];

    for (const v of verts){
      const neighbors = mesh.verticesOnVertex(v);
      const pos = v.position.clone();
      pos.multiplyScalar(-neighbors.length);

      for (let i = 0; i < neighbors.length; i++){
        pos.add(neighbors[i].position);
      }
      
      pos.multiplyScalar(delta/neighbors.length);

      offsets.push(pos);
    }
    for (let i = 0; i < verts.length; i++){
      verts[i].position.add(offsets[i]);
    }
  }
  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce("Smooth is not implemented yet");
  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Sharpen the mesh by moving selected vertices away from the average position
// of their neighbors (i.e. Laplacian smoothing in the negative direction)
Filters.sharpen = function(mesh, iter, delta) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 9 lines of code.
  Filters.smooth(mesh, iter, -delta);
  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce("Sharpen is not implemented yet");
  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Move every selected vertex along its normal direction
// Scale the amount by the provided factor and average edge length at that vertex
Filters.inflate = function(mesh, factor) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 16 lines of code.
  let newPos = [];
  for (let i = 0; i < verts.length; i++) {
    let scale = factor*mesh.averageEdgeLength(verts[i]);
    let norm = verts[i].normal;
    norm.multiplyScalar(scale);
    newPos.push(norm);
  }
  for (let i = 0; i < verts.length; i++){
    verts[i].position.add(newPos[i]);
  }
  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce("Inflate is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// rotate selected vertices around the Y axis by an amount
// proportional to its Y value times the scale factor.
Filters.twist = function(mesh, factor) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 8 lines of code.

  for (let i = 0; i < verts.length; ++i) {
    let t = new THREE.Euler(0, verts[i].position.y*factor, 0);
    verts[i].position.applyEuler(t);
  }

  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce("Twist is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// warp a mesh using a nonlinear mapping of your choice
Filters.wacky = function(mesh, factor) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 3 lines of code.
  const verts = mesh.getModifiableVertices();
  let bound = Math.floor(verts.length*factor);
  let choosen = [];
  for (let i = 0; i < bound; i++){
    
    let index;

    while (true){
      let j = Math.floor(Math.random()*verts.length);
      if (!choosen.includes(j)){
        choosen.push(j);
        index = j;
        break
      }
    }

    const pos = verts[index].normal.clone();

    pos.multiplyScalar(Math.random())
    verts[index].position.add(pos);
  }
  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce("Wacky is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Convert the selected faces from arbitrary polygons into all triangles
Filters.triangulate = function(mesh) {
  const faces = mesh.getModifiableFaces();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 4 lines of code.
  for( let i = 0; i < faces.length; i++){
    let f = faces[i];
    while (mesh.verticesOnFace(f).length != 3){
      const v1 = f.halfedge.opposite.vertex;
      const v2 = f.halfedge.next.vertex;
      f = mesh.splitFaceMakeEdge(f, v1, v2, v1);
      if(faces[i].selected){
        f.selected = true;
      }
    }
  }

  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce("triangulate is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for splitEdgeMakeVert in mesh.js
Filters.splitEdge = function(mesh) {
  const verts = mesh.getSelectedVertices();

  if (verts.length === 2) {
    mesh.splitEdgeMakeVert(verts[0], verts[1], 0.5);
  } else {
    console.log("ERROR: to use split edge, select exactly 2 adjacent vertices");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for joinEdgeKillVert in mesh.js
Filters.joinEdges = function(mesh) {
  const verts = mesh.getSelectedVertices();

  if (verts.length === 3) {
    let v0 = verts[0],
      v1 = verts[1],
      v2 = verts[2];

    const he01 = mesh.edgeBetweenVertices(v0, v1);
    const he12 = mesh.edgeBetweenVertices(v1, v2);

    if (he01) {
      if (he12) {
        mesh.joinEdgeKillVert(verts[0], verts[1], verts[2]);
      } else {
        mesh.joinEdgeKillVert(verts[1], verts[0], verts[2]);
      }
    } else {
      if (he12) {
        mesh.joinEdgeKillVert(verts[0], verts[2], verts[1]);
      } else {
        console.log(
          "ERROR: to use join edge, select exactly 3 vertices such that one only has edges to the other two"
        );
      }
    }
  } else {
    console.log("ERROR: to use join edge, select exactly 3 vertices");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for splitFaceMakeEdge in mesh.js
Filters.splitFace = function(mesh) {
  const verts = mesh.getSelectedVertices();
  const faces = mesh.getModifiableFaces();

  if (verts.length === 2 && faces.length === 1) {
    mesh.splitFaceMakeEdge(faces[0], verts[0], verts[1]);
  } else {
    console.log("ERROR: to use split face, select exactly 1 face and 2 nonadjacent vertices on it");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// wrapper for joinFaceKillEdge in mesh.js
Filters.joinFaces = function(mesh) {
  const verts = mesh.getSelectedVertices();
  const faces = mesh.getModifiableFaces();

  if (verts.length === 2 && faces.length === 2) {
    mesh.joinFaceKillEdge(faces[0], faces[1], verts[0], verts[1]);
  } else {
    console.log(
      "ERROR: to use split face, select exactly 2 adjacent faces the 2 vertices between them"
    );
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// extrude the selected faces from the mesh in the direction of their normal
// vector, scaled by the provided factor.
// See the spec for more detail.
Filters.extrude = function(mesh, factor) {
  const faces = mesh.getModifiableFaces();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 32 lines of code.

  const faceNum = faces.length;

  for (let i = 0; i < faceNum; i++){
    const fNorm = faces[i].normal.clone();
    fNorm.multiplyScalar(factor);

    let nv = [];
    let he = mesh.edgesOnFace(faces[i]);

    for (let j = 0; j < he.length; j++){
      const v1 = he[j].vertex;
      const v2 = he[j].opposite.vertex;
      let newV = mesh.splitEdgeMakeVert(v1, v2, 0);
      nv.push(newV);
      mesh.splitFaceMakeEdge(he[j].opposite.face, v1, v2, newV, true);
    }
    
    for (let j = 0; j < nv.length; j++){
      nv[j].position.add(fNorm);
    }

    nv.push(nv[0]);
    for (let j = 1; j < nv.length; j++){
      mesh.splitFaceMakeEdge(faces[i], nv[j-1], nv[j], nv[j+1]);
    }

    for(let j = 0; j < he.length; j++){
      mesh.joinFaceKillEdgeSimple(he[j]);
    }

  }

  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce("Extrude is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Truncate the selected vertices of the mesh by "snipping off" corners
// and replacing them with faces. factor specifies the size of the truncation.
// See the spec for more detail.
Filters.truncate = function(mesh, factor) {
  const verts = mesh.getModifiableVertices();

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 64 lines of code.
  let edges = [];
  let edgeVertPos = [];
  const len = verts.length;

  for (const v of verts){
    edges.push(mesh.edgesOnVertex(v));
  }
  for (let i = 0; i < len; i++){
    edgeVertPos.push([]);
    for (let j = 0; j < edges[i].length; j++){
      edgeVertPos[i].push(edges[i][j].vertex.position.clone());
    }
  }

  for (let i = 0; i < len; i++){
    for (let j = 1; j < edges[i].length; j++){
      const newV = mesh.splitEdgeMakeVert(verts[i], edges[i][j].vertex, 0);

      edgeVertPos[i][j].multiplyScalar(factor);
      newV.position.multiplyScalar(1-factor);
      newV.position.add(edgeVertPos[i][j]);
    }

    edgeVertPos[i][0].multiplyScalar(factor);
    verts[i].position.multiplyScalar(1-factor);
    verts[i].position.add(edgeVertPos[i][0]);

    for (let j = 2; j < edges[i].length; j++){
      mesh.splitFaceMakeEdge(edges[i][j].face, edges[i][j-1].vertex, edges[i][j].vertex);
    }

    for (let j = 2; j < edges[i].length - 1; j++){
      mesh.joinFaceKillEdgeSimple(edges[i][j]);
    }
  }
  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce("Truncate is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Apply the bevel operation to the mesh, scaling the degree of bevelling by factor
Filters.bevel = function ( mesh, factor ) {

    var verts = mesh.getModifiableVertices();

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 104 lines of code.

    const faceLen = mesh.faces.length;
    Filters.truncate(mesh, factor);

    let oldEdge = [];

    for (let i = faceLen; i < mesh.faces.length; i++){
      const edges = mesh.edgesOnFace(mesh.faces[i]);

      for (const he of edges){
        let tempEdges = mesh.edgesOnVertex(he.vertex);
        //longest edge will always be the edge that isnt the truncated faces' edge
        for(const t of tempEdges){
          if(t != he.next && t != he.opposite){
            oldEdge.push(t);
          }
        }
        mesh.splitEdgeMakeVert(he.opposite.vertex, he.vertex);
      }
    }

    for (const edge of oldEdge){
      mesh.splitFaceMakeEdge(edge.face, edge.next.vertex, edge.opposite.next.opposite.next.vertex);
    }
    for (const edge of oldEdge){
      mesh.joinFaceKillEdgeSimple(edge);
    }
    for (const edge of oldEdge){
      mesh.joinEdgeKillVert(edge.next.vertex, edge.next.opposite.vertex, edge.next.opposite.next.vertex);
    }
  
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('Bevel is not implemented yet');

    mesh.calculateFacesArea();
    mesh.updateNormals();
};

// Split the longest edges in the mesh into shorter edges.
// factor is a float in [0,1]. it tells the proportion
// of the total number of edges in the mesh that should be split.
Filters.splitLong = function(mesh, factor) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 35 lines of code.
  // ----------- STUDENT CODE END ------------
  Gui.alertOnce("Split Long Edges is not implemented yet");

  mesh.calculateFacesArea();
  mesh.updateNormals();
};


// making midpoints
midpoints = function(mesh, faces, vertsLen){
  for (const f of faces) {
    for (const he of mesh.edgesOnFace(f)) {
      if (he.vertex.id < vertsLen && he.opposite.vertex.id < vertsLen)
        mesh.splitEdgeMakeVert(he.vertex, he.opposite.vertex);
    }
  }
}

//subdividing
subdivide = function(mesh, faces, vertsLen)
{
  const numFaces = faces.length;
  let newFaces = [];

  for (let i = 0; i < numFaces; i ++){
    newFaces[i] = faces[i];
  }

  for (let i = 0; i < numFaces; i++) {
    const f = faces[i];
    let he = f.halfedge.next;
    const first = he;
    const midpoints = [];

    // new vertex at end of the original vertex array
    if(faces[i].halfedge.vertex.id < vertsLen){
      do {
        const midpoint = he.vertex;
        midpoints.push(midpoint);
        he = he.next.next;
      } while (he != first);
    }

    midpoints.push(midpoints[0]); // allow looping at end of array
    for (let j = 1; j < midpoints.length; j++){
      newFaces.push( mesh.splitFaceMakeEdge(f, midpoints[j-1], midpoints[j]) );
    }
  }

  return newFaces;
}

// Triangulate a mesh, and apply triangular subdivision to its faces.
// Repeat for the specified number of levels.
Filters.triSubdiv = function(mesh, levels) {
  Filters.triangulate(mesh);

  for (let l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 43 lines of code.
    let vertsLen = mesh.vertices.length;

    midpoints(mesh, faces, vertsLen);
    const newFaces = subdivide(mesh, faces, vertsLen);
  
    // make sure new faces are selected if original face was selected
    if (faces != mesh.faces) {
      let fIDs = [];
      console.log(newFaces);
      for (let i = 0; i < newFaces.length; i++){
        fIDs.push(newFaces[i].id);

      }
      mesh.setSelectedFaces(fIDs);
    }   

    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce("Triangle subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

boundaries = function(mesh, faces, v) {
  // all faces are selected
  if (faces == mesh.faces) {
    return null;
  }

  const edges = mesh.edgesOnVertex(v);
  const boundary = [];
  edges.push(edges[0]);
  let selected = edges[0].face.selected;

  for (let i = 1; i < edges.length; i++) {
    if (edges[i].face.selected != selected) {
      selected = !selected;
      boundary.push(edges[i-1]);
      //each vertex can only be connected to two boundary edges per face
      if (boundary.length == 2){
        return boundary;
      }
    }
  }
  return null;
}

// Triangulate the mesh and apply loop subdivision to the faces
// repeat for the specified number of levels.
Filters.loop = function(mesh, levels) {
  Filters.triangulate(mesh);

  for (let l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 123 lines of code.
    const vertsLen = mesh.vertices.length;
    const boundFaces = [];
    const newVPos = Array(2 * vertsLen);

    const verts = mesh.getModifiableVertices();

    for (const v of verts) {
      const newPos = new THREE.Vector3();
      const boundary = boundaries(mesh, faces, v);

      let weight;
      
      let vertFaces = mesh.facesOnVertex(v);

      //boundary case
      if (boundary != null) {
        for (const f of vertFaces) {
          if (!f.selected && !boundFaces.includes(f)) {
            boundFaces.push(f);
          }
        }
        newPos.add(boundary[0].vertex.position);
        newPos.add(boundary[1].vertex.position);
        newPos.multiplyScalar(1/8);

        weight = 3/4;

      } else {
        const neighborhood = mesh.verticesOnVertex(v);
        
        let beta = 3/16;
        if(neighborhood.length != 3){
          beta = 3/(8*neighborhood.length);
        }

        for (const neigh of neighborhood) { 
          newPos.add(neigh.position); 
        }

        newPos.multiplyScalar(beta);

        // k = neighborhood.length, weights = 1-k(beta)
        weight = 7/16;
        if (neighborhood.length != 3){
          weight = 5/8;
        }
      }

      const center = v.position.clone();
      center.multiplyScalar(weight);
      newPos.add(center);

      newVPos[v.id] = newPos;
    }

    // make midpoints
    midpoints(mesh, faces, vertsLen);

    // midpoint positions
    for (let i = vertsLen; i < mesh.vertices.length; i++) {
      let midpoint = mesh.vertices[i];
      halfEdge = midpoint.halfedge.opposite;

      if (boundaries(mesh, faces, midpoint) == null) {
        // adjacent
        const adjPos = halfEdge.next.vertex.position.clone();
        adjPos.add(halfEdge.opposite.vertex.position);
        adjPos.multiplyScalar(3/8);

        // opposite 
        const opPos = halfEdge.next.next.next.vertex.position.clone();
        opPos.add(halfEdge.opposite.next.next.vertex.position);
        opPos.multiplyScalar(1/8);

        opPos.add(adjPos);
        newVPos[midpoint.id] = opPos;
      }
    }

    // make midpoints with boundary faces
    midpoints(mesh, boundFaces, vertsLen);

    // subdivide
    let newFaces = subdivide(mesh, faces, vertsLen);
    subdivide(mesh, boundFaces, vertsLen)

    // set vertices to newPos
    for (let i = 0; i < newVPos.length; i++){
      if(newVPos[i] != undefined){
        mesh.vertices[i].position = newVPos[i];
      }
    }

    if (faces != mesh.faces) {
      let fIDs = [];
      console.log(newFaces);
      for (let i = 0; i < newFaces.length; i++){
        fIDs.push(newFaces[i].id);
      }
      mesh.setSelectedFaces(fIDs);
    }   
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce("Triangle subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Requires a quad mesh. Apply quad subdivision to the faces of the mesh.
// Repeat for the specified number of levels.
Filters.quadSubdiv = function(mesh, levels) {
  for (let l = 0; l < levels; l++) {
    const faces = mesh.getModifiableFaces();
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 55 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("Quad subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// Apply catmull clark subdivision to the faces of the mesh.
// Repeat for the specified number of levels.
Filters.catmullClark = function(mesh, levels) {
  for (let l = 0; l < levels; l++) {
    const faces = mesh.faces;
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 102 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce("Catmull-Clark subdivide is not implemented yet");
  }

  mesh.calculateFacesArea();
  mesh.updateNormals();
};

// ================= internal functions =======================

// internal function for selecting faces in the form of a loop
Filters.procFace = function(mesh, f) {
  const faceFlags = new Array(mesh.faces.length);
  for (let i = 0; i < mesh.faces.length; i++) {
    faceFlags[i] = 0;
  }
  let sum = f.area;
  const start_he = f.halfedge.opposite.next;
  let curr_he = start_he;
  do {
    if (faceFlags[curr_he.face.id] > 0) {
      break;
    }
    sum += curr_he.face.area;
    curr_he.face.selected = true;
    faceFlags[curr_he.face.id]++;
    const last_he = curr_he;
    curr_he = curr_he.opposite.next;
    if (curr_he.face == f) {
      curr_he = last_he.next.opposite.next;
    }
  } while (true);
};

Filters.parseSelected = function(sel) {
  if (sel === undefined || sel.replace === undefined) {
    return [];
  }
  if (typeof sel === "number") {
    return [sel];
  }
  // sel = sel.replace(/[\(\)]/g,'');
  sel = sel.split(",");
  const parsedSel = [];
  for (let i = 0; i < sel.length; i++) {
    const idx = parseInt(sel[i]);
    if (!isNaN(idx)) {
      parsedSel.push(idx);
    }
  }
  return parsedSel;
};

// internal filter for updating selection
Filters.selection = function(mesh, vertIdxs, faceIdxs) {
  mesh.setSelectedVertices(Filters.parseSelected(vertIdxs));
  mesh.setSelectedFaces(Filters.parseSelected(faceIdxs));
};

// internal filter for setting display settings
Filters.displaySettings = function(
  mesh,
  showLabels,
  showHalfedge,
  shading,
  showVN,
  showFN,
  showGrid,
  showVertDots,
  showAxes,
  showVC,
  meshColor
) {
  Main.displaySettings.showIdLabels = showLabels;
  Main.displaySettings.wireframe = showHalfedge;
  Main.displaySettings.shading = shading;
  Main.displaySettings.showVN = showVN;
  Main.displaySettings.showFN = showFN;
  Main.displaySettings.showGrid = showGrid;
  Main.displaySettings.showVertDots = showVertDots;

  Main.displaySettings.showAxes = showAxes;
  Main.displaySettings.showVC = showVC;
  // Main.displaySettings.meshColor = meshColor;

  // Main.refreshDisplaySettings();
};
