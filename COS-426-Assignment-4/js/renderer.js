"use strict";
var Reflection = Reflection || {
  ambient: new Pixel(0, 0, 0),
  diffuse: new Pixel(1.0, 1.0, 1.0),
  specular: new Pixel(1.0, 1.0, 1.0),
  shininess: 20,
};

Reflection.phongReflectionModel = function(vertex, view, normal, lightPos, phongMaterial) {
  var color = new Pixel(0, 0, 0);
  normal.normalize();

  // diffuse
  var light_dir = new THREE.Vector3().subVectors(lightPos, vertex).normalize();
  var ndotl = normal.dot(light_dir);
  color.plus(phongMaterial.diffuse.copy().multipliedBy(ndotl));

  // ----------- Our reference solution uses 12 lines of code.
  // Ambient color
  color.plus(phongMaterial.ambient);

  //specular color
  var reflection = light_dir.reflect(normal);
  var viewNorm = view.clone().normalize();
  var vdotr = Math.max(0,-viewNorm.dot(reflection));
  var vdotrpowAlpha = Math.pow(vdotr, phongMaterial.shininess);
  color.plus(phongMaterial.specular.copy().multipliedBy(vdotrpowAlpha));

  // ----------- STUDENT CODE END ------------

  return color;
};

var Renderer = Renderer || {
  meshInstances: new Set(),
  width: 1000,
  height: 750,
  negNear: 0.3,
  negFar: 1000,
  fov: 45,
  lightPos: new THREE.Vector3(10, 10, -10),
  shaderMode: "",
  cameraLookAtVector: new THREE.Vector3(0, 0, 0),
  cameraPosition: new THREE.Vector3(0, 0, -10),
  cameraUpVector: new THREE.Vector3(0, -1, 0),
  cameraUpdated: true,
};

Renderer.updateCameraParameters = function() {
  this.camera.position.copy(this.cameraPosition);
  this.camera.up.copy(this.cameraUpVector);
  this.camera.lookAt(this.cameraLookAtVector);
};

Renderer.initialize = function() {
  this.buffer = new Image(this.width, this.height);
  this.zBuffer = [];

  // set camera
  this.camera = new THREE.PerspectiveCamera(
    this.fov,
    this.width / this.height,
    this.negNear,
    this.negFar
  );
  this.updateCameraParameters();

  this.clearZBuffer();
  this.buffer.display(); // initialize canvas
};

Renderer.clearZBuffer = function() {
  for (var x = 0; x < this.width; x++) {
    this.zBuffer[x] = new Float32Array(this.height);
    for (var y = 0; y < this.height; y++) {
      this.zBuffer[x][y] = 1; // z value is in [-1 1];
    }
  }
};

Renderer.addMeshInstance = function(meshInstance) {
  assert(meshInstance.mesh, "meshInstance must have mesh to be added to renderer");
  this.meshInstances.add(meshInstance);
};

Renderer.removeMeshInstance = function(meshInstance) {
  this.meshInstances.delete(meshInstance);
};

Renderer.clear = function() {
  this.buffer.clear();
  this.clearZBuffer();
  Main.context.clearRect(0, 0, Main.canvas.width, Main.canvas.height);
};

Renderer.displayImage = function() {
  this.buffer.display();
};

Renderer.render = function() {
  this.clear();

  var eps = 0.01;
  if (
    !(
      this.cameraUpVector.distanceTo(this.camera.up) < eps &&
      this.cameraPosition.distanceTo(this.camera.position) < eps &&
      this.cameraLookAtVector.distanceTo(Main.controls.target) < eps
    )
  ) {
    this.cameraUpdated = false;
    // update camera position
    this.cameraLookAtVector.copy(Main.controls.target);
    this.cameraPosition.copy(this.camera.position);
    this.cameraUpVector.copy(this.camera.up);
  } else {
    // camera's stable, update url once
    if (!this.cameraUpdated) {
      Gui.updateUrl();
      this.cameraUpdated = true; //update one time
    }
  }

  this.camera.updateMatrixWorld();
  this.camera.matrixWorldInverse.getInverse(this.camera.matrixWorld);

  // light goes with the camera, COMMENT this line for debugging if you want
  this.lightPos = this.camera.position;

  for (var meshInst of this.meshInstances) {
    var mesh = meshInst.mesh;
    if (mesh !== undefined) {
      for (var faceIdx = 0; faceIdx < mesh.faces.length; faceIdx++) {
        var face = mesh.faces[faceIdx];
        var verts = [mesh.vertices[face.a], mesh.vertices[face.b], mesh.vertices[face.c]];
        var vert_normals = [
          mesh.vertex_normals[face.a],
          mesh.vertex_normals[face.b],
          mesh.vertex_normals[face.c],
        ];

        // camera's view matrix = K * [R | t] where K is the projection matrix and [R | t] is the inverse of the camera pose
        var viewMat = new THREE.Matrix4().multiplyMatrices(
          this.camera.projectionMatrix,
          this.camera.matrixWorldInverse
        );

        Renderer.drawTriangle(verts, vert_normals, mesh.uvs[faceIdx], meshInst.material, viewMat);
      }
    }
  }

  this.displayImage();
};

Renderer.getPhongMaterial = function(uv_here, material) {
  var phongMaterial = {};
  phongMaterial.ambient = Reflection.ambient;

  if (material.diffuse === undefined || uv_here === undefined) {
    phongMaterial.diffuse = Reflection.diffuse;
  } else if (Pixel.prototype.isPrototypeOf(material.diffuse)) {
    phongMaterial.diffuse = material.diffuse;
  } else {
    // note that this function uses point sampling. it would be better to use bilinear
    // subsampling and mipmaps for area sampling, but this good enough for now...
    phongMaterial.diffuse = material.diffuse.getPixel(
      Math.floor(uv_here.x * material.diffuse.width),
      Math.floor(uv_here.y * material.diffuse.height)
    );
  }

  if (material.specular === undefined || uv_here === undefined) {
    phongMaterial.specular = Reflection.specular;
  } else if (Pixel.prototype.isPrototypeOf(material.specular)) {
    phongMaterial.specular = material.specular;
  } else {
    phongMaterial.specular = material.specular.getPixel(
      Math.floor(uv_here.x * material.specular.width),
      Math.floor(uv_here.y * material.specular.height)
    );
  }

  phongMaterial.shininess = Reflection.shininess;

  return phongMaterial;
};

Renderer.projectVerticesNaive = function(verts) {
  // this is a naive orthogonal projection that does not even consider camera pose
  var projectedVerts = [];

  var orthogonalScale = 5;
  for (var i = 0; i < 3; i++) {
    projectedVerts[i] = new THREE.Vector4(verts[i].x, verts[i].y, verts[i].z, 1.0);

    projectedVerts[i].x /= orthogonalScale;
    projectedVerts[i].y /= (orthogonalScale * this.height) / this.width;

    projectedVerts[i].x = (projectedVerts[i].x * this.width) / 2 + this.width / 2;
    projectedVerts[i].y = (projectedVerts[i].y * this.height) / 2 + this.height / 2;
  }

  return projectedVerts;
};

Renderer.projectVertices = function(verts, viewMat) {
  // Vector3/Vector4 array of projected vertices in screen space coordinates
  // (you still need z for z buffering)
  var projectedVerts = [];

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 12 lines of code.
  for (var i = 0; i < 3; i++) {      
    projectedVerts[i] = new THREE.Vector4(verts[i].x, verts[i].y, verts[i].z, 1.0);

    projectedVerts[i].applyMatrix4(viewMat);

    projectedVerts[i].x /= projectedVerts[i].w;
    projectedVerts[i].y /= projectedVerts[i].w;
    projectedVerts[i].z /= projectedVerts[i].w;
  
    if (projectedVerts[i].z < -1 || projectedVerts[i].z > 1){
      return undefined;
    }

    projectedVerts[i].x = (projectedVerts[i].x * this.width) / 2 + this.width / 2;
    projectedVerts[i].y = (projectedVerts[i].y * this.height) / 2 + this.height / 2;

  }

  // ----------- STUDENT CODE END ------------

  return projectedVerts;
};

Renderer.computeBoundingBox = function(projectedVerts) {
  // Compute the screen-space bounding box for the triangle defined in projectedVerts[0-2].
  // We will need to call this helper function in the shading functions
  // to loop over pixel locations in the bounding box for rasterization.

  var box = {};

  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 14 lines of code.
  box.minX = Math.floor(Math.min(projectedVerts[0].x, projectedVerts[1].x, projectedVerts[2].x));
  box.minY = Math.floor(Math.min(projectedVerts[0].y, projectedVerts[1].y, projectedVerts[2].y));
  box.maxX = Math.ceil(Math.max(projectedVerts[0].x, projectedVerts[1].x, projectedVerts[2].x));
  box.maxY = Math.ceil(Math.max(projectedVerts[0].y, projectedVerts[1].y, projectedVerts[2].y));

  //constrained by screen size
  if (box.minX < 0){
    box.minX = 0;
  }
  if (box.maxX > this.width){
    box.maxX  = this.width;
  }
  if (box.minY < 0){
    box.minY = 0;
  }
  if (box.maxY > this.height){
    box.maxY = this.height;
  }

  // ----------- STUDENT CODE END ------------

  return box;
};

Renderer.computeBarycentric = function(projectedVerts, x, y) {
  var triCoords = [];
  // (see https://fgiesen.wordpress.com/2013/02/06/the-barycentric-conspirac/)
  // return undefined if (x,y) is outside the triangle
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 15 lines of code.

  var a = projectedVerts[0];
  var b = projectedVerts[1];
  var c = projectedVerts[2];

  var wANom = (b.x - c.x) * (c.y - y) - (c.x - x) * (b.y - c.y);
  var wADenom = (a.x - c.x) * (b.y - c.y) - (b.x - c.x) * (a.y - c.y);

  var wA = wANom / wADenom;
  if(wA < 0){
    return undefined;
  }

  var wBNom = (a.x - c.x) * (c.y - y) - (c.x - x) * (a.y - c.y);
  var wBDenom = (b.x - c.x) * (a.y - c.y) - (a.x - c.x) * (b.y - c.y);
  var wB = wBNom / wBDenom;

  if(wB < 0){
    return undefined;
  }

  var wC = 1 - wA - wB;

  if(wC < 0){
    return undefined;
  }

  triCoords.push(wA, wB, wC);
  // ----------- STUDENT CODE END ------------
  return triCoords;
};

Renderer.drawTriangleWire = function(projectedVerts) {
  var color = new Pixel(1.0, 0, 0);
  for (var i = 0; i < 3; i++) {
    var va = projectedVerts[(i + 1) % 3];
    var vb = projectedVerts[(i + 2) % 3];

    var ba = new THREE.Vector2(vb.x - va.x, vb.y - va.y);
    var len_ab = ba.length();
    ba.normalize();
    // draw line
    for (var j = 0; j < len_ab; j += 0.5) {
      var x = Math.round(va.x + ba.x * j);
      var y = Math.round(va.y + ba.y * j);
      this.buffer.setPixel(x, y, color);
    }
  }
};

Renderer.drawTriangleFlat = function(verts, projectedVerts, normals, uvs, material) {
  // Flat shader
  // Color of each face is computed based on the face normal
  // (average of vertex normals) and face centroid.
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 52 lines of code.

  // calculate face normal
  var faceNormal = (normals[0].clone().add(normals[1]).add(normals[2])).divideScalar(3);

  //calculate face centroid
  var faceCentroid = (verts[0].clone().add(verts[1]).add(verts[2])).divideScalar(3);

  var box = this.computeBoundingBox(projectedVerts);

  for (var x = box.minX; x < box.maxX; x++){
    for (var y = box.minY; y < box.maxY; y++){

      var baryCoords = this.computeBarycentric(projectedVerts, x, y);
      

      if (baryCoords !== undefined){
        var zCoord = baryCoords[0] * projectedVerts[0].z + baryCoords[1] * projectedVerts[1].z + baryCoords[2] * projectedVerts[2].z;
        if (zCoord < this.zBuffer[x][y]) {
          var uv = uvs;

          if (uvs !== undefined){
            uv = uvs[0].clone().multiplyScalar(baryCoords[0]).add(uvs[1].clone().multiplyScalar(baryCoords[1])).add(uvs[2].clone().multiplyScalar(baryCoords[2]));
          }
    
          var phongMaterial = this.getPhongMaterial(uv, material);
    
          var color = Reflection.phongReflectionModel(faceCentroid, this.camera.position, faceNormal, this.lightPos, phongMaterial);
                    
          this.buffer.setPixel(x, y, color);
          this.zBuffer[x][y] = zCoord;
        }        
      }
    }
  }

  // ----------- STUDENT CODE END ------------
};

Renderer.drawTriangleGouraud = function(verts, projectedVerts, normals, uvs, material) {
  // Gouraud shader
  // Interpolate the color for each pixel in the triangle using the barycentric coordinate.
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 49 lines of code.

  var phongMat0 = this.getPhongMaterial(undefined, material);
  var phongMat1 = phongMat0;
  var phongMat2 = phongMat0;

  if (uvs !== undefined){
    phongMat0 = this.getPhongMaterial(uvs[0], material);
    phongMat1 = this.getPhongMaterial(uvs[1], material);
    phongMat2 = this.getPhongMaterial(uvs[2], material);
  }

  var colors = [
    Reflection.phongReflectionModel(verts[0], this.camera.position, normals[0], this.lightPos, phongMat0),
    Reflection.phongReflectionModel(verts[1], this.camera.position, normals[1], this.lightPos, phongMat1),
    Reflection.phongReflectionModel(verts[2], this.camera.position, normals[2], this.lightPos, phongMat2)
  ];

  var box = this.computeBoundingBox(projectedVerts);

  for (let x = box.minX; x < box.maxX; x++) {
    for (let y = box.minY; y < box.maxY; y++) {
      var baryCoords = this.computeBarycentric(projectedVerts, x, y);
      
      if (baryCoords !== undefined) {

        var zCoord = baryCoords[0] * projectedVerts[0].z + baryCoords[1] * projectedVerts[1].z + baryCoords[2] * projectedVerts[2].z;
        
        if (zCoord < this.zBuffer[x][y]) {
          var color0 = colors[0].copy().multipliedBy(baryCoords[0]);
          var color1 = colors[1].copy().multipliedBy(baryCoords[1]);
          var color2 = colors[2].copy().multipliedBy(baryCoords[2]);

          var color = color0.plus(color1).plus(color2);
          
          this.buffer.setPixel(x, y, color);
          this.zBuffer[x][y] = zCoord;
        }
      }
    }
  }

  // ----------- STUDENT CODE END ------------
};

Renderer.drawTrianglePhong = function(verts, projectedVerts, normals, uvs, material) {
  // Phong shader
  // (1) Basic Phong shader: Interpolate the normal and vertex for each pixel in the triangle
  //                         using the barycentric coordinate.
  // (2) Texture mapping: If uvs is provided, compute interpolated uv coordinates
  //                      and map the phong material texture (if available)
  //                      at the uv coordinates to the pixel location.
  // (3) XYZ normal mapping: If xyz normal texture exists for the material,
  //                         convert the RGB value of the XYZ normal texture at the uv coordinates
  //                         to a normal vector and apply it at the pixel location.
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 62 lines of code.

  var box = this.computeBoundingBox(projectedVerts);

  for (var x = box.minX; x < box.maxX; x++){
    for (var y = box.minY; y < box.maxY; y++){

      var baryCoords = this.computeBarycentric(projectedVerts, x, y);


      if (baryCoords !== undefined) {
        var zCoord = baryCoords[0] * projectedVerts[0].z + baryCoords[1] * projectedVerts[1].z + baryCoords[2] * projectedVerts[2].z;

        if (zCoord < this.zBuffer[x][y]) {

          // interpolate normal
          var faceNorm0 = normals[0].clone().multiplyScalar(baryCoords[0]);
          var faceNorm1 = normals[1].clone().multiplyScalar(baryCoords[1]);
          var faceNorm2 = normals[2].clone().multiplyScalar(baryCoords[2]);

          var faceNormal = faceNorm0.add(faceNorm1).add(faceNorm2).normalize();

          var uv = undefined;

          if (uvs !== undefined){
            uv = uvs[0].clone().multiplyScalar(baryCoords[0]).add(uvs[1].clone().multiplyScalar(baryCoords[1])).add(uvs[2].clone().multiplyScalar(baryCoords[2]));

            // XYZ normal mapping
            if (material.xyzNormal !== undefined) {
              var rgbColor = material.xyzNormal.getPixel(Math.floor(uv.x*material.xyzNormal.width), Math.floor(uv.y*material.xyzNormal.height));
              faceNormal = (new THREE.Vector3(2*rgbColor.r-1, 2*rgbColor.g-1, 2*rgbColor.b-1)).normalize();
            }
          }

          // interpolate face center
          var faceCent0 = verts[0].clone().multiplyScalar(baryCoords[0]);
          var faceCent1 = verts[1].clone().multiplyScalar(baryCoords[1]);
          var faceCent2 = verts[2].clone().multiplyScalar(baryCoords[2]);
          var faceCentroid = faceCent0.add(faceCent1).add(faceCent2);

          var phongMaterial = this.getPhongMaterial(uv, material);
          var color = Reflection.phongReflectionModel(faceCentroid, this.camera.position, faceNormal, this.lightPos, phongMaterial);
                   
          this.buffer.setPixel(x, y, color);
          this.zBuffer[x][y] = zCoord;
        }    
      }
    }
  }
  // ----------- STUDENT CODE END ------------
};

Renderer.drawTriangle = function(verts, normals, uvs, material, viewMat) {
  var projectedVerts = this.projectVertices(verts, viewMat);
  if (projectedVerts === undefined) {
    // not within near and far plane
    return;
  } else if (projectedVerts.length <= 0) {
    projectedVerts = this.projectVerticesNaive(verts);
  }

  switch (this.shaderMode) {
    case "Wire":
      this.drawTriangleWire(projectedVerts);
      break;
    case "Flat":
      this.drawTriangleFlat(verts, projectedVerts, normals, uvs, material);
      break;
    case "Gouraud":
      this.drawTriangleGouraud(verts, projectedVerts, normals, uvs, material);
      break;
    case "Phong":
    this.drawTrianglePhong(verts, projectedVerts, normals, uvs, material);
      break;
    default:
  }
};
