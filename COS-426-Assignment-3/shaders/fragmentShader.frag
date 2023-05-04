// set the precision of the float values (necessary if using float)
#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif
precision mediump int;

// flag for using soft shadows (set to 1 only when using soft shadows)
#define SOFT_SHADOWS 0

// define number of soft shadow samples to take
#define SOFT_SAMPLING 25.0

// define constant parameters
// EPS is for the precision issue
#define INFINITY 1.0e+12
#define EPS 1.0e-3

// define maximum recursion depth for rays
#define MAX_RECURSION 8

// define constants for scene setting
#define MAX_LIGHTS 10

// define texture types
#define NONE 0
#define CHECKERBOARD 1
#define MYSPECIAL 2

// define material types
#define BASICMATERIAL 1
#define PHONGMATERIAL 2
#define LAMBERTMATERIAL 3

// define reflect types - how to bounce rays
#define NONEREFLECT 1
#define MIRRORREFLECT 2
#define GLASSREFLECT 3

struct Shape {
  int shapeType;
  vec3 v1;
  vec3 v2;
  float rad;
};

struct Material {
  int materialType;
  vec3 color;
  float shininess;
  vec3 specular;

  int materialReflectType;
  float reflectivity;
  float refractionRatio;
  int special;
};

struct Object {
  Shape shape;
  Material material;
};

struct Light {
  vec3 position;
  vec3 color;
  float intensity;
  float attenuate;
};

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct Intersection {
  vec3 position;
  vec3 normal;
};

// uniform
uniform mat4 uMVMatrix;
uniform int frame;
uniform float height;
uniform float width;
uniform vec3 camera;
uniform int numObjects;
uniform int numLights;
uniform Light lights[MAX_LIGHTS];
uniform vec3 objectNorm;

// varying
varying vec2 v_position;

// find then position some distance along a ray
vec3 rayGetOffset(Ray ray, float dist) {
  return ray.origin + (dist * ray.direction);
}

// if a newly found intersection is closer than the best found so far, record
// the new intersection and return true; otherwise leave the best as it was and
// return false.
bool chooseCloserIntersection(float dist, inout float best_dist,
                              inout Intersection intersect,
                              inout Intersection best_intersect) {
  if (best_dist <= dist)
    return false;
  best_dist = dist;
  best_intersect.position = intersect.position;
  best_intersect.normal = intersect.normal;
  return true;
}

// put any general convenience functions you want up here
// ----------- STUDENT CODE BEGIN ------------
// ----------- Our reference solution uses 118 lines of code.

// determin if point is inside the bounding points of a box
bool inBounds(vec2 pos, vec2 pmin, vec2 pmax){
  return (pmin.x <= pos.x && pos.x <= pmax.x && pmin.y <= pos.y && pos.y <+ pmax.y);
}

// citation: https://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl
// link gotten from https://www.cs.princeton.edu/courses/archive/fall22/cos426/assignments/A3/#(3.0)-soft-shadows
float rand(vec2 co){
  return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}
// ----------- STUDENT CODE END ------------

// forward declaration
float rayIntersectScene(Ray ray, out Material out_mat,
                        out Intersection out_intersect);

// Plane
// this function can be used for plane, triangle, and box
float findIntersectionWithPlane(Ray ray, vec3 norm, float dist,
                                out Intersection intersect) {
  float a = dot(ray.direction, norm);
  float b = dot(ray.origin, norm) - dist;

  if (a < EPS && a > -EPS)
    return INFINITY;

  float len = -b / a;
  if (len < EPS)
    return INFINITY;

  intersect.position = rayGetOffset(ray, len);
  intersect.normal = norm;
  return len;
}

// Triangle
float findIntersectionWithTriangle(Ray ray, vec3 t1, vec3 t2, vec3 t3,
                                   out Intersection intersect) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 28 lines of code.
  // currently reports no intersection

  Intersection tempIntersect;
  vec3 v1 = t2 - t1;
  vec3 v2 = t3 - t1;
  vec3 norm = normalize(cross(v1, v2));

  float dist = findIntersectionWithPlane(ray, norm, dot(norm, t1), tempIntersect); 
  
  if(dist == INFINITY){
    return INFINITY;
  }

  // first side of the triangle
  vec3 v3 = t1 - tempIntersect.position;
  vec3 v4 = t2 - tempIntersect.position;
  vec3 v5 = t3 - tempIntersect.position;
  vec3 n1 = cross(v4, v3);
 
  if(dot(ray.direction, n1) < EPS){
    return INFINITY;
  }
  
  // second side of the triangle
  vec3 n2 = cross(v3, v5);

  if(dot(ray.direction, n2) < EPS){
    return INFINITY;
  }
 
  // third side of the triangle
  vec3 n3 = cross(v5, v4);
  
  if(dot(ray.direction, n3) < EPS){
    return INFINITY;
  }

  intersect.position = tempIntersect.position;
  intersect.normal = tempIntersect.normal;

  return dist;
  // return INFINITY;
  // ----------- STUDENT CODE END ------------
}

// Sphere
float findIntersectionWithSphere(Ray ray, vec3 center, float radius,
                                 out Intersection intersect) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 25 lines of code.
  vec3 L = center - ray.origin;
  float tca = dot(L, ray.direction);
  if (tca < EPS){
    return INFINITY;
  }

  float dSquare = dot(L, L) - (tca * tca);
  float radSquare = (radius * radius);

  if (dSquare > radSquare){
    return INFINITY;
  }

  float thc = sqrt(radSquare - dSquare);
  float t1 = tca - thc;
  float t2 = tca + thc;

  if (t1 > EPS){
    intersect.position = rayGetOffset(ray, t1);
    intersect.normal = normalize(intersect.position - center);
    return t1;
  } 

  else if(t2 > EPS){
    intersect.position = rayGetOffset(ray, t2);
    intersect.normal = normalize(intersect.position - center);
    return t2;
  }

  return INFINITY;
  // ----------- STUDENT CODE END ------------
}

// Box
float findIntersectionWithBox(Ray ray, vec3 pmin, vec3 pmax,
                              out Intersection out_intersect) {
  // ----------- STUDENT CODE BEGIN ------------
  // pmin and pmax represent two bounding points of the box
  // pmin stores [xmin, ymin, zmin] and pmax stores [xmax, ymax, zmax]
  float closestLen = INFINITY;
  Intersection intersect;

  //xy face
  vec3 norm = vec3(0.0, 0.0, 1.0);
  float len = findIntersectionWithPlane(ray, norm, pmin.z, intersect);
  if(len <= closestLen && inBounds(intersect.position.xy, pmin.xy, pmax.xy)){
    chooseCloserIntersection(len, closestLen, intersect, out_intersect);
  }
  len = findIntersectionWithPlane(ray, norm, pmax.z, intersect);
  if(len <= closestLen && inBounds(intersect.position.xy, pmin.xy, pmax.xy)){
    chooseCloserIntersection(len, closestLen, intersect, out_intersect);
  }

  //xz face
  norm = vec3(0.0, 1.0, 0.0);
  len = findIntersectionWithPlane(ray, norm, pmin.y, intersect);
  if(len <= closestLen && inBounds(intersect.position.xz, pmin.xz, pmax.xz)){
    chooseCloserIntersection(len, closestLen, intersect, out_intersect);
  }
  len = findIntersectionWithPlane(ray, norm, pmax.y, intersect);
  if(len <= closestLen && inBounds(intersect.position.xz, pmin.xz, pmax.xz)){
    chooseCloserIntersection(len, closestLen, intersect, out_intersect);
  }

  //yz face
  norm = vec3(1.0, 0.0, 0.0);
  len = findIntersectionWithPlane(ray, norm, pmin.x, intersect);
  if(len <= closestLen && inBounds(intersect.position.yz, pmin.yz, pmax.yz)){
    chooseCloserIntersection(len, closestLen, intersect, out_intersect);
  }
  len = findIntersectionWithPlane(ray, norm, pmax.x, intersect);
  if(len <= closestLen && inBounds(intersect.position.yz, pmin.yz, pmax.yz)){
    chooseCloserIntersection(len, closestLen, intersect, out_intersect);
  }
  return closestLen;

  // ----------- Our reference solution uses 44 lines of code.
  // return INFINITY;
  // ----------- STUDENT CODE END ------------
}

// Cylinder
float getIntersectOpenCylinder(Ray ray, vec3 center, vec3 axis, float len,
                               float rad, out Intersection intersect) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 33 lines of code.
  // currently reports no intersection
  return INFINITY;
  // ----------- STUDENT CODE END ------------
}

float getIntersectDisc(Ray ray, vec3 center, vec3 norm, float rad,
                       out Intersection intersect) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 18 lines of code.
  float dist = dot(center, norm);
  Intersection planeIntersect;
  float len = findIntersectionWithPlane(ray, norm, dist, planeIntersect);

  vec3 distFromCenter = planeIntersect.position - center;
  
  if (dot(norm, distFromCenter) < EPS && pow(length(distFromCenter), 2.0) < (rad * rad)) {
    intersect.position = planeIntersect.position;
    intersect.normal = planeIntersect.normal;
    return len;
  } 

  return INFINITY;
  // ----------- STUDENT CODE END ------------
}

float findIntersectionWithCylinder(Ray ray, vec3 center, vec3 apex,
                                   float radius,
                                   out Intersection out_intersect) {
  vec3 axis = apex - center;
  float len = length(axis);
  axis = normalize(axis);

  Intersection intersect;
  float best_dist = INFINITY;
  float dist;

  // -- infinite cylinder
  dist = getIntersectOpenCylinder(ray, center, axis, len, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  // -- two caps
  dist = getIntersectDisc(ray, center, -axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);
  dist = getIntersectDisc(ray, apex, axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);
  return best_dist;
}

// Cone
float getIntersectOpenCone(Ray ray, vec3 apex, vec3 axis, float len,
                           float radius, out Intersection intersect) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 45 lines of code.
  vec3 vd = ray.origin - apex;
  float phi = dot(vd, axis);
  float theta = dot(ray.direction, axis);
  float alpha = atan(radius/len);
  float cosAlpha2 = pow(cos(alpha), 2.0);
  float sinAlpha2 = pow(sin(alpha), 2.0);
  
  float a = pow(length(ray.direction - theta * axis), 2.0) * cosAlpha2 - theta * theta * sinAlpha2;
  float b = 2.0 * (dot(ray.direction - theta * axis, vd - phi * axis) * cosAlpha2 - theta * phi * sinAlpha2);
  float c = pow(length(vd - phi * axis), 2.0) * cosAlpha2 - phi * phi * sinAlpha2;
  
  float quadForm = b * b - 4.0 * a * c;
  if (quadForm < -EPS){
    return INFINITY;
  }
  quadForm = sqrt(quadForm);
  
  float closest = INFINITY;
  float t = (-b + quadForm) / (2.0 * a);

  if (t < EPS){
    return INFINITY;
  }
  
  vec3 q = rayGetOffset(ray, t);
  vec3 center = apex + axis * len;
  
  if (dot(axis, q - apex) > EPS && dot(axis, q - center) < -EPS) {
    closest = t;
    intersect.position = q;
  }
  
  t = (-b - quadForm) / (2.0*a);

  if (t < EPS){
    return INFINITY;
  }
  
  q = rayGetOffset(ray, t);

  if (t < closest && dot(axis, q-apex) > EPS && dot(axis, q-center) < -EPS) {
    closest = t;
    intersect.position = q;
  }
  
  if (closest != INFINITY) {
    vec3 apexToIntersect = intersect.position - apex;
    float y = dot(apexToIntersect, axis);
    float x = length(apexToIntersect - y * axis);
    vec3 xVector = normalize(apexToIntersect - y * axis);

    intersect.normal = normalize(-x * axis + y * xVector);
  }
  return closest;
  
  // ----------- STUDENT CODE END ------------
}

float findIntersectionWithCone(Ray ray, vec3 center, vec3 apex, float radius,
                               out Intersection out_intersect) {
  vec3 axis = center - apex;
  float len = length(axis);
  axis = normalize(axis);

  // -- infinite cone
  Intersection intersect;
  float best_dist = INFINITY;
  float dist;

  // -- infinite cone
  dist = getIntersectOpenCone(ray, apex, axis, len, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  // -- caps
  dist = getIntersectDisc(ray, center, axis, radius, intersect);
  chooseCloserIntersection(dist, best_dist, intersect, out_intersect);

  return best_dist;
}

vec3 calculateSpecialDiffuseColor(Material mat, vec3 posIntersection,
                                  vec3 normalVector) {
  // ----------- STUDENT CODE BEGIN ------------
  if (mat.special == CHECKERBOARD) {
    // ----------- Our reference solution uses 7 lines of code.
    float scale = 0.25;
    float x = floor(scale * posIntersection.x + EPS);
    float y = floor(scale * posIntersection.y + EPS);
    float z = floor(scale * posIntersection.z + EPS);

    if(mod(x+y+z, 2.0) == 0.0){
      return mat.color;
    }
    else{
      float avg = (mat.color.x + mat.color.y + mat.color.z) / 3.0;
      return (mat.color / avg) - mat.color;
    }
  }
    
  else if (mat.special == MYSPECIAL) {
    // vec3 noise = vec3(rand(posIntersection.xy), rand(posIntersection.yz), rand(posIntersection.zx));

    // float red = posIntersection.x * noise.x;
    // float green = posIntersection.y * noise.y;
    // float blue = posIntersection.z * noise.z;

    // float avg = (red + green + blue) / 3.0;
    // return mat.color / (noise * avg);
    // ----------- Our reference solution uses 5 lines of code.
  }

  // If not a special material, just return material color.
  return mat.color;
  // ----------- STUDENT CODE END ------------
}

vec3 calculateDiffuseColor(Material mat, vec3 posIntersection,
                           vec3 normalVector) {
  // Special colors
  if (mat.special != NONE) {
    return calculateSpecialDiffuseColor(mat, posIntersection, normalVector);
  }
  return vec3(mat.color);
}

// check if position pos in in shadow with respect to a particular light.
// lightVec is the vector from that position to that light -- it is not
// normalized, so its length is the distance from the position to the light
bool pointInShadow(vec3 pos, vec3 lightVec) {
  // ----------- STUDENT CODE BEGIN ------------
  Ray ray = Ray(pos, normalize(lightVec));
  Material mat;
  Intersection intersect;

  float len = rayIntersectScene(ray, mat, intersect);

  if (EPS < length(lightVec) - len){
    return true;
  }

  // ----------- Our reference solution uses 15 lines of code.
  return false;
  // ----------- STUDENT CODE END ------------
}

// use random sampling to compute a ratio that represents the
// fractional contribution of the light to the position pos.
// lightVec is the vector from that position to that light -- it is not
// normalized, so its length is the distance from the position to the light
float softShadowRatio(vec3 pos, vec3 lightVec) {
  // ----------- STUDENT CODE BEGIN ------------
  // ----------- Our reference solution uses 19 lines of code.
  float sampleSize = SOFT_SAMPLING * SOFT_SAMPLING;
  float pi = 3.1415926535897932384626433832795;

  float seenTheLight = 0.0;
  Material hitMaterial;
  Intersection intersect;
  Ray ray;
  ray.origin = pos;
  
  for (float i = 0.0; i < SOFT_SAMPLING; i++) {
    for (float j = 0.0; j < SOFT_SAMPLING; j++) {
      vec2 noise = vec2( (i + rand(pos.xy + vec2(i,j)))/SOFT_SAMPLING, (j + rand(pos.yz + vec2(i,j)))/SOFT_SAMPLING);
      float theta = noise.x * 2.0 * pi;
      float u = noise.y * 2.0 - 1.0;
      float num = sqrt(1.0 - u * u);
      vec3 intersectPoint =  vec3(num * cos(theta), num * sin(theta), u);
      
      vec3 offsetLight = lightVec + intersectPoint;
      ray.direction = normalize(offsetLight);
      float dist = rayIntersectScene(ray, hitMaterial, intersect);
      if (length(offsetLight) - dist < EPS) seenTheLight++;
    }
  }
  
  return seenTheLight / sampleSize;

  // return 0.0;
  // ----------- STUDENT CODE END ------------
}

vec3 getLightContribution(Light light, Material mat, vec3 posIntersection,
                          vec3 normalVector, vec3 eyeVector, bool phongOnly,
                          vec3 diffuseColor) {
  vec3 lightVector = light.position - posIntersection;


  float ratio = 1.0; // default to 1.0 for hard shadows
  if (SOFT_SHADOWS == 1) {
    // if using soft shadows, call softShadowRatio to determine
    // fractional light contribution
    ratio = softShadowRatio(posIntersection, lightVector);
  }
  else {
    // check if point is in shadow with light vector
    if (pointInShadow(posIntersection, lightVector)) {
      return vec3(0.0, 0.0, 0.0);
    }
  }

  // Slight optimization for soft shadows
  if (ratio < EPS) {
    return vec3(0.0, 0.0, 0.0);
  }


  // normalize the light vector for the computations below
  float distToLight = length(lightVector);
  lightVector /= distToLight;

  if (mat.materialType == PHONGMATERIAL ||
      mat.materialType == LAMBERTMATERIAL) {
    vec3 contribution = vec3(0.0, 0.0, 0.0);

    // get light attenuation
    float attenuation = light.attenuate * distToLight;
    float diffuseIntensity =
        max(0.0, dot(normalVector, lightVector)) * light.intensity;

    // glass and mirror objects have specular highlights but no diffuse lighting
    if (!phongOnly) {
      contribution +=
          diffuseColor * diffuseIntensity * light.color / attenuation;
    }

    if (mat.materialType == PHONGMATERIAL) {
      // Start with just black by default (i.e. no Phong term contribution)
      vec3 phongTerm = vec3(0.0, 0.0, 0.0);
      // ----------- STUDENT CODE BEGIN ------------
      vec3 ks = mat.specular;
      float shininess = mat.shininess;
      vec3 reflection = reflect(normalize(lightVector), normalVector);
      phongTerm = ks * max(0.0, pow(dot(eyeVector, reflection), shininess)) * light.intensity * light.color/attenuation;
      // ----------- Our reference solution uses 4 lines of code.
      // ----------- STUDENT CODE END ------------
      contribution += phongTerm;
    }

    return ratio * contribution;
  } else {
    return ratio * diffuseColor;
  }
}

vec3 calculateColor(Material mat, vec3 posIntersection, vec3 normalVector,
                    vec3 eyeVector, bool phongOnly) {
  // The diffuse color of the material at the point of intersection
  // Needed to compute the color when accounting for the lights in the scene
  vec3 diffuseColor = calculateDiffuseColor(mat, posIntersection, normalVector);

  // color defaults to black when there are no lights
  vec3 outputColor = vec3(0.0, 0.0, 0.0);

  // Loop over the MAX_LIGHTS different lights, taking care not to exceed
  // numLights (GLSL restriction), and accumulate each light's contribution
  // to the point of intersection in the scene.
  // ----------- STUDENT CODE BEGIN ------------
  for (int i = 0; i < MAX_LIGHTS; i++){
    if (i >= numLights){
      break;
    }
    outputColor += getLightContribution(lights[i], mat, posIntersection, normalVector, eyeVector, phongOnly, diffuseColor);
  }
  return outputColor;
  // ----------- Our reference solution uses 9 lines of code.
  // Return diffuseColor by default, so you can see something for now.
  // return diffuseColor;
  // ----------- STUDENT CODE END ------------
}

// find reflection or refraction direction (depending on material type)
vec3 calcReflectionVector(Material material, vec3 direction, vec3 normalVector,
                          bool isInsideObj) {
  if (material.materialReflectType == MIRRORREFLECT) {
    return reflect(direction, normalVector);
  }
  // If it's not mirror, then it is a refractive material like glass.
  // Compute the refraction direction.
  // See lecture 13 slide (lighting) on Snell's law.
  
  // ----------- STUDENT CODE BEGIN ------------
  float eta =
      (isInsideObj) ? 1.0 / material.refractionRatio : material.refractionRatio;
  float thetaI = acos(dot(normalVector, -direction));
  float thetaR = asin(eta * sin(thetaI));
  vec3 snell = (eta * cos(thetaI) - cos(thetaR)) * normalVector + eta * direction;

  return snell;
  // ----------- Our reference solution uses 5 lines of code.
  // Return mirror direction by default, so you can see something for now.
  // return reflect(direction, normalVector);
  // ----------- STUDENT CODE END ------------
}

vec3 traceRay(Ray ray) {
  // Accumulate the final color from tracing this ray into resColor.
  vec3 resColor = vec3(0.0, 0.0, 0.0);

  // Accumulate a weight from tracing this ray through different materials
  // based on their BRDFs. Initially all 1.0s (i.e. scales the initial ray's
  // RGB color by 1.0 across all color channels). This captures the BRDFs
  // of the materials intersected by the ray's journey through the scene.
  vec3 resWeight = vec3(1.0, 1.0, 1.0);

  // Flag for whether the ray is currently inside of an object.
  bool isInsideObj = false;

  // Iteratively trace the ray through the scene up to MAX_RECURSION bounces.
  for (int depth = 0; depth < MAX_RECURSION; depth++) {
    // Fire the ray into the scene and find an intersection, if one exists.
    //
    // To do so, trace the ray using the rayIntersectScene function, which
    // also accepts a Material struct and an Intersection struct to store
    // information about the point of intersection. The function returns
    // a distance of how far the ray travelled before it intersected an object.
    //
    // Then, check whether or not the ray actually intersected with the scene.
    // A ray does not intersect the scene if it intersects at a distance
    // "equal to zero" or far beyond the bounds of the scene. If so, break
    // the loop and do not trace the ray any further.
    // (Hint: You should probably use EPS and INFINITY.)
    // ----------- STUDENT CODE BEGIN ------------
    Material hitMaterial;
    Intersection intersect;
    // ----------- Our reference solution uses 4 lines of code.
    float dist = rayIntersectScene(ray, hitMaterial, intersect);
    if(abs(dist) < EPS || abs(dist) >= INFINITY){
      break;
    }
    // ----------- STUDENT CODE END ------------

    // Compute the vector from the ray towards the intersection.
    vec3 posIntersection = intersect.position;
    vec3 normalVector    = intersect.normal;

    vec3 eyeVector = normalize(ray.origin - posIntersection);

    // Determine whether we are inside an object using the dot product
    // with the intersection's normal vector
    if (dot(eyeVector, normalVector) < 0.0) {
        normalVector = -normalVector;
        isInsideObj = true;
    } else {
        isInsideObj = false;
    }

    // Material is reflective if it is either mirror or glass in this assignment
    bool reflective = (hitMaterial.materialReflectType == MIRRORREFLECT ||
                       hitMaterial.materialReflectType == GLASSREFLECT);

    // Compute the color at the intersection point based on its material
    // and the lighting in the scene
    vec3 outputColor = calculateColor(hitMaterial, posIntersection,
      normalVector, eyeVector, reflective);

    // A material has a reflection type (as seen above) and a reflectivity
    // attribute. A reflectivity "equal to zero" indicates that the material
    // is neither reflective nor refractive.

    // If a material is neither reflective nor refractive...
    // (1) Scale the output color by the current weight and add it into
    //     the accumulated color.
    // (2) Then break the for loop (i.e. do not trace the ray any further).
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 4 lines of code.
    if (abs(hitMaterial.reflectivity) < EPS){
      resColor += outputColor * resWeight;
      break;
    }
    // ----------- STUDENT CODE END ------------

    // If the material is reflective or refractive...
    // (1) Use calcReflectionVector to compute the direction of the next
    //     bounce of this ray.
    // (2) Update the ray object with the next starting position and
    //     direction to prepare for the next bounce. You should modify the
    //     ray's origin and direction attributes. Be sure to normalize the
    //     direction vector.
    // (3) Scale the output color by the current weight and add it into
    //     the accumulated color.
    // (4) Update the current weight using the material's reflectivity
    //     so that it is the appropriate weight for the next ray's color.
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 8 lines of code.

    vec3 nextBounce = calcReflectionVector(hitMaterial, ray.direction, normalVector, isInsideObj);
    ray.origin = posIntersection;
    ray.direction = normalize(nextBounce);

    resColor += outputColor * resWeight;
    resWeight += outputColor * hitMaterial.reflectivity;
    // ----------- STUDENT CODE END ------------
  }

  return resColor;
}

void main() {
  float cameraFOV = 0.8;
  vec3 direction = vec3(v_position.x * cameraFOV * width / height,
                        v_position.y * cameraFOV, 1.0);

  Ray ray;
  ray.origin = vec3(uMVMatrix * vec4(camera, 1.0));
  ray.direction = normalize(vec3(uMVMatrix * vec4(direction, 0.0)));

  // trace the ray for this pixel
  vec3 res = traceRay(ray);

  // paint the resulting color into this pixel
  gl_FragColor = vec4(res.x, res.y, res.z, 1.0);
}
