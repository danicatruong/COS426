"use strict";

var Filters = Filters || {};

////////////////////////////////////////////////////////////////////////////////
// General utility functions
////////////////////////////////////////////////////////////////////////////////

// Constrain val to the range [min, max]
function clamp(val, min, max) {
  /* Shorthand for:
   * if (val < min) {
   *   return min;
   * } else if (val > max) {
   *   return max;
   * } else {
   *   return val;
   * }
   */
  return ((val < min) ? min : ((val > max) ? max : val));
}

// extract vertex coordinates from a URL string
function stringToCoords( vertsString ) {
  var centers = [];
  var coordStrings = vertsString.split("x");
  var coordsSoFar = 0;
  for (var i = 0; i < coordStrings.length; i++) {
    var coords = coordStrings[i].split("y");
    var x = parseInt(coords[0]);
    var y = parseInt(coords[1]);
    if (!isNaN(x) && !isNaN(y)) {
      centers.push({x: x, y: y})
    }
  }

  return centers;
}

////////////////////////////////////////////////////////////////////////////////
// Filters
////////////////////////////////////////////////////////////////////////////////

// Fill the entire image with color
Filters.fillFilter = function( image, color ) {
  for (var x = 0; x < image.width; x++) {
    for (var y = 0; y < image.height; y++) {
      // uncomment this line to enable this function
      image.setPixel(x, y, color);
    }
  }
  return image;
};

// At each center, draw a solid circle with the specified radius and color
Filters.brushFilter = function( image, radius, color, vertsString ) {
  // centers is an array of (x, y) coordinates that each defines a circle center
  var centers = stringToCoords(vertsString);

  // draw a filled circle centered at every location in centers[].
  // radius and color are specified in function arguments.
  // ----------- STUDENT CODE BEGIN ------------
  for (var i = 0; i < centers.length; i++) {
    // get current center from array of centers
    let current_center = centers[i];
    // looping to get each pixel in the desired drawn circle
    for (var x = current_center.x - radius; x < current_center.x + radius; x++){
      for (var y = current_center.y - radius; y < current_center.y + radius; y++){
        //check to make sure the pixel is inside the cirlce using distance formula
        let x_values = current_center.x - x;
        let y_values = current_center.y - y;
        let distance = Math.sqrt((x_values * x_values) + (y_values * y_values));
        if (distance < radius){image.setPixel(x, y, color)}
      }
    }
  }
  // ----------- Our reference solution uses 10 lines of code.
  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce ('brushFilter is not implemented yet');

  return image;
};

/*
 * At each center, draw a soft circle with the specified radius and color.
 * Pixel opacity should linearly decrease with the radius from alpha_at_center to 0.
 */
Filters.softBrushFilter = function( image, radius, color, alpha_at_center, vertsString ) {
  // centers is an array of (x, y) coordinates that each defines a circle center
  var centers = stringToCoords(vertsString);

  // draw a filled circle with opacity equals to alpha_at_center at the center of each circle
  // the opacity decreases linearly along the radius and becomes zero at the edge of the circle
  // radius and color are specified in function arguments.
  // ----------- STUDENT CODE BEGIN ------------
  for (var i = 0; i < centers.length; i++) {
    // get current center from array of centers
    let current_center = centers[i];
    // looping to get each pixel in the desired drawn circle
    for (var x = current_center.x - radius; x < current_center.x + radius; x++){
      for (var y = current_center.y - radius; y < current_center.y + radius; y++){
        //get distance of current pixel in the radius of the brush
        let x_values = current_center.x - x;
        let y_values = current_center.y - y;
        let distance = Math.sqrt((x_values * x_values) + (y_values * y_values));

        //calculate opacity which linearly decreases
        let opacity = alpha_at_center * (radius - distance) / radius;

        //get image color info
        let red_pixel_color = image.getPixel(x,y).data[0];
        let green_pixel_color = image.getPixel(x,y).data[1];
        let blue_pixel_color = image.getPixel(x,y).data[2];

        //get brush color info
        let red_brush_color = color.data[0];
        let green_brush_color = color.data[1];
        let blue_brush_color = color.data[2];

        //use alpha blending formula 
        let red_new_color = red_pixel_color * (1 - opacity) + red_brush_color * opacity;
        let green_new_color = green_pixel_color * (1 - opacity) + green_brush_color * opacity;
        let blue_new_color = blue_pixel_color * (1 - opacity) + blue_brush_color * opacity;

        //make a new pixel for the new color
        let new_pixel = new Pixel(red_new_color, green_new_color, blue_new_color, undefined, undefined)

        if (distance < radius){image.setPixel(x, y, new_pixel)}
      
      }
    }
  }

  // ----------- Our reference solution uses 21 lines of code.
  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce ('softBrushFilter is not implemented yet');

  return image;
};

Filters.customFilter = function( image, value ) {
  // You can use this filter to do whatever you want, for example
  // trying out some new idea or implementing something for the
  // art contest.
  // Currently the 'value' argument will be 1 or whatever else you set
  // it to in the URL. You could use this value to switch between
  // a bunch of different versions of your code if you want to
  // code up a bunch of different things for the art contest.
  // ----------- STUDENT CODE BEGIN ------------
  for (var x = 0; x < image.width; x++) {
    for (var y = 0; y < image.height; y++) {
      // black pixel
      let new_pixel = new Pixel(0, 0, 0, undefined, undefined)
      // randomize chance to set pixel to black
      // increase chance using value slider
      if (Math.random() * value > 0.1){
        image.setPixel(x,y, new_pixel);
      }
    }
  }
  // ----------- Our reference solution uses 0 lines of code.
  // ----------- STUDENT CODE END ------------
  // Gui.alertOnce ('customFilter is not implemented yet');
  return image;
};
