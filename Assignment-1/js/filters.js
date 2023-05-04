"use strict";

const Filters = {};

////////////////////////////////////////////////////////////////////////////////
// General utility functions
////////////////////////////////////////////////////////////////////////////////

// Hardcoded Pi value
// const pi = 3.14159265359;
const pi = Math.PI;

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
    return val < min ? min : val > max ? max : val;
}

// Extract vertex coordinates from a URL string
function stringToCoords(vertsString) {
    const centers = [];
    const coordStrings = vertsString.split("x");
    for (let i = 0; i < coordStrings.length; i++) {
        const coords = coordStrings[i].split("y");
        const x = parseInt(coords[0]);
        const y = parseInt(coords[1]);
        if (!isNaN(x) && !isNaN(y)) {
            centers.push({ x: x, y: y });
        }
    }

    return centers;
}

// Blend scalar start with scalar end. Note that for image blending,
// end would be the upper layer, and start would be the background
function blend(start, end, alpha) {
    return start * (1 - alpha) + end * alpha;
}

// ----------- STUDENT CODE BEGIN ------------
function warp(image, sourceLines, destLines, alpha, sampleMode){
    let newImg = new Image(image.width, image.height);
    
    const p = 0.5;
    const a = 0.01;
    const b = 2;

    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++){
            let xSum = 0;
            let ySum = 0;
            let weightSum = 0;
            for (let linePos = 0; linePos < destLines.length; linePos++){
                let curDestLine = destLines[linePos];
                let curSourceLine = sourceLines[linePos];

                // P = (destLines.x0, destLines.y0)
                // Q = (destLines.x1, destLines.y1)
                // P' = (sourceLines.x0, sourceLines.y0)
                // Q' = (sourceLines.x1, sourceLines.y1)
                // calculate u and v which are scalars
                const denom = ((curDestLine.x1 - curDestLine.x0)**2 + (curDestLine.y1 - curDestLine.y0)**2);
                const u = ((x - curDestLine.x0)*(curDestLine.x1 - curDestLine.x0) + (y - curDestLine.y0)*(curDestLine.y1 - curDestLine.y0)) / denom;
                const v = ((x - curDestLine.x0)*(curDestLine.y1 - curDestLine.y0) + (y - curDestLine.y0)*(-curDestLine.x1 + curDestLine.x0)) / Math.sqrt(denom);
                
                // calculate xPrime which is a unit vector 
                let xSource, ySource;
                const denomPrime = Math.sqrt((curSourceLine.x1 - curSourceLine.x0)**2 + (curSourceLine.y1 - curSourceLine.y0)**2);
                xSource = curSourceLine.x0 + u * (curSourceLine.x1 - curSourceLine.x0) + (v * (curSourceLine.y1 - curSourceLine.y0) / denomPrime);
                ySource = curSourceLine.y0 + u * (curSourceLine.y1 - curSourceLine.y0) + (v * (-curSourceLine.x1 + curSourceLine.x0) / denomPrime);
                
                // calculate displacement distance, d
                let d;
                if (u < 0) {
                    d = Math.sqrt((x - curDestLine.x0)**2 + (y - curDestLine.y0)**2);
                }
                else if (u > 1){
                    d = Math.sqrt((x - curDestLine.x1)**2 + (y - curDestLine.y1)**2);
                }
                else{
                    d = Math.abs(v);
                }
                let lengthPQ = Math.sqrt((curSourceLine.x1 - curSourceLine.x0)**2 + (curSourceLine.y1 - curSourceLine.y0)**2);
                let weight = (lengthPQ**p / (a+d))**b;

                xSum += xSource * weight;
                ySum += ySource * weight;
                weightSum += weight;
            }

            // calculate average
            let xAvg, yAvg;
            
            // check x bounds
            xAvg = Math.round(xSum / weightSum);
            if (x < 0){
                xAvg = 0;
            } 
            else if (x > image.width)  {
                xAvg = image.width - 1;
            }

            // check y bounds
            yAvg = Math.round(ySum / weightSum);
            if (y < 0) {
                yAvg = 0;
            }
            else if (y > image.height) {
                yAvg = image.height - 1;
            }


            // make new morphed pixel
            const pixel = Filters.samplePixel(image, xAvg, yAvg, sampleMode);
            pixel.a = alpha;
            newImg.setPixel(x, y, pixel);
        }
    }
    return newImg;
}

// ----------- STUDENT CODE END ------------

////////////////////////////////////////////////////////////////////////////////
// Filters
////////////////////////////////////////////////////////////////////////////////

// You've already implemented this in A0! Feel free to copy your code into here
Filters.fillFilter = function(image, color) {
    image.fill(color);
    return image;
};

// You've already implemented this in A0! Feel free to copy your code into here
Filters.brushFilter = function(image, radius, color, vertsString) {
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

// You've already implemented this in A0! Feel free to copy your code into here
Filters.softBrushFilter = function(image, radius, color, alpha_at_center, vertsString) {
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

// Ratio is a value in the domain [-1, 1]. When ratio is < 0, linearly blend the image
// with black. When ratio is > 0, linearly blend the image with white. At the extremes
// of -1 and 1, the image should be completely black and completely white, respectively.
Filters.brightnessFilter = function(image, ratio) {
    let alpha, dirLuminance;
    if (ratio < 0.0) {
        alpha = 1 + ratio;
        dirLuminance = 0; // blend with black
    } else {
        alpha = 1 - ratio;
        dirLuminance = 1; // blend with white
    }

    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);

            pixel.data[0] = alpha * pixel.data[0] + (1 - alpha) * dirLuminance;
            pixel.data[1] = alpha * pixel.data[1] + (1 - alpha) * dirLuminance;
            pixel.data[2] = alpha * pixel.data[2] + (1 - alpha) * dirLuminance;

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

// Reference at this:
//      https://en.wikipedia.org/wiki/Image_editing#Contrast_change_and_brightening
// value = (value - 0.5) * (tan ((contrast + 1) * PI/4) ) + 0.5;
// Note that ratio is in the domain [-1, 1]
Filters.contrastFilter = function(image, ratio) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 14 lines of code.
    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            // change each channel individually
            pixel.data[0] = (clamp(pixel.data[0]) - 0.5) * (Math.tan((ratio + 1) * Math.PI/4)) + 0.5;
            pixel.data[1] = (clamp(pixel.data[1]) - 0.5) * (Math.tan((ratio + 1) * Math.PI/4)) + 0.5;
            pixel.data[2] = (clamp(pixel.data[2]) - 0.5) * (Math.tan((ratio + 1) * Math.PI/4)) + 0.5;
            image.setPixel(x,y,pixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('contrastFilter is not implemented yet');
    return image;
};

// Note that the argument here is log(gamma)
Filters.gammaFilter = function(image, logOfGamma) {
    const gamma = Math.exp(logOfGamma);
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 9 lines of code.
    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            // change each channel individually
            pixel.data[0]= pixel.data[0]**gamma;
            pixel.data[1] = pixel.data[1]**gamma;
            pixel.data[2] = pixel.data[2]**gamma;
            image.setPixel(x,y,pixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('gammaFilter is not implemented yet');
    return image;
};

/*
* The image should be perfectly clear up to innerRadius, perfectly dark
* (black) at outerRadius and beyond, and smoothly increase darkness in the
* circular ring in between. Both are specified as multiples of half the length
* of the image diagonal (so 1.0 is the distance from the image center to the
* corner).
*
* Note that the vignette should still form a perfect circle!
*/
Filters.vignetteFilter = function(image, innerR, outerR) {
    // Let's ensure that innerR is at least 0.1 smaller than outerR
    innerR = clamp(innerR, 0, outerR - 0.1);
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 17 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('vignetteFilter is not implemented yet');
    return image;
};

/*
* You will want to build a normalized CDF of the L channel in the image.
*/
Filters.histogramEqualizationFilter = function(image) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 33 lines of code.
    // keep track of lumiance and occurance
    let luminanceMap = new Map();

    // pixels into hashmap by luminance value
    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++){
            let pixel = image.getPixel(x, y);
            let luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2];

            luminance = luminance.toFixed(2);

            if (luminanceMap.has(luminance)){
                luminanceMap.set(luminance, luminanceMap.get(luminance) + 1);
            }
            else{
                luminanceMap.set(luminance, 1);
            }
        }
    }

    //sort lumaninceMap
    let sortedLumMap = new Map([...luminanceMap.entries(0)].sort());

    let cdfMap = new Map();
    let cdf = 0;

    // calculate cdf
    for (const [lum, occurance] of sortedLumMap){
        cdf += occurance / (image.width * image.height);
        cdfMap.set(lum, cdf);
    }

    // equalize pixels
    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++){
            let pixel = image.getPixel(x, y);
            let luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2];

            luminance = luminance.toFixed(2);

            pixel = pixel.rgbToHsl();

            pixel.data[2] = (cdfMap.get(luminance) * (luminanceMap.size - 1)) / (luminanceMap.size - 1);
            pixel = pixel.hslToRgb();
            image.setPixel(x, y, pixel);
        }
    }

    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('histogramEqualizationFilter is not implemented yet');
    return image;
};

// Set each pixel in the image to its luminance
Filters.grayscaleFilter = function(image) {
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            const luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2];
            pixel.data[0] = luminance;
            pixel.data[1] = luminance;
            pixel.data[2] = luminance;

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

// Adjust each channel in each pixel by a fraction of its distance from the average
// value of the pixel (luminance).
// See: http://www.graficaobscura.com/interp/index.html
Filters.saturationFilter = function(image, ratio) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 13 lines of code.
    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            //get luminance 
            const luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2]
            // change each channel individually
            pixel.data[0]= clamp(pixel.data[0]) + (clamp(pixel.data[0]) - luminance) * ratio;
            pixel.data[1] = clamp(pixel.data[1]) + (clamp(pixel.data[1]) - luminance) * ratio;
            pixel.data[2] = clamp(pixel.data[2]) + (clamp(pixel.data[2]) - luminance) * ratio;
            image.setPixel(x,y,pixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('saturationFilter is not implemented yet');
    return image;
};

// Apply the Von Kries method: convert the image from RGB to LMS, divide by
// the LMS coordinates of the white point color, and convert back to RGB.
Filters.whiteBalanceFilter = function(image, white) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 23 lines of code.
    let lmsWhite = white.rgbToXyz().xyzToLms();
    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            let lms = pixel.rgbToXyz().xyzToLms();
            lms.data[0] = lms.data[0] / lmsWhite.data[0];
            lms.data[1] = lms.data[1] / lmsWhite.data[1];
            lms.data[2] = lms.data[2] / lmsWhite.data[2];
            let rgb = lms.lmsToXyz().xyzToRgb();
            image.setPixel(x,y,rgb);
        }
    }
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('whiteBalanceFilter is not implemented yet');
    return image;
};

// This is similar to the histogram filter, except here you should take the
// the CDF of the L channel in one image and
// map it to another
//
Filters.histogramMatchFilter = function(image, refImg) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 58 lines of code.
    // reference image section
    let refLumMap = new Map();

    // pixels into hashmap by luminance value
    for (let x = 0; x < refImg.width; x++){
        for (let y = 0; y < refImg.height; y++){
            let pixel = refImg.getPixel(x, y);
            let luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2];

            luminance = luminance.toFixed(2);

            if (refLumMap.has(luminance)){
                refLumMap.set(luminance, refLumMap.get(luminance) + 1);
            }
            else{
                refLumMap.set(luminance, 1);
            }
        }
    }

    //sort lumaninceMap
    let sortedRefLumMap = new Map([...refLumMap.entries(0)].sort());

    let refCdfMap = new Map();
    let refCdf = 0;

    // calculate cdf
    for (const [lum, occurance] of sortedRefLumMap){
        refCdf += occurance / (refImg.width * refImg.height);
        refCdfMap.set(lum, refCdf);
    }

    // image section
    let luminanceMap = new Map();

    // pixels into hashmap by luminance value
    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++){
            let pixel = image.getPixel(x, y);
            let luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2];

            luminance = luminance.toFixed(2);

            if (luminanceMap.has(luminance)){
                luminanceMap.set(luminance, luminanceMap.get(luminance) + 1);
            }
            else{
                luminanceMap.set(luminance, 1);
            }
        }
    }

    //sort lumaninceMap
    let sortedLumMap = new Map([...luminanceMap.entries(0)].sort());

    let cdfMap = new Map();
    let cdf = 0;

    // calculate cdf
    for (const [lum, occurance] of sortedLumMap){
        cdf += occurance / (image.width * image.height);
        cdfMap.set(lum, cdf);
    }

    // equalize pixels
    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++){
            let pixel = image.getPixel(x, y);
            let luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2];

            luminance = luminance.toFixed(2);

            pixel = pixel.rgbToHsl();

            let min = cdfMap.get(luminance);
            let closestLum = luminance;
            for (const [lum, cdf] of refCdfMap){
                if(Math.abs(cdfMap.get(luminance) - cdf) < min){
                    min = Math.abs(cdfMap.get(luminance) - cdf);
                    closestLum = lum;
                }
            }

            pixel.data[2] = closestLum * (luminanceMap.size - 1) / (luminanceMap.size - 1);
            pixel = pixel.hslToRgb();
            image.setPixel(x, y, pixel);
        }
    }

    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('histogramMatchFilter is not implemented yet');
    return image;
};

// Convolve the image with a gaussian filter.
// NB: Implement this as a seperable gaussian filter
Filters.gaussianFilter = function(image, sigma) {
    // note: this function needs to work in a new copy of the image
    //       to avoid overwriting original pixels values needed later
    // create a new image with the same size as the input image
    let newImg = image.createImg(image.width, image.height);
    // the filter window will be [-winR, winR] for a total diameter of roughly Math.round(3*sigma)*2+1;
    const winR = Math.round(sigma * 3);
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 58 lines of code.
    let newNewImg = image.createImg(image.width, image.height);
    let matrix = [];
    for (let i = 0; i < 2*winR; i++){
        matrix[i] = (1/(sigma*Math.sqrt(2*Math.PI))) * Math.exp(-(i-winR**2)/(2*sigma**2));
    }
    //1D Gaussian kernel horizontally
    for (let y = 0; y < image.height; y++){
        for (let x = 0; x < image.width; x++){
            let red = 0;
            let green = 0; 
            let blue = 0;
            let weight = 0;
            for (let pos = 0; pos < matrix.length; pos++){
                if(x+pos-matrix.length/2 >= 0 && x+pos-matrix.length/2 < image.width){
                    let pixel = image.getPixel(x+pos-matrix.length/2,y);
                    red = red + matrix[pos]*pixel.data[0];
                    green = green + matrix[pos]*pixel.data[1];
                    blue = blue + matrix[pos]*pixel.data[2];
                    weight = weight + matrix[pos];
                }
            }
            let newPixel = new Pixel(red/weight, green/weight, blue/weight);
            newNewImg.setPixel(x, y, newPixel);
        }
    }
    // 1D Gaussian kernel vertically
    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++){
            let red = 0;
            let green = 0; 
            let blue = 0;
            let weight = 0;
            for (let pos = 0; pos < matrix.length; pos++){
                if(y+pos-matrix.length/2 >= 0 && y+pos-matrix.length/2 < image.height){
                    let pixel = newNewImg.getPixel(x, y+pos-matrix.length/2);
                    red = red + matrix[pos]*pixel.data[0];
                    green = green + matrix[pos]*pixel.data[1];
                    blue = blue + matrix[pos]*pixel.data[2];
                    weight = weight + matrix[pos];
                }
            }
            let newPixel = new Pixel(red/weight, green/weight, blue/weight);
            newImg.setPixel(x, y, newPixel);
        }
    }
    // // ----------- STUDENT CODE END ------------
    // // Gui.alertOnce ('gaussianFilter is not implemented yet');
    return newImg;
};

/*
* First the image with the edge kernel and then add the result back onto the
* original image.
*/
Filters.sharpenFilter = function(image) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 33 lines of code.
    let matrix = [[-1, -1, -1], [-1, 9, -1], [-1, -1, -1]];
    let newImg = image.createImg(image.width, image.height);


    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++){
            let weight = 0;
            let newPixel = new Pixel(0,0,0);
            for (let posX = 0; posX < matrix.length; posX++){
                for (let posY = 0; posY < matrix[0].length; posY++){
                    let pixel = image.getPixel(x + posX - matrix.length/2, y + posY-matrix[0].length/2);
                    newPixel.data[0] += pixel.data[0] * matrix[posX][posY];
                    newPixel.data[1] += pixel.data[1] * matrix[posX][posY];
                    newPixel.data[2] += pixel.data[2] * matrix[posX][posY];
                    weight = weight + matrix[posX][posY];
                }
            }
            newPixel.data[0] /= weight;
            newPixel.data[1] /= weight;
            newPixel.data[2] /= weight;
            newImg.setPixel(x, y, newPixel);
        }
    }   
    image = newImg;
     // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('sharpenFilter is not implemented yet');
    return image;
};

/*
* Convolve the image with the edge kernel from class. You might want to define
* a convolution utility that convolves an image with some arbitrary input kernel
*
* For this filter, we recommend inverting pixel values to enhance edge visualization
*/
Filters.edgeFilter = function(image) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 57 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('edgeFilter is not implemented yet');
    return image;
};

// Set a pixel to the median value in its local neighbor hood. You might want to
// apply this seperately to each channel.
Filters.medianFilter = function(image, winR) {
    // winR: the window will be  [-winR, winR];
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 36 lines of code.
    let newImg = image.createImg(image.width, image.height);
    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++){
            let pixel = image.getPixel(x,y);
            
            let red = [];
            let green = [];
            let blue = [];

            for (let i = x - winR; i < x + winR; i++){
                for (let j = y - winR; j < y + winR; j++){
                    let newPixel = image.getPixel(i, j);
                    red.push(newPixel.data[0]);
                    green.push(newPixel.data[1]);
                    blue.push(newPixel.data[2]);
                }
            }
            red.sort();
            green.sort();
            blue.sort();

            const length = red.length;

            if(length%2 != 0){
                pixel.data[0] = red[length/2];
                pixel.data[1] = green[length/2];
                pixel.data[2] = blue[length/2];
            }
            else{
                pixel.data[0] = (red[length/2] + red[length/2 + 1])/2;
                pixel.data[1] = (green[length/2] + green[length/2 + 1])/2;
                pixel.data[2] = (blue[length/2] + blue[length/2 + 1])/2;
            }
            newImg.setPixel(x, y, pixel);
        }
    }
    image = newImg;
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('medianFilter is not implemented yet');
    return image;
};

// Apply a bilateral filter to the image. You will likely want to reference
// precept slides, lecture slides, and the assignments/examples page for help.
Filters.bilateralFilter = function(image, sigmaR, sigmaS) {
    // reference: https://en.wikipedia.org/wiki/Bilateral_filter
    // we first compute window size and preprocess sigmaR
    const winR = Math.round((sigmaR + sigmaS) * 1.5);
    sigmaR = sigmaR * (Math.sqrt(2) * winR);

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 53 lines of code.
    let newImg = image.createImg(image.width, image.height);

    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++){
            let xyPixel = image.getPixel(x,y);
            let newPixel = new Pixel(0,0,0);
            let weightSum = 0;
            for (let posX = 0; posX < winR*2; posX++){
                for (let posY = 0; posY < winR*2; posY++){
                    let pixel = image.getPixel(x + posX - winR, y + posY - winR);

                    let spatialDistance = -((posX - winR)**2 + (posY - winR)**2) / (2 * sigmaS**2);
                    let colorDistance = -(255**2)*((xyPixel.data[0] - pixel.data[0])**2 + (xyPixel.data[1] - pixel.data[1])**2 + (xyPixel.data[2] - pixel.data[2])**2) / (2*sigmaR**2);
                    let weight = Math.exp(spatialDistance + colorDistance);

                    newPixel.data[0] += pixel.data[0] * weight;
                    newPixel.data[1] += pixel.data[1] * weight;
                    newPixel.data[2] += pixel.data[2] * weight;

                    weightSum += weight;
                }
            }
            newPixel.data[0] /= weightSum;
            newPixel.data[1] /= weightSum;
            newPixel.data[2] /= weightSum;

            newPixel.clamp();
            newImg.setPixel(x, y, newPixel);
        }
    }   
    image = newImg;
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('bilateralFilter is not implemented yet');
    return image;
};

// Conver the image to binary
Filters.quantizeFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    // use center color
    for (let i = 0; i < image.height; i++) {
        for (let j = 0; j < image.width; j++) {
            const pixel = image.getPixel(j, i);
            for (let c = 0; c < 3; c++) {
                pixel.data[c] = Math.round(pixel.data[c]);
            }
            pixel.clamp();
            image.setPixel(j, i, pixel);
        }
    }
    return image;
};

// To apply random dithering, first convert the image to grayscale, then apply
// random noise, and finally quantize
Filters.randomFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 12 lines of code.
    for (let i = 0; i < image.height; i++) {
        for (let j = 0; j < image.width; j++) {
            const pixel = image.getPixel(j, i);
            for (let c = 0; c < 3; c++) {
                pixel.data[c] = pixel.data[c] + (Math.random() - 0.5);
                pixel.data[c] = Math.round(pixel.data[c]);
            }
            pixel.clamp();
            image.setPixel(j, i, pixel);
        }
    }
    image = Filters.grayscaleFilter(image);
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('randomFilter is not implemented yet');
    return image;
};

// Apply the Floyd-Steinberg dither with error diffusion
Filters.floydFilter = function(image) {
    // convert to grayscale
    image = Filters.grayscaleFilter(image);

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 27 lines of code.
    let quantizeError = 0;
    for (let i = 0; i < image.height; i++) {
        for (let j = 0; j < image.width; j++) {
            const pixel = image.getPixel(j, i);

            let aPixel = image.getPixel(j+1, i);
            let bPixel = image.getPixel(j-1, i+1);
            let cPixel = image.getPixel(j, i+1);
            let dPixel = image.getPixel(j+1, i+1);

            for (let c = 0; c < 3; c++) {
                let temp = pixel.data[c];
                pixel.data[c] = Math.round(pixel.data[c]);
                quantizeError = temp - pixel.data[c];
                
                aPixel.data[c] = aPixel.data[c] + (7/16)*quantizeError;
                bPixel.data[c] = bPixel.data[c] + (3/16)*quantizeError;
                cPixel.data[c] = cPixel.data[c] + (5/16)*quantizeError;
                dPixel.data[c] = dPixel.data[c] + (1/16)*quantizeError;
            }
            pixel.clamp();
            aPixel.clamp();
            bPixel.clamp();
            cPixel.clamp();
            dPixel.clamp();

            image.setPixel(j, i, pixel);
            image.setPixel(j+1, i, aPixel);
            image.setPixel(j-1, i+1, bPixel);
            image.setPixel(j, i+1, cPixel);
            image.setPixel(j+1, i+1, dPixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('floydFilter is not implemented yet');
    return image;
};

// Apply ordered dithering to the image. We recommend using the pattern from the
// examples page and precept slides.
Filters.orderedFilter = function(image) {
    // convert to gray scale
    image = Filters.grayscaleFilter(image);

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 31 lines of code.
    let d = [
        [15, 7, 13, 5],
        [3, 11, 1, 9],
        [12, 4, 14, 6],
        [0, 8, 2, 10]
    ];

    let m = 4;
    for (let x = 0; x < image.width; x++){
        for (let y = 0; y < image.height; y++){
            let i = x % m;
            let j = y % m;
            let pixel = image.getPixel(x,y);

            let err = (pixel.data[0] - Math.floor(pixel.data[0]));
            let threshold = (d[i][j] + 1) / (m**2 + 1);

            if (err > threshold){
                pixel.data[0] = Math.ceil(pixel.data[0]);
                pixel.data[1] = Math.ceil(pixel.data[1]);
                pixel.data[2] = Math.ceil(pixel.data[2]);
            }
            else{
                pixel.data[0] = Math.floor(pixel.data[0]);
                pixel.data[1] = Math.floor(pixel.data[1]);
                pixel.data[2] = Math.floor(pixel.data[2]);
            }
            pixel.clamp();
            image.setPixel(x, y, pixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('orderedFilter is not implemented yet');
    return image;
};

// Implement bilinear and Gaussian sampling (in addition to the basic point sampling).
// This operation doesn't appear on GUI and should be used as a utility function.
// Call this function from filters that require sampling (e.g. scale, rotate)
Filters.samplePixel = function(image, x, y, mode) {
    if (mode === "bilinear") {
        // ----------- STUDENT CODE BEGIN ------------
        // ----------- Our reference solution uses 21 lines of code.
   
        let ceilingX = Math.ceil(x);
        let ceilingY = Math.ceil(y);
        let floorX = Math.floor(x);
        let floorY = Math.floor(y);

        // corner cases
        // if (ceilingX >= image.width) ceilingX = image.width-1;
        // if (ceilingY > image.height) ceilingY = image.height-1;
        // if (floorX < 0) floorX = 0;
        // if (floorY < 0) floorY = 0;

        // combinations 
        let ccP = image.getPixel(ceilingX, ceilingY);
        let cfP = image.getPixel(ceilingX, floorY);
        let fcP = image.getPixel(floorX, ceilingY); 
        let ffP = image.getPixel(floorX, floorY);

        let a = Math.ceil(x)-x;
        let b = Math.ceil(y)-y;
        let pixel = new Pixel(
            (1 - a) * (1 - b) * ccP.data[0] + (1 - a)*b*cfP.data[0] + a*(1 - b) * fcP.data[0] + a * b * ffP.data[0],
            (1 - a) * (1 - b) * ccP.data[1] + (1 - a)*b*cfP.data[1] + a*(1 - b) * fcP.data[1] + a * b * ffP.data[1],
            (1 - a) * (1 - b) * ccP.data[2] + (1 - a)*b*cfP.data[2] + a*(1 - b) * fcP.data[2] + a * b * ffP.data[2],
            1, ccP.colorSpace);

        pixel.clamp();
        return pixel;

        // ----------- STUDENT CODE END ------------
        // Gui.alertOnce ('bilinear sampling is not implemented yet');
    } else if (mode === "gaussian") {
        // ----------- STUDENT CODE BEGIN ------------
        // ----------- Our reference solution uses 38 lines of code.
        let sigma = 1;
        const winR = Math.round(sigma * 3);
        let newPixel = new Pixel(0,0,0);
        let gaussianSum = 0;

        for (let xPos = -winR; xPos < winR; xPos++){
            for (let yPos = -winR; yPos < winR; yPos++){
                    if(xPos >= 0  && yPos >= 0 && xPos < image.width && yPos < image.height){
                        let pixel = image.getPixel(xPos + x, yPos + y);
                        let gaussian = Math.exp(-(xPos**2)/(2*sigma**2));

                        newPixel.data[0] += pixel.data[0]*gaussian;
                        newPixel.data[1] += pixel.data[1]*gaussian;
                        newPixel.data[2] += pixel.data[2]*gaussian;
                        gaussianSum += gaussian;
                    }
                
            }
        }
        newPixel.data[0] /= gaussianSum;
        newPixel.data[1] /= gaussianSum;
        newPixel.data[2] /= gaussianSum;

        newPixel.clamp();
        return newPixel;

        // ----------- STUDENT CODE END ------------
        // Gui.alertOnce ('gaussian sampling is not implemented yet');
    } else {
        // point sampling
        y = Math.max(0, Math.min(Math.round(y), image.height - 1));
        x = Math.max(0, Math.min(Math.round(x), image.width - 1));
        return image.getPixel(x, y);
    }
};

// Translate the image by some x, y and using a requested method of sampling/resampling
Filters.translateFilter = function(image, x, y, sampleMode) {
    // Note: set pixels outside the image to RGBA(0,0,0,0)
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 21 lines of code.
    let newImg = image.createImg(image.width, image.height);
    Filters.fillFilter(newImg, new Pixel("rgba(0,0,0,0)"));

    for (let i = 0; i < image.width; i++){
        for (let j = 0; j < image.height; j++){
            if (i-x >= 0 && j-y >= 0 && i-x < image.width && j-y < image.height){
                newImg.setPixel(i, j, Filters.samplePixel(image, i-x, j-y, sampleMode));
            }
        }
    }
    image = newImg;

    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('translateFilter is not implemented yet');
    return image;
};

// Scale the image by some ratio and using a requested method of sampling/resampling
Filters.scaleFilter = function(image, ratio, sampleMode) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 19 lines of code.
    let newImg = image.createImg(image.width, image.height);
    Filters.fillFilter(newImg, new Pixel("rgba(0,0,0,0)"));
    
    let width = Math.round(image.width*ratio);
    let height = Math.round(image.height*ratio);
    
    for (let x = 0; x < width; x++){
        for (let y = 0; y < height; y++){
            newImg.setPixel(x,y,Filters.samplePixel(image, x/ratio, y/ratio, sampleMode));
        }
    }
    image = newImg;

    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('scaleFilter is not implemented yet');
    return image;
};

// Rotate the image by some angle and using a requested method of sampling/resampling
Filters.rotateFilter = function(image, radians, sampleMode) {
    // Note: set pixels outside the image to RGBA(0,0,0,0)
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 29 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('rotateFilter is not implemented yet');
    return image;
};

// Swirl the filter about its center. The rotation of the swirl should be in linear increase
// along the radial axis up to radians
Filters.swirlFilter = function(image, radians, sampleMode) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 26 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('swirlFilter is not implemented yet');
    return image;
};

// Set alpha from luminance
Filters.getAlphaFilter = function(backgroundImg, foregroundImg) {
    for (let i = 0; i < backgroundImg.height; i++) {
        for (let j = 0; j < backgroundImg.width; j++) {
            const pixelBg = backgroundImg.getPixel(j, i);
            const pixelFg = foregroundImg.getPixel(j, i);
            const luminance =
            0.2126 * pixelFg.data[0] + 0.7152 * pixelFg.data[1] + 0.0722 * pixelFg.data[2];
            pixelBg.a = luminance;
            backgroundImg.setPixel(j, i, pixelBg);
        }
    }

    return backgroundImg;
};

// Composites the foreground image over the background image, using the alpha
// channel of the foreground image to blend two images.
Filters.compositeFilter = function(backgroundImg, foregroundImg) {
    // Assume the input images are of the same sizes.
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 14 lines of code.
    for (let x = 0; x < backgroundImg.width; x++){
        for (let y = 0; y < backgroundImg.height; y++){
            let backgroundPixel = backgroundImg.getPixel(x,y);
            let foregroundPixel = foregroundImg.getPixel(x,y);
            let alpha = foregroundPixel.a;
            if (alpha > 0){
                backgroundPixel.data[0] = alpha * foregroundPixel.data[0] + (1-alpha) * backgroundPixel.data[0];
                backgroundPixel.data[1] = alpha * foregroundPixel.data[1] + (1-alpha) * backgroundPixel.data[1];
                backgroundPixel.data[2] = alpha * foregroundPixel.data[2] + (1-alpha) * backgroundPixel.data[2];
            }
            backgroundPixel.clamp();
            backgroundImg.setPixel(x, y, backgroundPixel);
        }
    }
    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('compositeFilter is not implemented yet');
    return backgroundImg;
};

// Morph two images according to a set of correspondance lines
Filters.morphFilter = function(initialImg, finalImg, alpha, sampleMode, linesFile) {
    const lines = Parser.parseJson("images/" + linesFile);

    // The provided linesFile represents lines in a flipped x, y coordinate system
    //  (i.e. x for vertical direction, y for horizontal direction).
    // Therefore we first fix the flipped x, y coordinates here.
    for (let i = 0; i < lines.initial.length; i++) {
        [lines.initial[i].x0, lines.initial[i].y0] = [lines.initial[i].y0, lines.initial[i].x0];
        [lines.initial[i].x1, lines.initial[i].y1] = [lines.initial[i].y1, lines.initial[i].x1];
        [lines.final[i].x0, lines.final[i].y0] = [lines.final[i].y0, lines.final[i].x0];
        [lines.final[i].x1, lines.final[i].y1] = [lines.final[i].y1, lines.final[i].x1];
    }

    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 114 lines of code.
    let warpedLines = new Array();
    for (let i = 0; i < lines.initial.length; i++){
        warpedLines.push({
            x0: lines.initial[i].x0 * (1-alpha) + lines.final[i].x0 * alpha,
            y0: lines.initial[i].y0 * (1-alpha) + lines.final[i].y0 * alpha,
            x1: lines.initial[i].x1 * (1-alpha) + lines.final[i].x1 * alpha,
            y1: lines.initial[i].y1 * (1-alpha) + lines.final[i].y1 * alpha,
        })
    }

    const initialWarp = warp(initialImg, lines.initial, warpedLines, 1, sampleMode);
    const finalWarp = warp(finalImg, lines.final, warpedLines, alpha, sampleMode);
    let image = Filters.compositeFilter(initialWarp, finalWarp);

    // ----------- STUDENT CODE END ------------
    // Gui.alertOnce ('morphFilter is not implemented yet');
    return image;
};

// Use k-means to extract a pallete from an image
Filters.paletteFilter = function(image, colorNum) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 89 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('paletteFilter is not implemented yet');
    return image;
};

// Read the following paper and implement your own "painter":
//      http://mrl.nyu.edu/publications/painterly98/hertzmann-siggraph98.pdf
Filters.paintFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 59 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('paintFilter is not implemented yet');
    return image;
};

/*
* Read this paper for background on eXtended Difference-of-Gaussians:
*      http://www.cs.princeton.edu/courses/archive/spring19/cos426/papers/Winnemoeller12.pdf
* Read this paper for an approach that develops a flow field based on a bilateral filter
*      http://www.cs.princeton.edu/courses/archive/spring19/cos426/papers/Kang09.pdf
*/
Filters.xDoGFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 70 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('xDoGFilter is not implemented yet');
    return image;
};

// You can use this filter to do whatever you want, for example
// trying out some new idea or implementing something for the
// art contest.
// Currently the 'value' argument will be 1 or whatever else you set
// it to in the URL. You could use this value to switch between
// a bunch of different versions of your code if you want to
// code up a bunch of different things for the art contest.
Filters.customFilter = function(image, value) {
    // ----------- STUDENT CODE BEGIN ------------
    // ----------- Our reference solution uses 0 lines of code.
    // ----------- STUDENT CODE END ------------
    Gui.alertOnce ('customFilter is not implemented yet');
    return image;
};
