uniform sampler2D Texture;

uniform int scrWidth;
uniform int scrHeight;

const float sigma = 2.0;
const float pi = 3.14159265;
const float numBlurPixelsPerSide = 2.0;

const vec2 vertical = vec2(0.0, 1.0);
const vec2 horizontal = vec2(1.0, 0.0);

vec4 gaussianBlur(sampler2D tex, vec2 tCoord, vec2 multiplyVec, float blurSize) {
	// Incremental Gaussian Coefficent Calculation (See GPU Gems 3 pp. 877 - 889)
	vec3 incrementalGaussian;
  	incrementalGaussian.x = 1.0 / (sqrt(2.0 * pi) * sigma);
  	incrementalGaussian.y = exp(-0.5 / (sigma * sigma));
  	incrementalGaussian.z = incrementalGaussian.y * incrementalGaussian.y;

  	vec4 avgValue = vec4(0.0, 0.0, 0.0, 0.0);
  	float coefficientSum = 0.0;

  	// Take the central sample first...
  	avgValue += texture2D(tex, tCoord) * incrementalGaussian.x;
  	coefficientSum += incrementalGaussian.x;
  	incrementalGaussian.xy *= incrementalGaussian.yz;

  	// Go through the remaining 8 vertical samples (4 on each side of the center)
  	for (float i = 1.0; i <= numBlurPixelsPerSide; i++) {
  		vec2 offset = vec2((i * blurSize * multiplyVec/float(scrWidth)));
  		if (tCoord.x - offset.x < 0 || tCoord.x + offset.x > 1 ||
  			tCoord.y - offset.y < 0 || tCoord.y + offset.y > 1) {
  			avgValue += 2 * texture2D(tex, tCoord) * incrementalGaussian.x;
  		} else {
  			avgValue += texture2D(tex, tCoord - offset) * incrementalGaussian.x;  		         
	    	avgValue += texture2D(tex, tCoord + offset) * incrementalGaussian.x;
	    }
	             
	    coefficientSum += 2.0 * incrementalGaussian.x;
	    incrementalGaussian.xy *= incrementalGaussian.yz;
  	}
  	
  	return vec4(avgValue / coefficientSum);
}

void main() {
	vec2 tCoord   = vec2(gl_FragCoord.x/float(scrWidth), gl_FragCoord.y/float(scrHeight));
		
	float blurSize = 2.0;
  	gl_FragColor = gaussianBlur(Texture, tCoord, vertical, blurSize);
}
