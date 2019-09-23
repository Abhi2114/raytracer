#pragma once

#include <glm.hpp>

using glm::vec3;
using glm::mat3;
using glm::normalize;

class Ray {

	// a ray has an origin and a direction
	vec3 origin, direction;

public:

	// direct towards the x-axis if no args specified
	Ray(vec3 origin = vec3(),
		vec3 direction = vec3(1, 0, 0)) :
		origin(origin),
		direction(normalize(direction))
	{ }

	// rotate by angle and set the new origin to the given intersection point
	void rotate(float angle, vec3 intersection, vec3 axis) {
		// rotate the ray by the given angle
		direction = normalize(mat3(glm::rotate(angle, axis)) * direction);
		// set the new origin
		origin = intersection + 0.05f * direction;  // move the origin away from the surface of the object
	}

	vec3 getOrigin()	const { return origin; }
	vec3 getDirection() const { return direction; }
};