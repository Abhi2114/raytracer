#pragma once

#include <glm.hpp>

using glm::vec3;
using glm::normalize;
using glm::cross;

class Camera {

	// position and co-ordinate frame for the camera
	// center is the look at point
	vec3 position, center, viewDirection, right, down;

public:

	// default position at the centre with the camera co-ordinate frame aligned with the world
	// co-ordinate frame and the camera facing the negative z-direction
	Camera(vec3 position, vec3 center) : position(position), center(center)
	{
		// the y-axis in world co-ordinates
		vec3 y(0, 1, 0);

		// construct a co-ordinate frame for the camera
		viewDirection = normalize(center - position);
		right = normalize(cross(y, viewDirection));
		down = normalize(cross(right, viewDirection));
	}

	vec3 getPosition()		const { return position; }
	vec3 getViewDirection() const { return viewDirection; }
	vec3 getRight()			const { return right; }
	vec3 getDown()			const { return down; }
};