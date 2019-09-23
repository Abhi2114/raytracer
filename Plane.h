#pragma once

#include<glm.hpp>
#include "Object.h"
#include "Ray.h"
#include "Color.h"

using glm::vec3;
using glm::dot;
using glm::normalize;

class Plane : public Object {

	vec3 normal;		// gives the orientation of the plane
	double distance;    // distance from the origin
	Color color;

public:

	Plane(vec3 normal, double distance, Color color) :
		normal(normalize(normal)),
		distance(distance),
		color(color) {}

	// getters
	vec3 getNormalAt(vec3 intersection) const {
		// ignore the intersection point since the normal to the plane will always stay the same
		return normal;
	}
	double getDistance()		 const { return distance; }
	virtual Color getColor()	 const { return color; }

	virtual Type getType() { return Type::PLANE; }

	// intersection between the ray of light and the plane
	// return value: distance between the ray origin and the intersection point
	virtual double findIntersection(Ray ray) {

		vec3 rayDirection = ray.getDirection();

		double a = dot(ray.getDirection(), normal);

		if (a == 0) {
			// ray is parallel to the plane
			return -1;
		}
		else {
			double b = dot(normal, ray.getOrigin() - (float)distance * normal);
			return -1 * b / a;
		}

	}
};