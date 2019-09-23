#pragma once

#include <glm.hpp>

#include "Object.h"
#include "Color.h"

using glm::vec3;
using glm::dot;
using glm::sqrt;
using glm::normalize;

class Sphere : public Object {

	double radius;
	vec3 center;  // position of the center of the sphere
	Color color;
	bool transparent;
	float index;  // refractive index of the sphere

public:

	Sphere(double radius, vec3 center, Color color, bool transparent = false, float index = 0.0) :
		radius(radius), center(center), color(color), transparent(transparent), index(index)
	{}

	// getters
	double getRadius()		   const { return radius; }
	vec3 getCenter()		   const { return center; }
	bool isTransparent()	   const { return transparent; }
	float getIndex()		   const { return index; }

	virtual Color getColor()   const { return color; }

	virtual Type getType() { return Type::SPHERE; }

	virtual vec3 getNormalAt(vec3 intersection) const { return normalize(intersection - center); }

	virtual double findIntersection(Ray ray) {

		// get the origin and direction of the ray
		vec3 origin = ray.getOrigin();
		vec3 direction = ray.getDirection();

		// get the dot products that form the co-efficients of our quadratic equation
		double a = dot(direction, direction);
		double b = dot((float)2.0f * direction, origin - center);
		double c = dot(origin - center, origin - center) - radius * radius;

		// check the discriminant first
		double d = b * b - 4 * a * c;

		double distance = -1;  // if there is no intersection

		if (d > 0)  // there is an intersection
		{
			double r1 = (-b - sqrt(d)) / (2 * a) - 0.001;
			double r2 = (-b + sqrt(d)) / (2 * a) - 0.001;

			// return the smaller distance
			if (r1 > 0)  // r1 is the smallest + root
				distance = r1;
			else  // r2 is the smallest + root
				distance = r2;
		}

		return distance;
	}
};