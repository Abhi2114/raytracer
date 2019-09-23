#pragma once

#include <glm.hpp>
#include "Object.h"
#include "Color.h"

using glm::normalize;
using glm::cross;
using glm::dot;
using glm::abs;

class Triangle : public Object {

	vec3 a, b, c; // the 3 vertices of the triangle
	Color color;
	double u, v, w;
	vec3 normA, normB, normC;  // the vertex normal

public:

	Triangle() {}

	// the order the vertices are specified in matters
	// also specify all the vertex normals
	Triangle(vec3 a, vec3 b, vec3 c, 
			 vec3 normA, vec3 normB, vec3 normC, 
			 Color color) : a(a), b(b), c(c), 
							normA(normA), normB(normB), normC(normC), 
							color(color) {
	}

	virtual Type getType() { return Type::TRIANGLE; }

	virtual Color getColor() const { return color; }

	// the fast moller trumbore ray-triangle intersection test
	virtual double findIntersection(Ray ray) {

		vec3 rayOrigin = ray.getOrigin();
		vec3 rayDirection = ray.getDirection();

		vec3 e1 = b - a;
		vec3 e2 = c - a;

		// get the barycentric co-ordinates of the intersection point
		u = dot((rayOrigin - a), cross(rayDirection, e2)) / dot(e1, cross(rayDirection, e2));
		v = dot(rayDirection, (cross(rayOrigin - a, e1)) / dot(e1, cross(rayDirection, e2)));
		w = 1 - u - v;

		if ((u < 0) || (u > 1) || (v < 0) || (u + v > 1))
			return -1;
		
		return dot(e2, cross(rayOrigin - a, e1)) / dot(e1, cross(rayDirection, e2));
	}

	// do some kickass phong shading for a smoother shading effect
	virtual vec3 getNormalAt(vec3 intersection) const {
		// compute the interpolated normal at this point
		return (float)w * normA + (float)u * normB + (float)v * normC;
	}
};
