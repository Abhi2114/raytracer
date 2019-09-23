#pragma once

#include "Ray.h"
#include "Color.h"

// the types of objects that we will make
enum Type {
	OBJECT, PLANE, SPHERE, TRIANGLE
};

// the root of the object hierarchy, an abstract class
class Object {

public:

	Object() {}

	virtual Color getColor() const { return Color(0, 0, 0, 0); }

	virtual double findIntersection(Ray ray) { return 0; }

	virtual Type getType() { return OBJECT; }

	virtual vec3 getNormalAt(vec3 intersection) const { return vec3(); }
};