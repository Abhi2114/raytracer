#pragma once

#include <glm.hpp>
#include "Source.h"
#include "Color.h"

using glm::vec3;

// model a point light, it has a position and a color
class Light : public Source {

	vec3 position;
	Color color;

public:

	Light(vec3 position, Color color) : position(position), color(color) {}

	// getters
	virtual vec3 getPosition() const { return position; }
	virtual Color getColor()   const { return color; }
};