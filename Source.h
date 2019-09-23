#pragma once

#include <glm.hpp>
#include "Color.h"

class Source {

public:

	Source() { }

	virtual vec3 getPosition() const { return vec3(); }
	virtual Color getColor()   const { return Color(1, 1, 1, 0); }
};