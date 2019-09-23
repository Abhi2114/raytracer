#pragma once

#include "Vertex.h"
#include <iostream>

struct Vertex;
struct Color;
class ShapeGenerator;

typedef unsigned int uint;
typedef unsigned short ushort;

class ShapeData {

	Vertex* vertices;
	uint numVertices;  // number of vertices = number of colors

	ushort* indices;
	uint numIndices;  // only 256 unique vertices

public:

	ShapeData() : vertices(nullptr), numVertices(0), 
		indices(nullptr), numIndices(0) {}

	size_t vertexBufferSize() const { return numVertices * sizeof(Vertex); }
		
	size_t indexBufferSize() const { return numIndices * sizeof(ushort); }

	Vertex* vertexData() const { return vertices; }

	ushort* indexData() const { return indices; }

	uint getNumIndices() const { return numIndices; }

	~ShapeData() {
		delete[] vertices;
		delete[] indices;
	}

	friend class ShapeGenerator;
};
