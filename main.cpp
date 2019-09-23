#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <gtx/transform.hpp>

#include "Color.h"
#include "Ray.h"
#include "Camera.h"
#include "Light.h"
#include "Object.h"
#include "Source.h"
#include "Sphere.h"
#include "Plane.h"
#include "Triangle.h"

#include <omp.h>

#include "Primitives/Vertex.h"
#include "Primitives/ShapeGenerator.h"

// the system uses 12 threads by default
#define NUM_THREADS 12

// used for vector and matrix operations
using glm::vec3;
using glm::mat4;
using glm::vec4;
using glm::length;
using glm::pow;
using glm::floor;

vec3 eyePosition(3, 3.5, -7.0);  // position of the eye of the camera

// save a bitmap image
void saveBMP(
	const char* filename,
	int width,
	int height,
	int dpi,
	ColorArray& data) {

	// some old C style file handling code, update to modern C++ if possible
	FILE* file;
	int numPixels = width * height;  // total number of pixels in the image
	// 4 color channels for every pixel, (r, g, b, a)
	// every channel contains 1 byte of information
	int sizeBytes = 4 * numPixels;  // size occupied by image data in bytes
	int extra = 54;  // extra bytes required for metadata

	int totalFileSize = extra + sizeBytes;

	double factor = 39.375;
	int ppm = dpi * static_cast<int>(factor);

	// setup the buffers for storing the meta-info
	unsigned char bmpFileHeader[14] = { 'B', 'M', 0, 0, 0, 0,  0, 0, 0, 0,  (unsigned char)extra, 0, 0, 0 };
	unsigned char bmpInfoHeader[40] = { 40, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  1, 0, 24, 0 };

	// get the size of an individual element in the above buffers
	size_t unitSize = sizeof(unsigned char);

	// get the size of the 2 buffers declared above
	size_t fileHeaderSize = sizeof(bmpFileHeader) / sizeof(unsigned char);
	size_t infoHeaderSize = sizeof(bmpInfoHeader) / sizeof(unsigned char);

	assert(fileHeaderSize + infoHeaderSize == extra);

	// store the total file size, image dimensions, etc into the headers.
	{
		// store the integer values 8 bits at a time 4 times, right shift to move the next batch to the LSB
		bmpFileHeader[2] = (unsigned char)(totalFileSize);
		bmpFileHeader[3] = (unsigned char)(totalFileSize >> 8);
		bmpFileHeader[4] = (unsigned char)(totalFileSize >> 16);
		bmpFileHeader[5] = (unsigned char)(totalFileSize >> 24);

		bmpInfoHeader[4] = (unsigned char)(width);
		bmpInfoHeader[5] = (unsigned char)(width >> 8);
		bmpInfoHeader[6] = (unsigned char)(width >> 16);
		bmpInfoHeader[7] = (unsigned char)(width >> 24);

		bmpInfoHeader[8] = (unsigned char)(height);
		bmpInfoHeader[9] = (unsigned char)(height >> 8);
		bmpInfoHeader[10] = (unsigned char)(height >> 16);
		bmpInfoHeader[11] = (unsigned char)(height >> 24);

		bmpInfoHeader[21] = (unsigned char)(sizeBytes);
		bmpInfoHeader[22] = (unsigned char)(sizeBytes >> 8);
		bmpInfoHeader[23] = (unsigned char)(sizeBytes >> 16);
		bmpInfoHeader[24] = (unsigned char)(sizeBytes >> 24);

		bmpInfoHeader[25] = (unsigned char)(ppm);
		bmpInfoHeader[26] = (unsigned char)(ppm >> 8);
		bmpInfoHeader[27] = (unsigned char)(ppm >> 16);
		bmpInfoHeader[28] = (unsigned char)(ppm >> 24);

		bmpInfoHeader[29] = (unsigned char)(ppm);
		bmpInfoHeader[30] = (unsigned char)(ppm >> 8);
		bmpInfoHeader[31] = (unsigned char)(ppm >> 16);
		bmpInfoHeader[32] = (unsigned char)(ppm >> 24);
	}

	// buffers setup done, save data to the file now
	file = fopen(filename, "wb");

	if (file) {

		fwrite(bmpFileHeader, unitSize, fileHeaderSize, file);
		fwrite(bmpInfoHeader, unitSize, infoHeaderSize, file);

		// loop over every pixel color and save it in the file
		for (int i = 0; i < numPixels; ++i) {

			unsigned char color[3];

			for (int c = 0; c < 3; ++c)
				color[c] = (unsigned char)(data.colors[i][c] * 255);

			fwrite(color, unitSize, 3, file);
		}

		fclose(file);
	}
}

Ray constructRay(Camera & camera, double moveHorizontal, double moveVertical) {

	// get the camera co-ordinate frame
	vec3 viewDirection = camera.getViewDirection();
	vec3 right = camera.getRight();
	vec3 down = camera.getDown();

	vec3 rayDirection = viewDirection + ((float)moveHorizontal - 0.5f) * right + ((float)moveVertical - 0.5f) * down;

	Ray ray(camera.getPosition(), rayDirection);

	return ray;
}

// (i, j) is the pixel
void positionCameraForPixel(int i, int j, int width, int height,
	double& moveHorizontal, double& moveVertical,
	double add
) {

	double aspect = (double)width / (double)height;

	if (aspect > 1) {
		// the image is wider than it is tall
		moveHorizontal = ((i + add) / width) * aspect - ((((double)width - (double)height) / (double)height) / 2);
		moveVertical = (((double)height - (double)j) + add) / height;
	}
	else if (aspect < 1) {
		// the imager is taller than it is wide
		moveHorizontal = (i + add) / width;
		moveVertical = ((((double)height - (double)j) + add) / height) / aspect - ((((double)height - (double)width) / (double)width) / 2);
	}
	else {
		// the image is square
		moveHorizontal = (i + 0.5) / width;
		moveVertical = (((double)height - (double)j) + add) / height;
	}

}

// get the angle by which to rotate the ray
// index is the refractive index of the object
float getAngle(vec3 v1, vec3 v2, float index) {

	// get the angle of incidence
	float incident = glm::acos(dot(normalize(v1), normalize(v2)));

	// now get the angle of refraction
	float refraction = glm::asin(glm::sin(incident) / index);

	return incident - refraction;
}

// intersection is the intersection point on the object
void getRefractedRay(Sphere *sphere, vec3 intersection, Ray &ray) {

	// get the normal at this intersection point
	vec3 normal = sphere->getNormalAt(intersection);

	// get the angle by which the ray needs to rotate
	float angle = getAngle(normal, ray.getOrigin() - intersection, sphere->getIndex());

	vec3 axis = cross(-normal, intersection - ray.getOrigin());
	ray.rotate(-angle, intersection, axis);

	double closest = sphere->findIntersection(ray);

	// get the new intersection point
	vec3 newIntersection = ray.getOrigin() + (float)closest * ray.getDirection();
	// get the new normal at this new point
	vec3 newNormal = sphere->getNormalAt(newIntersection);

	ray.rotate(angle, newIntersection, axis);
}

// find the closest object that the ray shot from the camera intersects with
void findClosestObject(Ray &ray, std::vector<Object*> & objects, double& closest, int& closestIndex, 
					   bool refract, bool shadowed = false) {

	closest = std::numeric_limits<double>::max();
	closestIndex = -1; // index of the object closest to ray

	for (int k = 0; k < objects.size(); ++k) {

		Object* object = objects[k];

		double intersection = object->findIntersection(ray);

		if (object->getType() == Type::PLANE && intersection > 30)
			continue;

		// update the index of the closest object if this object is closer
		if (intersection > 0 && intersection < closest) {
			closest = intersection;
			closestIndex = k;
		}
	}

	if (closestIndex != -1 && refract) {

		// save closest and closestIndex
		double oldclosest = closest;
		int oldindex = closestIndex;

		// is the closest object transparent?
		Object* object = objects[closestIndex];

		// check if there is an intersection with a transparent object
		if (object->getType() == Type::SPHERE && ((Sphere*)object)->isTransparent()) {
			// refract baby
			// compute a new refracted ray coming out of the sphere and check if that intersects with 
			// any other object in the scene or not
			getRefractedRay((Sphere*)object, ray.getOrigin() + (float)closest * ray.getDirection(), ray);
			// the ray argument will be modified by the above call
			findClosestObject(ray, objects, closest, closestIndex, refract);

			if (closestIndex == -1 && !shadowed) {
				closest = oldclosest;
				closestIndex = oldindex;
			}

		}
	}
}

// given the ray from the camera/object and the normal to the intersection point
// compute the reflection
vec3 getReflectedVector(vec3 objectNormal, vec3 cameraRayDirection) {

	double dot1 = dot(objectNormal, -cameraRayDirection);
	vec3 scalar1 = (float)dot1 * objectNormal;
	vec3 add1 = scalar1 + cameraRayDirection;
	vec3 scalar2 = 2.0f * add1;
	vec3 add2 = scalar2 - cameraRayDirection;
	vec3 reflectionDirection = normalize(add2);

	return reflectionDirection;
}

// compute the color value at the intersection point by taking into
// account all the light sources in the scene
Color getColorAt(vec3 intersection, vec3 cameraRayDirection, Object * object,
				 std::vector<Source*> & lights, std::vector<Object*> & scene,
				 double accuracy, double ambient, int depth) {

	if (depth == 0)
		return Color();

	if (object->getType() == Type::PLANE) {
		// move the point over the plane a little bit to prevent any 
		// floating point mistakes in computing shadows
		intersection = intersection + vec3(0, 0.01, 0);
	}
	else if (object->getType() == Type::TRIANGLE) {
		intersection = intersection - vec3(0.01, 0.01, 0.01);
	}

	// get the normal to the object at this intersection point
	vec3 objectNormal = object->getNormalAt(intersection);
	// get the color of the object as well
	Color objectColor = object->getColor();

	// check if the object needs a checkered pattern (for the plane)
	if (objectColor.a == 2) {

		// determine the color at this intersection point
		int square = (int)(floor(intersection.x) + floor(intersection.z));  // since the plane lies on the x-z plane

		if (square % 2 == 0) {
			// black tile
			objectColor.g = objectColor.b = 0;
			objectColor.r = 1;
		}
		else {
			// white tile
			objectColor.r = objectColor.g = 1;
			objectColor.b = 0;
		}
	}

	// final color of the pixel, will be updated if needed in the following code
	Color finalColor = ambient * objectColor;

	// now we will perform some reflections calculations
	// first check if shininess is on or not
	vec3 reflectionDirection = getReflectedVector(objectNormal, cameraRayDirection);

	if (objectColor.a > 0 && objectColor.a <= 1) {

		// construct a ray for this reflection
		Ray reflectionRay(intersection, reflectionDirection);
		// determine what this ray intersects with
		double closest;
		int closestIndex;

		findClosestObject(reflectionRay, scene, closest, closestIndex, true);

		if (closestIndex != -1 && closest > accuracy) {

			// get the new intersection point
			vec3 newIntersection = reflectionRay.getOrigin() + (float)closest * reflectionRay.getDirection();

			// let the recursive ray-tracing begin
			finalColor = finalColor + objectColor.a * getColorAt(newIntersection, reflectionDirection,
				scene[closestIndex], lights, scene, accuracy, ambient, --depth);
		}
	}

	// apply all the light sources in the coloring
	for (size_t i = 0; i < lights.size(); ++i) {

		Source* light = lights[i];

		// get the position of the light
		vec3 lightPosition = light->getPosition();

		// get the shadow ray direction from the intersection point to the light source
		vec3 shadowRayDirection = lightPosition - intersection;
		// compute the intensity of light falling on the point
		double intensity = dot(normalize(shadowRayDirection), objectNormal);

		if (intensity > 0) {

			// get the magnitude of the ray direction
			double rayLength = length(shadowRayDirection);

			// create a ray that originates at the intersection point and hits the light source
			Ray shadowRay(intersection, shadowRayDirection);
			// check if this ray intersects with any other object in the scene
			double closest;
			int closestIndex;

			findClosestObject(shadowRay, scene, closest, closestIndex, true, true);

			bool shadowed = false;  // turn off the shadows initially

			if (closestIndex != -1 && closest > accuracy && closest < rayLength) {
				// the shadow ray is obstructed
				shadowed = true;
			}

			// technically, the 2 if conditions here can be merged into one, but I think
			// this setting looks more intuitive. pro tip: when negating the above condition use demorgan's law :)
			if (!shadowed) {

				finalColor = finalColor + intensity * light->getColor() * objectColor;

				// do some specular lighting if shininess is set for the object
				if (objectColor.a > 0 && objectColor.a <= 1) {
					// shininess parameter set

					double specular = dot(reflectionDirection, normalize(shadowRayDirection));

					if (specular > 0) {
						// set the tightness of the reflected light
						specular = pow(specular, 10);
						finalColor = finalColor + specular * objectColor.a * light->getColor();
					}
				}
			}
		}

	}

	return finalColor.clip();
}

int main() {

	int dpi = 72;
	int width = 12 * 70;  // high res image
	int height = 9 * 70;

	const int aaDepth = 3;  // anti-aliasing depth
	// double threshold = 0.1;     // threshold for anti-aliasing

	// setup the timers to measure the time taken for rendering
	clock_t t1, t2;
	t1 = clock();

	double ambient = 0.4;
	double accuracy = 0.01;

	ColorArray arr(width, height);

	// origin
	vec3 origin(0, 0, 0);

	// setup the camera
	vec3 center = origin;  // look at the origin
	Camera camera(eyePosition, center);  // sets up the co-oridnate frame for the camera

	// setup the light for the scene
	vec3 lightPosition(-1, 4, -5);  // -7, 10, -10
	Color lightColor(1.0, 1.0, 1.0, 0.0);  // white light

	Light sceneLight(lightPosition, lightColor);
	Light light2(vec3(3, 1.5, 4), Color(1.0, 1.0, 1.0, 0.0));

	// setup the objects for the scene
	// starting with the sphere
	vec3 spherePosition = origin + vec3(-3.0, -0.4, -2.0);
	double sphereRadius = 1.5;
	Color sphereColor(0.1, 0.8, 0.1, 0.7);  // green

	Sphere sphere(sphereRadius, spherePosition, sphereColor);
	Sphere sphere1(0.5, origin + vec3(1.2, -1.5, -3), Color(0.4, 0.0, 0.3, 0.2));
	Sphere sphere2(1.4, origin + vec3(3.8, -0.2, 0), Color(0.3, 0.3, 1.0, 0.2));

	vec3 a(1.0, 3.0, -1.5);
	vec3 b(-0.9, -1.8, -1.0);
	vec3 c(1.9, -1.8, -1.0);
	vec3 normal = normalize(cross(c - a, b - a));
	Triangle triangle(a, b, c, normal, normal, normal, Color(0.1, 0.2, 0.6, 0.5));
	normal = -normal;
	Triangle tr2(a, b, c, normal, normal, normal, Color(0.1, 0.2, 0.6, 0.5));

	// now the plane
	Plane plane(vec3(0, 1, 0), -2.0, Color(0.5, 0.25, 0.25, 2.0));

	// create a vector holding all objects in the scene
	std::vector<Object*> sceneObjects = { &sphere, &plane, &sphere1, &sphere2 };

	// make a teapot as well
	ShapeData* teapot = ShapeGenerator::makeTeapot(40);
	// get all the triangles in the mesh that make up the teapot and add all of them to the sceneObjects
	Vertex *vertices = teapot->vertexData();
	unsigned short* indices = teapot->indexData();

	float minX, minY, minZ, maxX, maxY, maxZ;
	minX = minY = minZ = std::numeric_limits<float>::max();
	maxX = maxY = maxZ = std::numeric_limits<float>::min();

	mat4 T = glm::translate(vec3(0.0f, -1.0f, 0.0f)) *
			 glm::rotate(glm::radians(70.0f), vec3(0, 1, 0)) * 
			 glm::rotate(glm::radians(-90.0f), vec3(1, 0, 0)) * 
			 glm::scale(vec3(0.8f, 0.8f, 0.8f));

	for (size_t i = 0; i < teapot->getNumIndices(); i = i + 3) {

		Vertex v = vertices[indices[i]];
		// extract the position and normal to the vertex
		vec3 a = v.position;
		vec3 normA = v.normal;

		// do the same for vertices b and c as well
		v = vertices[indices[i + 1]];
		vec3 b = v.position;
		vec3 normB = v.normal;

		v = vertices[indices[i + 2]];
		vec3 c = v.position;
		vec3 normC = v.normal;

		// model to world transformations for the vertices and the normals

		a = T * vec4(a, 1.0);
		b = T * vec4(b, 1.0);
		c = T * vec4(c, 1.0);

		normA = mat3(glm::transpose(glm::inverse(T))) * normA;
		normB = mat3(glm::transpose(glm::inverse(T))) * normB;
		normC = mat3(glm::transpose(glm::inverse(T))) * normC;

		// make the triangle
		Triangle *triangle = new Triangle(a, b, c, normA, normB, normC, Color(0.2, 0.5, 0.3, 0.5));
		sceneObjects.push_back(triangle);
	}

	// create a vector holding all the light sources in the scene as well
	std::vector<Source*> sceneLights = { &sceneLight };

	double moveHorizontal, moveVertical;  // for adjusting the camera based on the aspect ratio

	for (int i = 0; i < width; ++i) {  // initial...complete

		for (int j = 0; j < height; ++j) {

			// get the pixel index in the color buffer
			// int pixel = CHUNK_SIZE * j + i - initial;
			int pixel = width * j + i;

			// Color& color = buffers[id]->colors[pixel];  // set the color value at the current pixel later
			Color& color = arr.colors[pixel];

			for (int x = 0; x < aaDepth; ++x) {
				for (int y = 0; y < aaDepth; ++y) {

					double add = 0.5;  // anti-aliasiang is off by default
					if (aaDepth > 1) {
						// anti aliasing enabled
						add = (double)x / (double)aaDepth;
					}

					positionCameraForPixel(i, j, width, height, moveHorizontal, moveVertical, add);

					// construct a ray that will go through the pixel into the screen now
					Ray ray = constructRay(camera, moveHorizontal, moveVertical);

					// find the intersections of this ray with each and every object in the scene
					// and choose the closest one
					int closestIndex;
					double closest;

					findClosestObject(ray, sceneObjects, closest, closestIndex, true);

					if (closestIndex != -1 && closest > accuracy) {
						// get the intersection point
						vec3 intersection = ray.getOrigin() + (float)closest * ray.getDirection();

						// set the recursion depth
						int depth = 5;

						// get the color at the intersection
						Color pointColor = getColorAt(intersection, ray.getDirection(), sceneObjects[closestIndex],
							sceneLights, sceneObjects, accuracy, ambient, depth);
						// color of the object at the intersection point

						// average out the color
						color = color + pointColor / ((double)aaDepth * (double)aaDepth);
					}

				}
			}

		}

	}

	saveBMP("image_alias.bmp", width, height, dpi, arr);

	// stop the clock
	t2 = clock();
	float difference = ((float)t2 - (float)t1) / 1000.0f;

	std::cout << "Time taken: " << difference << " seconds.\n";

	return 0;
}