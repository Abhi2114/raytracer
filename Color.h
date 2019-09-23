#pragma once
#include <stdexcept>

struct Color {

	double r, g, b, a;  // red, green, blue and alpha channels

	Color() { r = g = b = a = 0; /* black */ }

	Color(double r, double g, double b, double a) : r(r), g(g), b(b), a(a) {}

	double operator[] (int i) {
		if (i == 0)
			return b;
		else if (i == 1)
			return g;
		else if (i == 2)
			return r;
		else if (i == 3)
			return a;
		else
			throw std::out_of_range("Out of bounds mate...");
	}

	double getBrightness() {
		return (r + g + b) / 3;
	}

	// scalar * color, scale the color
	friend Color operator * (double scalar, Color color) {
		return Color(scalar * color.r, scalar * color.g, scalar * color.b, color.a);
	}

	// add 2 colors together
	Color operator + (Color color) {
		return Color(r + color.r, g + color.g, b + color.b, a);
	}

	// multiply 2 colors
	Color operator * (Color color) {
		return Color(r * color.r, g * color.g, b * color.b, a);
	}

	Color average(Color color) {
		return Color((r + color.r) / 2, (g + color.g) / 2, (b + color.b) / 2, a);
	}

	// divide the color value
	Color operator / (double scalar) {
		return Color(r / scalar, g / scalar, b / scalar, a);
	}

	Color clip() {

		double total = r + g + b;
		double excess = total - 3;
		if (excess > 0) {
			r = r + excess * (r / total);
			g = g + excess * (g / total);
			b = b + excess * (b / total);
		}
		if (r > 1) r = 1;
		if (g > 1) g = 1;
		if (b > 1) b = 1;

		if (r < 0) r = 0;
		if (g < 0) g = 0;
		if (b < 0) b = 0;

		return Color(r, g, b, a);
	}
};

// simple smart pointer to manage memory for the Color buffer
struct ColorArray {

	Color* colors;

	ColorArray(int width, int height) {
		colors = new Color[width * height];
	}

	~ColorArray() {
		delete[] colors;
		colors = nullptr;
	}
};