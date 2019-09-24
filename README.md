# raytracer
Raytraced rendering of images with various shapes using realistic reflections and shadows.

Intersections supported for Planes, Spheres and Triangles. 

The recursion depth for computing reflections is currently set to 5, increase it at your own risk :)

The render times can be greatly improved just by using some OpenMP directives to make use of all the cores on your machine. Would 
really appreciate it if someone could help me port this to a CUDA runtime to make use of the GPU, since making CUDA work with C++ classes 
is a real pain.

Following are some of the results I was able to produce with this implementation:

![](images/ellipsoids.PNG)
![](images/image.bmp)
![](images/ray_tracer.PNG)
![](images/refraction1.PNG)
![](images/teapot.PNG)
![](images/triangle.PNG)
![](images/two_spheres_refraction2.PNG)