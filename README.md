LM-ICP

# Caution
This code is not sophisticated, impractical, and it may contain mistakes.
It is not beyond my personal demonstration.

# What's this
Experimental implementation of LM-ICP(Levenberg-Marquard Iterative Closest Point).
In sample code, it will perform camera transformation for 4 3D points into image plane,
with parameters which makes projected points close to 4 ground-truth 2D points.

The camera model is from Guosheng Hu et al.(2017)

# Result
![Result](https://i.imgur.com/jfVImFw.jpg "Result")
Green points are default 3D points which are projected with default camera parameters.
Blue circles denote "ground-truth" or "objective" points. LM-ICP aims to obtain specific camera parameters which project green points as close as possible to blue circles.
Red points are actual projection with parameters obtained in LM-ICP.

# Requirements
Eigen, OpenCV

# Reference
Robust Registration of 2D and 3D Point Sets
https://www.robots.ox.ac.uk/~cvrg/michaelmas2004/fitzgibbon01c.pdf

Efficient 3D Morphable Face Model Fitting
http://epubs.surrey.ac.uk/813634/1/2017_PR.pdf