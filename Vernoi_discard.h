// this is using brutes force algorithm but for testcases we are keeping this file for now to doublecheck
#ifndef VERNOITRY_H
#define VERNOITRY_H

#include <vector>
#include <string>
#include <limits>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "QuadTree.h"

// Constants
const double PI = 3.14159265358979323846;
const double EARTH_RADIUS = 6371.0; // Earth's radius in kilometers
const double EPSILON = 1e-9;        // Epsilon for double comparison

// struct Point defined in Quadtree.h

struct Edge
{
    Point start, end;
    Point site1, site2;

    Edge(const Point &s1, const Point &s2) : site1(s1), site2(s2) {}
};

struct Region
{
    Point site;
    std::vector<Edge> edges;
    std::vector<Point> vertices;
};

// Function declarations
double haversineDistanceVernoi(const Point &p1, const Point &p2);
std::vector<Region> constructVoronoiDiagram(std::vector<Point> &sites, double minLat, double maxLat, double minLon, double maxLon);
Point findNearestSiteInVoronoi(const std::vector<Region> &regions, double queryLat, double queryLon);
bool pointInPolygon(const Point &point, const std::vector<Point> &vertices);

#endif