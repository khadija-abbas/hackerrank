#ifndef VORNOI_DIAGRAM_H
#define VORNOI_DIAGRAM_H

#include <vector>
#include <string>
#include <limits>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "QuadTree.h"
#include <queue>
#include <set>
#include <map>

// Constants
#define PI 3.14159265358979323846
#define EARTH_RADIUS 6371.0 // Earth radius in kilometers
#define EPSILON 1e-9        // Small value for floating-point comparisons

// Edge structure for Voronoi diagram
struct Edge
{
    Point start; // Start point of the edge
    Point end;   // End point of the edge
    Point site1; // First site
    Point site2; // Second site

    Edge() {}
    Edge(const Point &s1, const Point &s2) : site1(s1), site2(s2) {}
};

// Region structure for Voronoi diagram
struct Region
{
    Point site;                  // Site point
    std::vector<Edge> edges;     // Edges of the region
    std::vector<Point> vertices; // Vertices of the region
};

// Point with additional data for Fortune's algorithm
struct FortunePoint
{
    double x, y;    // Projected coordinates
    Point original; // Original lat/lon point
    int siteIndex;  // Index of the site point

    FortunePoint() : x(0), y(0), siteIndex(-1) {}

    FortunePoint(double x, double y, const Point &orig, int idx)
        : x(x), y(y), original(orig), siteIndex(idx) {}
};

// Event types for Fortune's algorithm
enum EventType
{
    SITE_EVENT,
    CIRCLE_EVENT
};

// Event for Fortune's algorithm
struct Event
{
    EventType type;
    FortunePoint point;
    double y;   // Y-coordinate of the event (sweepline position)
    void *data; // Additional data (e.g., pointer to arc for circle events)

    Event(EventType t, const FortunePoint &p, double eventY, void *d = nullptr)
        : type(t), point(p), y(eventY), data(d) {}

    // Comparison operator for priority queue (sorts by y-coordinate in descending order)
    bool operator<(const Event &other) const
    {
        return y < other.y; // Note: This will make the priority queue a max-heap
    }
};

// Arc in the beachline
struct Arc
{
    FortunePoint site;  // Site that defines this arc
    Arc *prev;          // Previous arc in the beachline
    Arc *next;          // Next arc in the beachline
    Event *circleEvent; // Circle event that this arc is part of, if any

    Arc(const FortunePoint &s) : site(s), prev(nullptr), next(nullptr), circleEvent(nullptr) {}
};

// Half-edge data structure for Voronoi diagram
struct HalfEdge
{
    FortunePoint start; // Start point of the half-edge
    FortunePoint end;   // End point of the half-edge (may be null for infinite edges)
    int site1, site2;   // Indices of the two sites separated by this edge
    HalfEdge *twin;     // Twin half-edge

    HalfEdge(int s1, int s2) : site1(s1), site2(s2), twin(nullptr) {}
};

// Calculate Haversine distance between two points (in km)
double haversineDistanceVernoi(const Point &p1, const Point &p2);

// Main function to construct Voronoi diagram using Fortune's algorithm
std::vector<Region> fortunesVoronoiDiagram(std::vector<Point> &sites, double minLat, double maxLat, double minLon, double maxLon);

// Wrapper function for compatibility
std::vector<Region> constructVoronoiDiagram(std::vector<Point> &sites, double minLat, double maxLat, double minLon, double maxLon);

// Find the nearest site to a query point using the Voronoi diagram
Point findNearestSiteInVoronoi(const std::vector<Region> &regions, double queryLat, double queryLon);

// Check if a point is inside a polygon using ray casting algorithm
bool pointInPolygon(const Point &point, const std::vector<Point> &vertices);

#endif // VORNOI_DIAGRAM_H