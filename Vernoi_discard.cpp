// this is using brutes force algorithm but for testcases we are keeping this file for now to doublecheck

#include "Vernoi_discard.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

// Calculate Haversine distance between two points (in km)
double haversineDistanceVernoi(const Point &p1, const Point &p2)
{
    const double PI = 3.14159265358979323846;

    double lat1Rad = p1.lat * PI / 180.0;
    double lon1Rad = p1.lon * PI / 180.0;
    double lat2Rad = p2.lat * PI / 180.0;
    double lon2Rad = p2.lon * PI / 180.0;

    double dLat = lat2Rad - lat1Rad;
    double dLon = lon2Rad - lon1Rad;
    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(lat1Rad) * cos(lat2Rad) *
                   sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    double earthRadius = 6371.0; // Earth's radius in kilometers

    return earthRadius * c;
}

// Generate a simplified Voronoi diagram -Brute force algorithm
// Instead of computing exact Voronoi cells, we'll create approximate regions
std::vector<Region> constructVoronoiDiagram(std::vector<Point> &sites, double minLat, double maxLat, double minLon, double maxLon)
{
    std::cout << "Starting Voronoi diagram construction..." << std::endl;

    int n = sites.size();
    std::vector<Region> regions(n);

    // Initialize regions with their sites
    for (int i = 0; i < n; ++i)
    {
        regions[i].site = sites[i];
        // std::cout << "Initialized region " << i + 1 << " of " << n << std::endl;
    }

    // For each site, create a simplified polygonal region
    for (int i = 0; i < n; ++i)
    {
        // std::cout << "Processing region " << i + 1 << " of " << n << std::endl;

        std::vector<Point> approximateVertices;
        const int numDirections = 12; // Number of rays to cast out from the site

        for (int angle = 0; angle < 360; angle += 360 / numDirections)
        {
            const double PI = 3.14159265358979323846;
            double radians = angle * PI / 180.0;
            double maxDist = haversineDistanceVernoi(Point(minLat, minLon), Point(maxLat, maxLon)); // Max possible distance

            // Scale down max distance to avoid numerical issues
            maxDist *= 0.5;

            // Start point is the site
            Point startPoint = sites[i];

            // End point is far away in the direction of angle
            double endLat = startPoint.lat + maxDist * sin(radians) / 111.0; // Approx 111km per degree of latitude
            double endLon = startPoint.lon + maxDist * cos(radians) / (111.0 * cos(startPoint.lat * PI / 180.0));

            // Clamp to boundaries
            endLat = std::max(minLat, std::min(maxLat, endLat));
            endLon = std::max(minLon, std::min(maxLon, endLon));

            // Binary search to find the boundary point
            double low = 0.0;
            double high = 1.0;
            Point boundaryPoint;

            for (int iter = 0; iter < 10; ++iter) // 10 iterations should be enough for good precision
            {
                double mid = (low + high) / 2.0;

                double testLat = startPoint.lat + mid * (endLat - startPoint.lat);
                double testLon = startPoint.lon + mid * (endLon - startPoint.lon);
                Point testPoint(testLat, testLon);

                // Find the closest site to this test point
                int closestSite = i;
                double minDist = haversineDistanceVernoi(testPoint, sites[i]);

                for (int j = 0; j < n; ++j)
                {
                    if (j != i)
                    {
                        double dist = haversineDistanceVernoi(testPoint, sites[j]);
                        if (dist < minDist)
                        {
                            minDist = dist;
                            closestSite = j;
                            break; // Optimization: we just need to know if any other site is closer
                        }
                    }
                }

                if (closestSite == i)
                {
                    // This point is still in our region, move further out
                    low = mid;
                    boundaryPoint = testPoint;
                }
                else
                {
                    // This point is in another region, move closer in
                    high = mid;
                }
            }

            // Add the boundary point to our vertices
            if (low > 0) // Only add if we found a valid boundary
            {
                approximateVertices.push_back(boundaryPoint);
            }
        }

        // Sort vertices in clockwise order around the site
        if (!approximateVertices.empty())
        {
            sort(approximateVertices.begin(), approximateVertices.end(), [&sites, i](const Point &a, const Point &b)
                 {
                    double angleA = atan2(a.lat - sites[i].lat, a.lon - sites[i].lon);
                    double angleB = atan2(b.lat - sites[i].lat, b.lon - sites[i].lon);
                    return angleA < angleB; });

            // Remove duplicate points
            approximateVertices.erase(
                unique(approximateVertices.begin(), approximateVertices.end()),
                approximateVertices.end());

            // std::cout << "Region " << i + 1 << " has " << approximateVertices.size() << " vertices" << std::endl;
        }

        // Add vertices to the region
        regions[i].vertices = approximateVertices;

        // Create edges between consecutive vertices
        if (!approximateVertices.empty())
        {
            for (size_t j = 0; j < approximateVertices.size(); ++j)
            {
                // Use the current site for both endpoints of the edge
                Edge edge(sites[i], sites[i]);
                edge.start = approximateVertices[j];
                edge.end = approximateVertices[(j + 1) % approximateVertices.size()];
                regions[i].edges.push_back(edge);
            }
        }
    }

    std::cout << "Voronoi diagram construction completed successfully." << std::endl;
    return regions;
}

// // Wrapper for backward compatibility
// std::vector<Region> constructVoronoiDiagram(std::vector<Point> &sites, double minLat, double maxLat, double minLon, double maxLon)
// {
//     return simplifiedVoronoiDiagram(sites, minLat, maxLat, minLon, maxLon);
// }

Point findNearestSiteInVoronoi(const std::vector<Region> &regions, double queryLat, double queryLon)
{
    Point queryPoint(queryLat, queryLon);

    // Use the Voronoi diagram to find which region contains the query point
    for (const auto &region : regions)
    {
        if (pointInPolygon(queryPoint, region.vertices))
        {
            return region.site;
        }
    }

    // Fallback to linear search if point is outside all polygons
    // std::std::cout << "linear" << std::std::endl;
    double minDist = std::numeric_limits<double>::max();
    Point nearest;

    for (const auto &region : regions)
    {
        double dist = haversineDistanceVernoi(queryPoint, region.site);
        if (dist < minDist)
        {
            minDist = dist;
            nearest = region.site;
        }
    }

    return nearest;
}

// Check if a point is inside a polygon using ray casting algorithm
bool pointInPolygon(const Point &point, const std::vector<Point> &vertices)
{
    if (vertices.size() < 3)
        return false;

    bool inside = false;
    size_t j = vertices.size() - 1;

    for (size_t i = 0; i < vertices.size(); i++)
    {
        if ((vertices[i].lat > point.lat) != (vertices[j].lat > point.lat) &&
            (point.lon < (vertices[j].lon - vertices[i].lon) * (point.lat - vertices[i].lat) /
                                 (vertices[j].lat - vertices[i].lat) +
                             vertices[i].lon))
        {
            inside = !inside;
        }
        j = i;
    }

    return inside;
}