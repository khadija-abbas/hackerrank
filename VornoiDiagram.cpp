#include "VornoiDiagram.h"

// Calculate Haversine distance between two points (in km)
double haversineDistanceVernoi(const Point &p1, const Point &p2)
{
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

    return EARTH_RADIUS * c;
}

// Project lat/lon coordinates to Cartesian for easier computation
std::vector<FortunePoint> projectToCartesian(const std::vector<Point> &sites)
{
    std::vector<FortunePoint> projectedPoints;
    projectedPoints.reserve(sites.size());

    // Use simple equirectangular projection centered on the mean of the points
    double meanLat = 0, meanLon = 0;
    for (const auto &site : sites)
    {
        meanLat += site.lat;
        meanLon += site.lon;
    }
    meanLat /= sites.size();
    meanLon /= sites.size();

    double cosLat = cos(meanLat * PI / 180.0);

    for (size_t i = 0; i < sites.size(); i++)
    {
        const Point &site = sites[i];
        // Scale longitude differences by cosine of latitude to account for spherical distortion
        double x = (site.lon - meanLon) * cosLat;
        double y = site.lat - meanLat;

        // Scale up for more precision in calculations
        x *= 1000.0;
        y *= 1000.0;

        projectedPoints.emplace_back(x, y, site, i);
    }

    return projectedPoints;
}

// Convert Cartesian coordinates back to lat/lon
Point unprojectToLatLon(const FortunePoint &point, double meanLat, double meanLon, double cosLat)
{
    // Undo the scaling
    double x = point.x / 1000.0;
    double y = point.y / 1000.0;

    double lon = meanLon + x / cosLat;
    double lat = meanLat + y;

    return Point(lat, lon);
}

// Calculate the intersection of two parabolas
FortunePoint getParabolaIntersection(const FortunePoint &f1, const FortunePoint &f2, double sweepline)
{
    // Handle edge cases
    if (fabs(f1.y - f2.y) < EPSILON)
    {
        return FortunePoint((f1.x + f2.x) / 2, (f1.y + f2.y) / 2, Point(0, 0), -1);
    }

    // Focus points (sites) of the two parabolas
    double x1 = f1.x, y1 = f1.y;
    double x2 = f2.x, y2 = f2.y;

    // Calculate intersection of two parabolas with directrix at sweepline
    double a1 = 1 / (2 * (y1 - sweepline));
    double a2 = 1 / (2 * (y2 - sweepline));

    double b1 = -2 * x1 * a1;
    double b2 = -2 * x2 * a2;

    double c1 = a1 * x1 * x1 + (y1 + sweepline) / 2;
    double c2 = a2 * x2 * x2 + (y2 + sweepline) / 2;

    // Solve quadratic equation for intersection
    double A = a1 - a2;
    double B = b1 - b2;
    double C = c1 - c2;

    double disc = B * B - 4 * A * C;
    if (disc < 0)
        disc = 0; // Handle numerical issues

    double sqrtDisc = sqrt(disc);
    double x_1 = (-B + sqrtDisc) / (2 * A);
    double x_2 = (-B - sqrtDisc) / (2 * A);

    // Choose the correct intersection point
    double x, y;
    if (y1 > y2)
    {
        x = (x1 <= x2) ? x_1 : x_2;
    }
    else
    {
        x = (x1 > x2) ? x_1 : x_2;
    }

    // Calculate y from the parabola equation
    y = a1 * x * x + b1 * x + c1;

    return FortunePoint(x, y, Point(0, 0), -1);
}

// Check if three consecutive arcs form a circle event
void checkCircleEvent(Arc *arc, double sweepline, std::priority_queue<Event> &eventQueue)
{
    if (!arc || !arc->prev || !arc->next)
        return;
    if (arc->prev->site.siteIndex == arc->next->site.siteIndex)
        return;

    // Get the three sites for the arcs
    FortunePoint a = arc->prev->site;
    FortunePoint b = arc->site;
    FortunePoint c = arc->next->site;

    // Check if the arcs form a valid triple for a circle event
    // Skip if points are collinear
    double det = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    if (fabs(det) < EPSILON)
        return;

    // Calculate the center of the circle through the three points
    double A = b.x - a.x;
    double B = b.y - a.y;
    double C = c.x - a.x;
    double D = c.y - a.y;
    double E = A * (a.x + b.x) + B * (a.y + b.y);
    double F = C * (a.x + c.x) + D * (a.y + c.y);
    double G = 2 * (A * D - B * C);

    if (fabs(G) < EPSILON)
        return; // Avoid division by near-zero

    double centerX = (D * E - B * F) / G;
    double centerY = (A * F - C * E) / G;

    // Calculate the radius
    double radius = sqrt(pow(a.x - centerX, 2) + pow(a.y - centerY, 2));

    // Calculate the bottom point of the circle
    double bottomY = centerY - radius;

    // Only consider the circle event if it's below the sweep line
    if (bottomY < sweepline + EPSILON)
        return;

    // Create circle event
    FortunePoint circlePoint(centerX, bottomY, Point(0, 0), -1);

    // If this arc already has a circle event, remove it
    if (arc->circleEvent)
    {
        // In a real implementation, we'd remove from the priority queue
        // Here we'll just mark it as invalid
        arc->circleEvent = nullptr;
    }

    // Create a new circle event
    Event *event = new Event(CIRCLE_EVENT, circlePoint, bottomY, arc);
    arc->circleEvent = event;
    eventQueue.push(*event);
}

// Add a site to the beachline and handle events
void addSiteToBeachline(const FortunePoint &site, Arc *&beachline, double sweepline,
                        std::vector<HalfEdge *> &edges, std::priority_queue<Event> &eventQueue)
{
    // If beachline is empty, create first arc
    if (!beachline)
    {
        beachline = new Arc(site);
        return;
    }

    // Find the arc above the new site
    Arc *arc = beachline;
    while (arc->next && getParabolaIntersection(arc->site, arc->next->site, sweepline).x < site.x)
    {
        arc = arc->next;
    }

    // If this arc has a circle event, remove it
    if (arc->circleEvent)
    {
        // In a real implementation, we'd remove from the priority queue
        // Here we'll just mark it as invalid
        arc->circleEvent = nullptr;
    }

    // Split the arc
    Arc *newArc = new Arc(site);
    Arc *splitArc = new Arc(arc->site);

    // Connect arcs properly
    newArc->prev = arc;
    newArc->next = splitArc;
    splitArc->prev = newArc;
    splitArc->next = arc->next;

    if (arc->next)
    {
        arc->next->prev = splitArc;
    }
    arc->next = newArc;

    // Create new half-edges for the new breakpoints
    HalfEdge *edge1 = new HalfEdge(arc->site.siteIndex, site.siteIndex);
    HalfEdge *edge2 = new HalfEdge(site.siteIndex, arc->site.siteIndex);
    edge1->twin = edge2;
    edge2->twin = edge1;

    // Start point is where the new site intersects the existing arc's parabola
    FortunePoint intersection = getParabolaIntersection(arc->site, site, sweepline);
    edge1->start = intersection;
    edge2->start = intersection;

    edges.push_back(edge1);
    edges.push_back(edge2);

    // Check for new circle events
    checkCircleEvent(arc, sweepline, eventQueue);
    checkCircleEvent(splitArc, sweepline, eventQueue);
}

// Handle a circle event
void handleCircleEvent(Event *event, Arc *&beachline, double sweepline,
                       std::vector<HalfEdge *> &edges, std::priority_queue<Event> &eventQueue)
{
    Arc *arc = static_cast<Arc *>(event->data);

    if (!arc || !arc->prev || !arc->next)
    {
        return; // Invalid arc
    }

    // Get the arcs involved
    Arc *leftArc = arc->prev;
    Arc *rightArc = arc->next;

    // Create a new edge at the center of the circle
    HalfEdge *edge1 = new HalfEdge(leftArc->site.siteIndex, rightArc->site.siteIndex);
    HalfEdge *edge2 = new HalfEdge(rightArc->site.siteIndex, leftArc->site.siteIndex);
    edge1->twin = edge2;
    edge2->twin = edge1;

    // The center of the circle is the Voronoi vertex
    edge1->start = event->point;
    edge2->start = event->point;

    edges.push_back(edge1);
    edges.push_back(edge2);

    // Remove the disappearing arc
    leftArc->next = rightArc;
    rightArc->prev = leftArc;

    delete arc;

    // Check for new circle events
    checkCircleEvent(leftArc, sweepline, eventQueue);
    checkCircleEvent(rightArc, sweepline, eventQueue);
}

// Clean up the beachline
void cleanupBeachline(Arc *beachline)
{
    Arc *current = beachline;
    while (current)
    {
        Arc *next = current->next;
        delete current;
        current = next;
    }
}

// Main Fortune's algorithm implementation
std::vector<Region> fortunesVoronoiDiagram(std::vector<Point> &sites, double minLat, double maxLat, double minLon, double maxLon)
{
    std::cout << "Starting Fortune's algorithm for Voronoi diagram construction..." << std::endl;

    if (sites.empty())
    {
        return {};
    }

    // Step 1: Project lat/lon coordinates to Cartesian
    std::vector<FortunePoint> projectedSites = projectToCartesian(sites);

    // Calculate mean lat/lon for later unprojection
    double meanLat = 0, meanLon = 0;
    for (const auto &site : sites)
    {
        meanLat += site.lat;
        meanLon += site.lon;
    }
    meanLat /= sites.size();
    meanLon /= sites.size();
    double cosLat = cos(meanLat * PI / 180.0);

    // Step 2: Initialize data structures
    std::priority_queue<Event> eventQueue;
    Arc *beachline = nullptr;
    std::vector<HalfEdge *> edges;

    // Add all site events to the queue
    for (const auto &site : projectedSites)
    {
        eventQueue.push(Event(SITE_EVENT, site, site.y));
    }

    // Step 3: Process events
    while (!eventQueue.empty())
    {
        Event currentEvent = eventQueue.top();
        eventQueue.pop();

        double sweepline = currentEvent.y;

        if (currentEvent.type == SITE_EVENT)
        {
            addSiteToBeachline(currentEvent.point, beachline, sweepline, edges, eventQueue);
        }
        else
        { // CIRCLE_EVENT
            handleCircleEvent(&currentEvent, beachline, sweepline, edges, eventQueue);
        }
    }

    // Step 4: Convert edges to Voronoi regions
    std::vector<Region> regions(sites.size());

    // Initialize regions with their sites
    for (size_t i = 0; i < sites.size(); i++)
    {
        regions[i].site = sites[i];
    }

    // Add edges to regions
    for (HalfEdge *edge : edges)
    {
        if (edge->site1 >= 0 && edge->site1 < static_cast<int>(regions.size()))
        {
            // Convert edge points back to lat/lon
            Point start = unprojectToLatLon(edge->start, meanLat, meanLon, cosLat);

            // For infinite edges, create a point far in the appropriate direction
            Point end;
            if (edge->end.x == 0 && edge->end.y == 0)
            { // Uninitialized end point (infinite edge)
                // Calculate direction perpendicular to the line between the two sites
                FortunePoint &s1 = projectedSites[edge->site1];
                FortunePoint &s2 = projectedSites[edge->site2];
                double dx = s2.x - s1.x;
                double dy = s2.y - s1.y;

                // Rotate 90 degrees and normalize
                double length = sqrt(dx * dx + dy * dy);
                double nx = -dy / length;
                double ny = dx / length;

                // Create a point far in this direction
                double scale = 1000.0; // Large enough to reach the boundary
                FortunePoint farPoint(edge->start.x + nx * scale, edge->start.y + ny * scale, Point(0, 0), -1);
                end = unprojectToLatLon(farPoint, meanLat, meanLon, cosLat);

                // Clip to boundaries
                end.lat = std::max(minLat, std::min(maxLat, end.lat));
                end.lon = std::max(minLon, std::min(maxLon, end.lon));
            }
            else
            {
                end = unprojectToLatLon(edge->end, meanLat, meanLon, cosLat);
            }

            // Create edge and add to region
            Edge voronoiEdge(sites[edge->site1], sites[edge->site2]);
            voronoiEdge.start = start;
            voronoiEdge.end = end;
            regions[edge->site1].edges.push_back(voronoiEdge);

            // Add vertices to the region
            bool startFound = false;
            for (const auto &vertex : regions[edge->site1].vertices)
            {
                if (fabs(vertex.lat - start.lat) < EPSILON && fabs(vertex.lon - start.lon) < EPSILON)
                {
                    startFound = true;
                    break;
                }
            }
            if (!startFound)
            {
                regions[edge->site1].vertices.push_back(start);
            }

            bool endFound = false;
            for (const auto &vertex : regions[edge->site1].vertices)
            {
                if (fabs(vertex.lat - end.lat) < EPSILON && fabs(vertex.lon - end.lon) < EPSILON)
                {
                    endFound = true;
                    break;
                }
            }
            if (!endFound)
            {
                regions[edge->site1].vertices.push_back(end);
            }
        }
    }

    // Sort vertices in clockwise order for each region
    for (auto &region : regions)
    {
        if (!region.vertices.empty())
        {
            // Calculate center of the vertices
            double centerLat = 0, centerLon = 0;
            for (const auto &vertex : region.vertices)
            {
                centerLat += vertex.lat;
                centerLon += vertex.lon;
            }
            centerLat /= region.vertices.size();
            centerLon /= region.vertices.size();

            // Sort vertices clockwise around the center
            Point center(centerLat, centerLon);
            sort(region.vertices.begin(), region.vertices.end(), [&center](const Point &a, const Point &b)
                 {
                double angleA = atan2(a.lat - center.lat, a.lon - center.lon);
                double angleB = atan2(b.lat - center.lat, b.lon - center.lon);
                return angleA < angleB; });
        }
    }

    // Clean up memory
    cleanupBeachline(beachline);
    for (auto edge : edges)
    {
        delete edge;
    }

    std::cout << "Fortune's algorithm completed successfully." << std::endl;
    return regions;
}

// Wrapper function for compatibility
std::vector<Region> constructVoronoiDiagram(std::vector<Point> &sites, double minLat, double maxLat, double minLon, double maxLon)
{
    return fortunesVoronoiDiagram(sites, minLat, maxLat, minLon, maxLon);
}

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
    std::cout << "No service found." << std::endl;
    return;
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