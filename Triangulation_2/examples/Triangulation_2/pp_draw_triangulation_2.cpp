#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/IO/write_VTU.h>
#include <CGAL/IO/io.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/draw_voronoi_diagram_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K>      Triangulation;
typedef Triangulation::Point                                Point;
// typedefs to get a Site_2 type (used by some CGAL adaptors)
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<Triangulation>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<Triangulation> AP;
typedef CGAL::Voronoi_diagram_2<Triangulation,AT,AP>                                    VD;

// typedef for the result type of the point location
typedef AT::Site_2                    Site_2;

int main(int argc, char* argv[]) {
  std::ifstream in((argc>1)?argv[1]:"data/triangulation_prog1_copy.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;

  // Read points into a vector so we can compute the convex hull
  std::vector<Point> points;
  std::copy(begin, end, std::back_inserter(points));

    Triangulation t;
  if(!points.empty()) t.insert(points.begin(), points.end());

  VD vd;
  // Convert points to Site_2 (example)
  std::vector<Site_2> sites;
  sites.reserve(points.size());
  for(const Point& p : points) {
    // construct Site_2 from Point (Site_2 may be a Point_2 type)
    Site_2 site(p);
    sites.emplace_back(site);
    vd.insert(site);
  }
  // print the sites inserted
  for(const Site_2& site : sites) {
    std::cout << "Inserted site: " << site << std::endl;
  }
  assert( vd.is_valid() );

  // Compute convex hull of the input points and insert constraints along the hull
  if(points.size() >= 2) {
    std::vector<Point> hull;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(hull));
    if(hull.size() >= 2) {
      std::cout << "Inserting " << hull.size() << " constraint edges along convex hull" << std::endl;
      for(std::size_t i = 0, n = hull.size(); i < n; ++i) {
        const Point& p = hull[i];
        const Point& q = hull[(i+1) % n];
        t.insert_constraint(p, q);
      }
    }
  }

  std::cout << "=== Triangulation Properties ===" << std::endl;
  std::cout << "Number of vertices: " << t.number_of_vertices() << std::endl;
  std::cout << "Number of faces: " << t.number_of_faces() << std::endl;
  
  // Count edges
  int num_edges = 0;
  for(auto eit = t.finite_edges_begin(); eit != t.finite_edges_end(); ++eit) {
    num_edges++;
  }
  std::cout << "Number of edges: " << num_edges << std::endl;
  
  std::cout << "Dimension: " << t.dimension() << std::endl;
  std::cout << "Is valid: " << (t.is_valid() ? "true" : "false") << std::endl;
  
  std::cout << "\nPoints (vertices):" << std::endl;
  int index = 0;
  for(auto vit = t.finite_vertices_begin(); vit != t.finite_vertices_end(); ++vit) {
    std::cout << "  [" << index << "] " << vit->point() << std::endl;
    index++;
  }

  std::cout << "\nConnections (edges):" << std::endl;
  for(auto eit = t.finite_edges_begin(); eit != t.finite_edges_end(); ++eit) {
    auto face = eit->first;
    int index = eit->second;
    auto v1 = face->vertex((index + 1) % 3);
    auto v2 = face->vertex((index + 2) % 3);
    std::cout << "  " << v1->point() << " <--> " << v2->point() << std::endl;
  }

  std::cout << "\nAdjacency list (vertex neighbors):" << std::endl;
  for(auto vit = t.finite_vertices_begin(); vit != t.finite_vertices_end(); ++vit) {
    std::cout << "  " << vit->point() << " connects to: ";
    auto vc = t.incident_vertices(vit);
    auto done = vc;
    if(vc != 0) {
      do {
        if(!t.is_infinite(vc)) {
          std::cout << vc->point() << ", ";
        }
        ++vc;
      } while(vc != done);
    }
    std::cout << std::endl;
  }
  


  // Write VTU using CGAL's Mesh_2 writer (ASCII mode)
  {
    std::ofstream ofs((argc>2)? argv[2] : "triangulation.vtu");
    if(ofs) {
      CGAL::IO::write_VTU(ofs, t, CGAL::IO::ASCII);
    } else {
      std::cerr << "Could not open output file for VTU." << std::endl;
    }
  }

  CGAL::draw(t);
  CGAL::draw(vd);
  return EXIT_SUCCESS;
}
