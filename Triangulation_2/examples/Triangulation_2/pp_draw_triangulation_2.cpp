#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_2<K>                            Triangulation;
typedef Triangulation::Point                                Point;

int main(int argc, char* argv[]) {
  std::ifstream in((argc>1)?argv[1]:"data/triangulation_prog1_copy.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;

  Triangulation t;
  t.insert(begin, end);

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

  CGAL::draw(t);

  return EXIT_SUCCESS;
}
