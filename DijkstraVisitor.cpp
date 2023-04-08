// 1232387
#include <iostream>
#include <set>
#include <unordered_map>
#include <vector>

const int kMax = 2009000999;

template <typename Vertex, typename Weight>
struct DefaultEdge {
  Vertex src;
  Vertex dst;
  Weight weight;
  DefaultEdge(const Vertex& vertex_a, const Vertex& vertex_b,
              const Weight& value)
      : src(vertex_a), dst(vertex_b), weight(value) {}
  const Vertex& Source() const { return src; }
  const Vertex& Destination() const { return dst; }
  const Weight& WeightValue() const { return weight; }
};

template <typename Vertex, typename Weight,
          typename Edge = DefaultEdge<Vertex, Weight> >
class AbstractGraph {
 public:
  using VertexType = Vertex;
  using EdgeType = Edge;
  explicit AbstractGraph(size_t vertices_num, size_t edges_num = 0)
      : vertices_number_(vertices_num), edges_number_(edges_num) {}
  size_t GetVerticesNumber() const { return vertices_number_; }
  size_t GetEdgesNumber() const { return edges_number_; }
  virtual std::vector<std::pair<Vertex, Weight> > GetNeighbours(
      const Vertex& vertex) = 0;

 protected:
  size_t vertices_number_ = 0;
  size_t edges_number_ = 0;
};

template <typename Vertex, typename Weight,
          typename Edge = DefaultEdge<Vertex, Weight> >
class Graph : public AbstractGraph<Vertex, Weight, Edge> {
 public:
  std::unordered_map<Vertex, std::vector<std::pair<Vertex, Weight> > >
      connections;
  Graph(size_t vertices_num, const std::vector<Edge>& edges)
      : AbstractGraph<Vertex, Weight, Edge>(vertices_num, edges.size()) {
    for (const auto& edge : edges) {
      connections[edge.Source()].push_back(
          std::make_pair(edge.Destination(), edge.WeightValue()));
      connections[edge.Destination()].push_back(
          std::make_pair(edge.Source(), edge.WeightValue()));
    }
  }
  std::vector<std::pair<Vertex, Weight> > GetNeighbours(
      const Vertex& vertex) final {
    return connections[vertex];
  }
};

template <typename Vertex, typename Weight,
          typename Edge = DefaultEdge<Vertex, Weight> >
class DijkstraVisitor {
 private:
  Graph<Vertex, Weight, Edge> graph_;
  Vertex start_point_;

 public:
  DijkstraVisitor(Vertex start, Graph<Vertex, Weight, Edge> main_graph)
      : start_point_(start), graph_(main_graph) {}
  std::vector<Vertex> Dijkstra() {
    std::vector<Vertex> answer;
    answer.resize(graph_.connections.size());
    for (size_t i = 0; i < graph_.connections.size(); ++i) {
      answer[i] = kMax;
    }
    answer[start_point_] = 0;
    std::set<std::pair<Vertex, Vertex> > distance;
    distance.insert(std::make_pair(answer[start_point_], start_point_));
    while (!distance.empty()) {
      int index = distance.begin()->second;
      distance.erase(distance.begin());
      for (size_t j = 0; j < graph_.connections[index].size(); ++j) {
        int sum = answer[index] + graph_.connections[index].at(j).second;
        if (sum < answer[graph_.connections[index].at(j).first]) {
          distance.erase(
              std::make_pair(answer[graph_.connections[index].at(j).first],
                             graph_.connections[index].at(j).first));
          answer[graph_.connections[index].at(j).first] = sum;
          distance.insert(
              std::make_pair(answer[graph_.connections[index].at(j).first],
                             graph_.connections[index].at(j).first));
        }
      }
    }
    return answer;
  }
};

int main() {
  int maps;
  std::cin >> maps;
  std::vector<std::vector<int> > total_answer;
  while (maps != 0) {
    size_t rooms;
    int connections;
    std::cin >> rooms >> connections;
    std::vector<DefaultEdge<int, int> > edges;
    while (connections != 0) {
      int room_out;
      int room_in;
      int weight;
      std::cin >> room_out >> room_in >> weight;
      edges.push_back(DefaultEdge<int, int>(room_out, room_in, weight));
      --connections;
    }
    Graph<int, int, DefaultEdge<int, int> > graph(rooms, edges);
    int start_point;
    std::cin >> start_point;
    std::cout << std::endl;
    DijkstraVisitor<int, int, DefaultEdge<int, int> > dijkstra(start_point,
                                                               graph);
    std::vector<int> tmp = dijkstra.Dijkstra();
    total_answer.push_back(tmp);
    --maps;
  }
  for (int i = 0; i < (int)total_answer.size(); ++i) {
    for (int j = 0; j < (int)total_answer[i].size(); ++j) {
      std::cout << total_answer[i][j] << " ";
    }
    std::cout << std::endl;
  }
  return 0;
}
