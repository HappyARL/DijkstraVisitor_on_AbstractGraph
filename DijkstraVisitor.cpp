// 1232387
#include <iostream>
#include <set>
#include <unordered_map>
#include <vector>

const int kMax = 2009000999;

// структура ребра
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

// абстрактный граф
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

// реализация на списке смежности
template <typename Vertex, typename Weight,
          typename Edge = DefaultEdge<Vertex, Weight> >
class Graph : public AbstractGraph<Vertex, Weight, Edge> {
 public:
  std::unordered_map<Vertex, std::vector<std::pair<Vertex, Weight> > >
      list_neighbors;

  Graph(size_t vertices_num, const std::vector<Edge>& edges)
      : AbstractGraph<Vertex, Weight, Edge>(vertices_num, edges.size()) {
    for (const auto& edge : edges) {
      list_neighbors[edge.Source()].push_back(
          std::make_pair(edge.Destination(), edge.WeightValue()));
      list_neighbors[edge.Destination()].push_back(
          std::make_pair(edge.Source(), edge.WeightValue()));
    }
  }

  std::vector<std::pair<Vertex, Weight> > GetNeighbours(
      const Vertex& vertex) final {
    return list_neighbors[vertex];
  }
};

// реализация на матрице смежности
template <typename Vertex, typename Weight,
          typename Edge = DefaultEdge<Vertex, Weight> >
class GraphMatrix : public AbstractGraph<Vertex, Weight, Edge> {
 public:
  std::vector<std::vector<Weight> > matrix_neighbour;

  GraphMatrix(size_t vertices_num, const std::vector<Edge>& edges)
      : AbstractGraph<Vertex, Weight, Edge>(vertices_num, edges.size()) {
    for (size_t i = 0; i < vertices_num; ++i) {
      std::vector<Weight> tmp(vertices_num, 0);
      matrix_neighbour.push_back(tmp);
    }
    for (const auto& edge : edges) {
      matrix_neighbour[edge.Source()][edge.Destination()] = edge.WeightValue();
      matrix_neighbour[edge.Destination()][edge.Source()] = edge.WeightValue();
    }
  }

  std::vector<std::pair<Vertex, Weight> > GetNeighbours(
      const Vertex& vertex) final {
    std::vector<std::pair<Vertex, Weight> > result;
    for (size_t i = 0; i < this->GetVerticesNumber(); ++i) {
      if (matrix_neighbour[vertex][i] > 0) {
        result.push_back(std::make_pair(i, matrix_neighbour[vertex][i]));
      }
    }
    return result;
  }
  /*
    void Print() {
      for (size_t i = 0; i < this->GetVerticesNumber(); ++i) {
        for (size_t j = 0; j < this->GetVerticesNumber(); ++j) {
          std::cout << matrix_neighbour[i][j] << " ";
        }
        std::cout << '\n';
      }
    }
  */
};

// Алгоритм Дейкстры
template <typename Vertex, typename Weight,
          typename Edge = DefaultEdge<Vertex, Weight> >
class DijkstraVisitor {
 public:
  DijkstraVisitor(Vertex start, Graph<Vertex, Weight, Edge> main_graph)
      : start_point_(start), graph_(main_graph) {}

  std::vector<Vertex> Dijkstra() {
    std::vector<Vertex> answer;
    answer.resize(graph_.list_neighbors.size());
    for (size_t i = 0; i < graph_.list_neighbors.size(); ++i) {
      answer[i] = kMax;
    }
    answer[start_point_] = 0;
    std::set<std::pair<Vertex, Vertex> > distance;
    distance.insert(std::make_pair(answer[start_point_], start_point_));
    while (!distance.empty()) {
      int index = distance.begin()->second;
      distance.erase(distance.begin());
      for (size_t j = 0; j < graph_.list_neighbors[index].size(); ++j) {
        int sum = answer[index] + graph_.list_neighbors[index].at(j).second;
        if (sum < answer[graph_.list_neighbors[index].at(j).first]) {
          distance.erase(
              std::make_pair(answer[graph_.list_neighbors[index].at(j).first],
                             graph_.list_neighbors[index].at(j).first));
          answer[graph_.list_neighbors[index].at(j).first] = sum;
          distance.insert(
              std::make_pair(answer[graph_.list_neighbors[index].at(j).first],
                             graph_.list_neighbors[index].at(j).first));
        }
      }
    }
    return answer;
  }

 private:
  Graph<Vertex, Weight, Edge> graph_;
  Vertex start_point_;
};

int main() {
  int maps;
  std::cin >> maps;
  std::vector<std::vector<int> > total_answer;
  while (maps != 0) {
    size_t rooms;
    int list_neighbors;
    std::cin >> rooms >> list_neighbors;
    std::vector<DefaultEdge<int, int> > edges;
    while (list_neighbors != 0) {
      int room_out;
      int room_in;
      int weight;
      std::cin >> room_out >> room_in >> weight;
      edges.push_back(DefaultEdge<int, int>(room_out, room_in, weight));
      --list_neighbors;
    }
    Graph<int, int, DefaultEdge<int, int> > graph(rooms, edges);
    // GraphMatrix<int, int, DefaultEdge<int, int> > graph_matrix(rooms, edges);
    // graph_matrix.Print(); - для проверки работоспособности
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
