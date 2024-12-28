#include <iostream>
#include <cmath>

#include <fstream>
#include <sstream>
#include <string>

#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <limits>
#include <queue>
#include <stack>
#include <chrono>

// Структура для представления узла графа
struct Node {
    int id;
    double longitude;
    double latitude;
    std::vector<std::pair<Node*, double>> neighbors;

    Node(int id, double lon, double lat) : id(id), longitude(lon), latitude(lat) {}
};

// Класс для представления графа
class Graph {
public:
    std::vector<Node*> nodes;
    std::unordered_map<std::string, Node*> node_map;

    // Генерируем уникю ключ по координатам
    std::string generateKey(double lon, double lat) {
        return std::to_string(lon) + "," + std::to_string(lat);
    }

    // Добавляем ребро
    void insertEdge(Node* from, Node* to, double weight) {
        from->neighbors.emplace_back(to, weight);
        to->neighbors.emplace_back(from, weight);
    }

    // Ищем узел по координатам к кеше
    Node* getNode(double lon, double lat) {
        std::string key = generateKey(lon, lat);
        auto it = node_map.find(key);
        return it != node_map.end() ? it->second : nullptr;
    }

    // Ближайший узел к заданным координатам
    Node* getClosestNode(double lon, double lat) {
        double minDist = std::numeric_limits<double>::max();
        Node* closestNode = nullptr;
        for (auto& node : nodes) {
            double distance = std::pow(node->longitude - lon, 2) + std::pow(node->latitude - lat, 2);
            if (distance < minDist) {
                minDist = distance;
                closestNode = node;
            }
        }
        return closestNode;
    }

    // Парсинг данных из файла и создание узлов и рёбер
    void loadFromFile(const std::string& filename) {
        std::ifstream file(filename);
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            char separator;
            double lon, lat;
            iss >> lon >> separator >> lat >> separator;
            Node* currentNode = getNode(lon, lat);
            if (!currentNode) {
                currentNode = new Node(nodes.size(), lon, lat);
                nodes.push_back(currentNode);
                node_map[generateKey(lon, lat)] = currentNode;
            }
            while (!iss.eof()) {
                double neighborLon, neighborLat, weight;
                iss >> neighborLon >> separator >> neighborLat >> separator >> weight >> separator;
                Node* neighborNode = getNode(neighborLon, neighborLat);
                if (!neighborNode) {
                    neighborNode = new Node(nodes.size(), neighborLon, neighborLat);
                    nodes.push_back(neighborNode);
                    node_map[generateKey(neighborLon, neighborLat)] = neighborNode;
                }
                insertEdge(currentNode, neighborNode, weight);
            }
        }
        file.close();
    }

    // Поиск в глубину (DFS)
    double depthFirstSearch(Node* start, Node* end) {
        std::stack<Node*> stack; // O(1) - создание стека
        std::unordered_map<Node*, Node*> previousNodes; // O(V) - память для хранения предыдущих узлов
        std::unordered_set<Node*> visitedNodes; // O(V) - память для хранения посещенных узлов
        std::unordered_map<Node*, double> distances; // O(V) - память для хранения расстояний
        stack.push(start);
        distances[start] = 0.0;
        while (!stack.empty()) {
            Node* current = stack.top(); // O(1) - извлечение узла из стека
            stack.pop();
            if (visitedNodes.count(current)) continue; // O(1) - проверка на посещение
            visitedNodes.insert(current); // O(1) - пометка узла как посещенного
            if (current == end) return distances[current]; // O(1) - проверка на конец
            for (const auto& neighbor : current->neighbors) { // O(E) - обход соседей
                if (!visitedNodes.count(neighbor.first)) { // O(1) - проверка на посещение соседа
                    double newDistance = distances[current] + neighbor.second; // O(1)
                    if (distances.find(neighbor.first) == distances.end() || newDistance < distances[neighbor.first]) { // O(1)
                        distances[neighbor.first] = newDistance; // O(1)
                        previousNodes[neighbor.first] = current; // O(1)
                        stack.push(neighbor.first); // O(1)
                    }
                }
            }
        }
        return -1.0; // Путь не найден
    }
    /*Итог: время O(V+E) V-вершины E-ребра; память O(V)*/


    // Поиск в ширину (BFS)
    double breadthFirstSearch(Node* start, Node* end) {
        std::queue<Node*> queue; // O(1) - создание очереди
        std::unordered_map<Node*, double> distances; // O(V) - память для хранения расстояний
        std::unordered_set<Node*> visitedNodes; // O(V) - память для хранения посещенных узлов
        queue.push(start);
        distances[start] = 0.0;
        while (!queue.empty()) {
            Node* current = queue.front(); // O(1) - извлечение узла из очереди
            queue.pop();
            if (current == end) return distances[current]; // O(1)
            for (const auto& neighbor : current->neighbors) { // O(E)
                if (!visitedNodes.count(neighbor.first)) { // O(1)
                    visitedNodes.insert(neighbor.first); // O(1)
                    distances[neighbor.first] = distances[current] + neighbor.second; // O(1)
                    queue.push(neighbor.first); // O(1)
                }
            }
        }
        return -1.0; // Путь не найден
    }
    /*Итог: время O(V+E) V-вершины E-ребра; память O(V)*/

    // Алгоритм Дейкстры для нахождения кратчайшего пути
    double dijkstra(Node* start, Node* end) {
        std::unordered_map<Node*, double> distances; // O(V)
        for (Node* node : nodes) {
            distances[node] = std::numeric_limits<double>::infinity(); // O(V)
        }
        distances[start] = 0.0;
        auto compare = [](const std::pair<Node*, double>& a, const std::pair<Node*, double>& b) {
            return a.second > b.second;
            };
        std::priority_queue<std::pair<Node*, double>,
            std::vector<std::pair<Node*, double>>,
            decltype(compare)> priorityQueue(compare); // O(V log V)
        priorityQueue.push({ start, 0.0 });
        while (!priorityQueue.empty()) {
            Node* current = priorityQueue.top().first; // O(1)
            priorityQueue.pop();
            if (current == end) return distances[current]; // O(1)
            for (const auto& neighbor : current->neighbors) {
                double newDistance = distances[current] + neighbor.second;
                if (newDistance < distances[neighbor.first]) {
                    distances[neighbor.first] = newDistance;
                    priorityQueue.push({ neighbor.first, newDistance }); // O(log V)
                }
            }
        }
        return -1.0; // Путь не найден
    }
    /*Итог: время O(E*log(V)) V-вершины E-ребра; память O(V)*/

};

// Тестикииии
void ttests() {
    Graph graph;
    graph.loadFromFile("C:\\Users\\ya\\Test_files\test.txt");
    Node* startNode = graph.getNode(0.0, 0.0);
    Node* endNode = graph.getNode(1.0, 1.0);
    if (startNode && endNode) {
        std::cout << "Testing DFS:\n";
        std::cout << "DFS Distance: " << graph.depthFirstSearch(startNode, endNode) << "\n";

        std::cout << "Testing BFS:\n";
        std::cout << "BFS Distance: " << graph.breadthFirstSearch(startNode, endNode) << "\n";

        std::cout << "Testing Dijkstra:\n";
        std::cout << "Dijkstra Distance: " << graph.dijkstra(startNode, endNode) << "\n";
    }
    else {
        std::cout << "!S||!E\n";
    }
}

// Основная функция программы
int main() {
    Graph graph;
    // Начало отсчета времени
    auto start_time_load = std::chrono::high_resolution_clock::now();
    graph.loadFromFile("C:\\Users\\ya\\Test_files\\spb_graph.txt");
    // Конец отсчета времени
    auto end_time_load = std::chrono::high_resolution_clock::now();
    auto elapsedTimeLoad = std::chrono::duration<double>(end_time_load - start_time_load).count();
    std::cout << "Time to load graph: " << elapsedTimeLoad << " seconds\n";
    double startLongitude = 30.309057;
    double startLatitude = 59.956577;
    double endLongitude = 30.303048;
    double endLatitude = 59.971823;
    Node* startPoint = graph.getClosestNode(startLongitude, startLatitude);
    Node* endPoint = graph.getClosestNode(endLongitude, endLatitude);
    if (startPoint && endPoint) {
        //Тут запускаем функции и считаем время
        auto start_time_dfs = std::chrono::high_resolution_clock::now();
        double dfsPathLength = graph.depthFirstSearch(startPoint, endPoint);
        auto end_time_dfs = std::chrono::high_resolution_clock::now();
        auto elapsedTimeDFS = std::chrono::duration<double>(end_time_dfs - start_time_dfs).count();
        std::cout << "DFS_path: " << dfsPathLength << "\n";
        std::cout << "DFS_t: " << elapsedTimeDFS << " seconds\n";
        auto start_time_bfs = std::chrono::high_resolution_clock::now();
        double bfsPathLength = graph.breadthFirstSearch(startPoint, endPoint);
        auto end_time_bfs = std::chrono::high_resolution_clock::now();
        auto elapsedTimeBFS = std::chrono::duration<double>(end_time_bfs - start_time_bfs).count();
        std::cout << "BFS_path: " << bfsPathLength << "\n";
        std::cout << "BFS_t: " << elapsedTimeBFS << " seconds\n";
        auto start_time_dijkstra = std::chrono::high_resolution_clock::now();
        double dijkstraPathLength = graph.dijkstra(startPoint, endPoint);
        auto end_time_dijkstra = std::chrono::high_resolution_clock::now();
        auto elapsedTimeDijkstra = std::chrono::duration<double>(end_time_dijkstra - start_time_dijkstra).count();
        std::cout << "Dijkstra_min_path: " << dijkstraPathLength << "\n";
        std::cout << "Dijkstra_t: " << elapsedTimeDijkstra << " seconds\n";
    }
    else {
        std::cout << "!S||!E\n";
    }

    ttests();
    return 0;
}
