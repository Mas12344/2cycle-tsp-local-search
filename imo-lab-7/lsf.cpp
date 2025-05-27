#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <time.h>
#include <array>
#include <stdint.h>
#include <list>
#include <limits>
#include <set>
#include <chrono>
#include <numeric>
#include <unordered_set>
#include <fstream>
#include <string>
#include <sstream>
#include <thread>
#include <mutex>


#if 1
#define NODE_COUNT 200
#else
#define NODE_COUNT 12
#endif
#define CYCLE_SIZE (NODE_COUNT / 2)

#define MAX_THREADS 6

extern "C" {
    // MARK: Helpers
struct Instance { 
    int dist[NODE_COUNT][NODE_COUNT];
};

struct SolutionOut {
    int cycle1[CYCLE_SIZE];
    int cycle2[CYCLE_SIZE];
    int value;
}; 

struct Point {
    double x;
    double y;
};

std::vector<Point> loadData(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<Point> points;
    bool dataNow = false;

    while (std::getline(file, line)) {
        if (line == "EOF") break;
        if (line == "NODE_COORD_SECTION") {
            dataNow = true;
            continue;
        }
        if (dataNow) {
            std::istringstream iss(line);
            int idx;
            double x, y;
            iss >> idx >> x >> y;
            points.push_back({x, y});
        }
    }

    return points;
}

int euclidDistance(const Point& a, const Point& b) {
    return static_cast<int>(std::round(std::hypot(a.x - b.x, a.y - b.y)));
}

Instance genInst(const std::vector<Point>& points) {
    Instance inst;
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = 0; j < points.size(); ++j) {
            if (i == j)
                inst.dist[i][j] = 0;
            else
                inst.dist[i][j] = euclidDistance(points[i], points[j]);
        }
    }
    return inst;
}


int absMod(int val, int mod) {
    while(val<0) val+=mod;
    return val%mod;
}

bool sameCycle(int i1, int i2) {
    if (i1 < CYCLE_SIZE) return i2 < CYCLE_SIZE;
    else return i2 >= CYCLE_SIZE;
}

int indexOf(const int c[NODE_COUNT], int vertex) {
    for (int i = 0; i < NODE_COUNT; i++) {
        if (c[i] == vertex) return i;
    }
    return -1;
}


class Solution {
public:
    int cycles[NODE_COUNT];
    int value;

    Solution(SolutionOut* so) {
        for (int i=0; i < CYCLE_SIZE; i++) {
            this->cycles[i] = so->cycle1[i];
            this->cycles[i+CYCLE_SIZE] = so->cycle2[i];
        }
        this->value = so->value;
    };

    void backToOut(SolutionOut* so) {
        for (int i=0; i < CYCLE_SIZE; i++) {
            so->cycle1[i] = this->cycles[i];
            so->cycle2[i] = this->cycles[i+CYCLE_SIZE];
        }
        so->value = this->value;
    };

    void eval(const Instance& instance) {
        int delta = 0;
        for (int i=1; i < CYCLE_SIZE; i++) {
            int a = this->cycles[i-1];
            int b = this->cycles[i];
            int c = this->cycles[i-1+CYCLE_SIZE];
            int d = this->cycles[i+CYCLE_SIZE];
            delta += instance.dist[a][b];
            delta += instance.dist[c][d];
        }
        delta += instance.dist[this->cycles[CYCLE_SIZE-1]][this->cycles[0]];
        delta += instance.dist[this->cycles[NODE_COUNT-1]][this->cycles[CYCLE_SIZE]];
        this->value = delta;
    };   

};

void weighted_regret_heuristic(const Instance& instance, SolutionOut* outSolution, double alpha = 1.0) {
    int indices[NODE_COUNT];
    for (int i = 0; i < NODE_COUNT; ++i) indices[i] = i;

    static std::random_device rd;
    static std::mt19937 g(rd());
    std::shuffle(indices, indices + NODE_COUNT, g);

    int target_size1 = (NODE_COUNT + 1) / 2;
    int target_size2 = NODE_COUNT / 2;

    int cycle1[CYCLE_SIZE];
    int cycle2[CYCLE_SIZE];
    int size1 = 2, size2 = 2;

    bool visited[NODE_COUNT] = {false};
    cycle1[0] = indices[0];
    cycle1[1] = indices[1];
    cycle2[0] = indices[2];
    cycle2[1] = indices[3];
    visited[indices[0]] = true;
    visited[indices[1]] = true;
    visited[indices[2]] = true;
    visited[indices[3]] = true;

    while (true) {
        bool all_visited = true;
        for (int i = 0; i < NODE_COUNT; ++i) {
            if (!visited[i]) {
                all_visited = false;
                break;
            }
        }
        if (all_visited) break;

        double best_score = -std::numeric_limits<double>::infinity();
        int best_node = -1;
        int best_cycle = -1;
        int best_pos = -1;

        for (int idx = 0; idx < NODE_COUNT; ++idx) {
            int node = indices[idx];
            if (visited[node]) continue;

            double best_cost = std::numeric_limits<double>::infinity();
            double second_best_cost = std::numeric_limits<double>::infinity();
            int best_position = -1;
            int best_cycle_choice = -1;

            if (size1 < target_size1) {
                for (int i = 0; i < size1; ++i) {
                    int from = cycle1[i];
                    int to = cycle1[(i + 1) % size1];
                    int cost = instance.dist[from][node] + instance.dist[node][to] - instance.dist[from][to];

                    if (cost < best_cost) {
                        second_best_cost = best_cost;
                        best_cost = cost;
                        best_position = i + 1;
                        best_cycle_choice = 0;
                    } else if (cost < second_best_cost) {
                        second_best_cost = cost;
                    }
                }
            } else if (size2 < target_size2) {
                for (int i = 0; i < size2; ++i) {
                    int from = cycle2[i];
                    int to = cycle2[(i + 1) % size2];
                    int cost = instance.dist[from][node] + instance.dist[node][to] - instance.dist[from][to];

                    if (cost < best_cost) {
                        second_best_cost = best_cost;
                        best_cost = cost;
                        best_position = i + 1;
                        best_cycle_choice = 1;
                    } else if (cost < second_best_cost) {
                        second_best_cost = cost;
                    }
                }
            } else {
                for (int c = 0; c < 2; ++c) {
                    int* cycle = (c == 0) ? cycle1 : cycle2;
                    int size = (c == 0) ? size1 : size2;
                    for (int i = 0; i < size; ++i) {
                        int from = cycle[i];
                        int to = cycle[(i + 1) % size];
                        int cost = instance.dist[from][node] + instance.dist[node][to] - instance.dist[from][to];

                        if (cost < best_cost) {
                            second_best_cost = best_cost;
                            best_cost = cost;
                            best_position = i + 1;
                            best_cycle_choice = c;
                        } else if (cost < second_best_cost) {
                            second_best_cost = cost;
                        }
                    }
                }
            }

            double regret = (second_best_cost < std::numeric_limits<double>::infinity()) ? second_best_cost - best_cost : best_cost;
            double score = -alpha * best_cost + alpha * regret;

            if (score > best_score) {
                best_score = score;
                best_node = node;
                best_cycle = best_cycle_choice;
                best_pos = best_position;
            }
        }

        if (best_node != -1) {
            if (best_cycle == 0) {
                for (int i = size1; i > best_pos; --i) {
                    cycle1[i] = cycle1[i - 1];
                }
                cycle1[best_pos] = best_node;
                ++size1;
            } else if (best_cycle == 1) {
                for (int i = size2; i > best_pos; --i) {
                    cycle2[i] = cycle2[i - 1];
                }
                cycle2[best_pos] = best_node;
                ++size2;
            }
            visited[best_node] = true;
        } else {
            break;
        }
    }


    for (int i = 0; i < CYCLE_SIZE; ++i) {
        outSolution->cycle1[i] = cycle1[i];
        outSolution->cycle2[i] = cycle2[i];
    }
    Solution solution(outSolution);
    solution.eval(instance);
    solution.backToOut(outSolution);
}

void przesunMinusJedynki(std::vector<int>& cykl) {
    auto it = std::stable_partition(cykl.begin(), cykl.end(), [](int x) {
        return x != -1;
    });
    std::fill(it, cykl.end(), -1);
}

void weighted_regret_heuristic_5(const Instance& instance, SolutionOut* outSolution, double alpha = 1.0) {
    int indices[NODE_COUNT];
    for (int i = 0; i < NODE_COUNT; ++i) indices[i] = i;

    static std::random_device rd;
    static std::mt19937 g(rd());
    std::shuffle(indices, indices + NODE_COUNT, g);

    auto it1 = std::stable_partition(outSolution->cycle1, outSolution->cycle1+CYCLE_SIZE, [](int x) {
        return x != -1;
    });
    std::fill(it1, outSolution->cycle1+CYCLE_SIZE, -1);
    
    auto it2 = std::stable_partition(outSolution->cycle2, outSolution->cycle2+CYCLE_SIZE, [](int x) {
        return x != -1;
    });
    std::fill(it2, outSolution->cycle2+CYCLE_SIZE, -1);
    

    int target_size = (NODE_COUNT ) / 2;

    int cycle1[CYCLE_SIZE];
    int cycle2[CYCLE_SIZE];
    int size1 = 0, size2 = 0;

    bool visited[NODE_COUNT] = {false};
    for (int i = 0; i < CYCLE_SIZE; ++i) {
        int node1 = outSolution->cycle1[i];
        if (node1 != -1) {
            cycle1[size1++] = node1;
            visited[node1] = true;
        }
        int node2 = outSolution->cycle2[i];
        if (node2 != -1) {
            cycle2[size2++] = node2;
            visited[node2] = true;
        }
    }

    if (size1 == 0) {
        for (int i = 0; i < NODE_COUNT; i++){
            if (!visited[i]) {
                cycle1[size1++] = i;
                visited[i] = true;
            }
            if (size1==2) {
                break;
            }
        }
    }

    if (size2 == 0) {
        for (int i = 0; i < NODE_COUNT; i++){
            if (!visited[i]) {
                cycle2[size2++] = i;
                visited[i] = true;
            }
            if (size2==2) {
                break;
            }
        }
    }



    while (true) {
        bool all_visited = true;
        for (int i = 0; i < NODE_COUNT; ++i) {
            if (!visited[i]) {
                all_visited = false;
                break;
            }
        }
        if (all_visited) break;

        double best_score = -std::numeric_limits<double>::infinity();
        int best_node = -1;
        int best_cycle = -1;
        int best_pos = -1;

        for (int idx = 0; idx < NODE_COUNT; ++idx) {
            int node = indices[idx];
            if (visited[node]) continue;

            double best_cost = std::numeric_limits<double>::infinity();
            double second_best_cost = std::numeric_limits<double>::infinity();
            int best_position = -1;
            int best_cycle_choice = -1;

            if (size1 < target_size) {
                for (int i = 0; i < size1; ++i) {
                    int from = cycle1[i];
                    int to = cycle1[(i + 1) % size1];
                    int cost = instance.dist[from][node] + instance.dist[node][to] - instance.dist[from][to];

                    if (cost < best_cost) {
                        second_best_cost = best_cost;
                        best_cost = cost;
                        best_position = i + 1;
                        best_cycle_choice = 0;
                    } else if (cost < second_best_cost) {
                        second_best_cost = cost;
                    }
                }
            } else if (size2 < target_size) {
                for (int i = 0; i < size2; ++i) {
                    int from = cycle2[i];
                    int to = cycle2[(i + 1) % size2];
                    int cost = instance.dist[from][node] + instance.dist[node][to] - instance.dist[from][to];

                    if (cost < best_cost) {
                        second_best_cost = best_cost;
                        best_cost = cost;
                        best_position = i + 1;
                        best_cycle_choice = 1;
                    } else if (cost < second_best_cost) {
                        second_best_cost = cost;
                    }
                }
            } else {
                for (int c = 0; c < 2; ++c) {
                    int* cycle = (c == 0) ? cycle1 : cycle2;
                    int size = (c == 0) ? size1 : size2;
                    for (int i = 0; i < size; ++i) {
                        int from = cycle[i];
                        int to = cycle[(i + 1) % size];
                        int cost = instance.dist[from][node] + instance.dist[node][to] - instance.dist[from][to];

                        if (cost < best_cost) {
                            second_best_cost = best_cost;
                            best_cost = cost;
                            best_position = i + 1;
                            best_cycle_choice = c;
                        } else if (cost < second_best_cost) {
                            second_best_cost = cost;
                        }
                    }
                }
            }

            double regret = (second_best_cost < std::numeric_limits<double>::infinity()) ? second_best_cost - best_cost : best_cost;
            double score = -alpha * best_cost + alpha * regret;

            if (score > best_score) {
                best_score = score;
                best_node = node;
                best_cycle = best_cycle_choice;
                best_pos = best_position;
            }
        }
        if (best_node < 0 || best_node > NODE_COUNT) __debugbreak();

        if (best_node != -1) {
            if (best_cycle == 0) {
                for (int i = size1; i > best_pos; --i) {
                    cycle1[i] = cycle1[i - 1];
                }
                cycle1[best_pos] = best_node;
                ++size1;
            } else if (best_cycle == 1) {
                for (int i = size2; i > best_pos; --i) {
                    cycle2[i] = cycle2[i - 1];
                }
                cycle2[best_pos] = best_node;
                ++size2;
            }
            visited[best_node] = true;
        } else {
            break;
        }
    }


    for (int i = 0; i < CYCLE_SIZE; ++i) {
        outSolution->cycle1[i] = cycle1[i];
        outSolution->cycle2[i] = cycle2[i];
    }
    Solution solution(outSolution);
    solution.eval(instance);
    solution.backToOut(outSolution);
}

// MARK: weighted_regret_heuristic
void weighted_regret_heuristic_4(const Instance& instance, SolutionOut* outSolution, double alpha = 1.0) {
    int indices[NODE_COUNT];
    for (int i = 0; i < NODE_COUNT; ++i) indices[i] = i;

    static std::random_device rd;
    static std::mt19937 g(rd());
    std::shuffle(indices, indices + NODE_COUNT, g);

    int target_size = (NODE_COUNT ) / 2;

    int cycle1[CYCLE_SIZE];
    int cycle2[CYCLE_SIZE];
    int size1 = 0, size2 = 0;

    bool visited[NODE_COUNT] = {false};
    for (int i = 0; i < CYCLE_SIZE; ++i) {
        int node1 = outSolution->cycle1[i];
        if (node1 != -1) {
            cycle1[size1++] = node1;
            visited[node1] = true;
        }
        int node2 = outSolution->cycle2[i];
        if (node2 != -1) {
            cycle2[size2++] = node2;
            visited[node2] = true;
        }
    }


    while (true) {
        bool all_visited = true;
        for (int i = 0; i < NODE_COUNT; ++i) {
            if (!visited[i]) {
                all_visited = false;
                break;
            }
        }
        if (all_visited) break;

        double best_score = -std::numeric_limits<double>::infinity();
        int best_node = -1;
        int best_cycle = -1;
        int best_pos = -1;

        for (int idx = 0; idx < NODE_COUNT; ++idx) {
            int node = indices[idx];
            if (visited[node]) continue;

            double best_cost = std::numeric_limits<double>::infinity();
            double second_best_cost = std::numeric_limits<double>::infinity();
            int best_position = -1;
            int best_cycle_choice = -1;

            if (size1 < target_size) {
                for (int i = 0; i < size1; ++i) {
                    int from = cycle1[i];
                    int to = cycle1[(i + 1) % size1];
                    int cost = instance.dist[from][node] + instance.dist[node][to] - instance.dist[from][to];

                    if (cost < best_cost) {
                        second_best_cost = best_cost;
                        best_cost = cost;
                        best_position = i + 1;
                        best_cycle_choice = 0;
                    } else if (cost < second_best_cost) {
                        second_best_cost = cost;
                    }
                }
            } else if (size2 < target_size) {
                for (int i = 0; i < size2; ++i) {
                    int from = cycle2[i];
                    int to = cycle2[(i + 1) % size2];
                    int cost = instance.dist[from][node] + instance.dist[node][to] - instance.dist[from][to];

                    if (cost < best_cost) {
                        second_best_cost = best_cost;
                        best_cost = cost;
                        best_position = i + 1;
                        best_cycle_choice = 1;
                    } else if (cost < second_best_cost) {
                        second_best_cost = cost;
                    }
                }
            } else {
                for (int c = 0; c < 2; ++c) {
                    int* cycle = (c == 0) ? cycle1 : cycle2;
                    int size = (c == 0) ? size1 : size2;
                    for (int i = 0; i < size; ++i) {
                        int from = cycle[i];
                        int to = cycle[(i + 1) % size];
                        int cost = instance.dist[from][node] + instance.dist[node][to] - instance.dist[from][to];

                        if (cost < best_cost) {
                            second_best_cost = best_cost;
                            best_cost = cost;
                            best_position = i + 1;
                            best_cycle_choice = c;
                        } else if (cost < second_best_cost) {
                            second_best_cost = cost;
                        }
                    }
                }
            }

            double regret = (second_best_cost < std::numeric_limits<double>::infinity()) ? second_best_cost - best_cost : best_cost;
            double score = -alpha * best_cost + alpha * regret;

            if (score > best_score) {
                best_score = score;
                best_node = node;
                best_cycle = best_cycle_choice;
                best_pos = best_position;
            }
        }

        if (best_node != -1) {
            if (best_cycle == 0) {
                for (int i = size1; i > best_pos; --i) {
                    cycle1[i] = cycle1[i - 1];
                }
                cycle1[best_pos] = best_node;
                ++size1;
            } else if (best_cycle == 1) {
                for (int i = size2; i > best_pos; --i) {
                    cycle2[i] = cycle2[i - 1];
                }
                cycle2[best_pos] = best_node;
                ++size2;
            }
            visited[best_node] = true;
        } else {
            break;
        }
    }


    for (int i = 0; i < CYCLE_SIZE; ++i) {
        outSolution->cycle1[i] = cycle1[i];
        outSolution->cycle2[i] = cycle2[i];
    }
    Solution solution(outSolution);
    solution.eval(instance);
    solution.backToOut(outSolution);
}


int evalEdgeSwap(const Instance& instance, const int cycles[NODE_COUNT], int i1, int i2) {
    int delta = 0, a, b, c, d;
    if (i1 < CYCLE_SIZE) {
        a = cycles[i1];
        b = cycles[absMod(i1+1, CYCLE_SIZE)];
        c = cycles[i2];
        d = cycles[absMod(i2+1, CYCLE_SIZE)];
    } else {

        int i1n = absMod(i1+1, CYCLE_SIZE) + CYCLE_SIZE;
        int i2n = absMod(i2+1, CYCLE_SIZE) + CYCLE_SIZE;
        a = cycles[i1];
        b = cycles[i1n];
        c = cycles[i2];
        d = cycles[i2n];
    }
    delta += instance.dist[a][c];
    delta += instance.dist[b][d];
    delta -= instance.dist[a][b];
    delta -= instance.dist[c][d];
    return delta;
}

// i1 i i2 należą do tego samego cyklu, i1 != i2, bo to nie ma sensu
void swapEdge(int* c, int i1, int i2) {
    int f = i1, s = i2;
    if (f >= CYCLE_SIZE) {
        c += CYCLE_SIZE;
        f %= CYCLE_SIZE;
        s %= CYCLE_SIZE;
    }
    f = absMod(f+1, CYCLE_SIZE);

    int len = (s - f + CYCLE_SIZE) % CYCLE_SIZE + 1;
    for (int i = 0; i < len / 2; ++i) {
        int a = (f + i) % CYCLE_SIZE;
        int b = (s - i + CYCLE_SIZE) % CYCLE_SIZE;
        std::swap(c[a], c[b]);
    }
}

int evalVertexSwap(const Instance& instance, const int cycles[NODE_COUNT], int i1, int i2) {
    int temp1n = absMod(i1+1, CYCLE_SIZE);
    int temp2n = absMod(i2+1, CYCLE_SIZE);
    int temp1p = absMod(i1-1, CYCLE_SIZE);
    int temp2p = absMod(i2-1, CYCLE_SIZE);
    int i1n = i1 < CYCLE_SIZE ? temp1n : temp1n + CYCLE_SIZE;
    int i2n = i2 < CYCLE_SIZE ? temp2n : temp2n + CYCLE_SIZE;
    int i1p = i1 < CYCLE_SIZE ? temp1p : temp1p + CYCLE_SIZE;
    int i2p = i2 < CYCLE_SIZE ? temp2p : temp2p + CYCLE_SIZE;

    int delta = 0;
    int a = cycles[i1p];
    int b = cycles[i1];
    int c = cycles[i1n];
    int d = cycles[i2p];
    int e = cycles[i2];
    int f = cycles[i2n];



    if(sameCycle(i1, i2) && (i1==i2p || i2==i1p)) {
        if(i1==i2p) {
            delta += instance.dist[a][c];
            delta += instance.dist[b][f];
            delta -= instance.dist[a][b];
            delta -= instance.dist[e][f];

        } else {
            delta += instance.dist[d][b];
            delta += instance.dist[e][c];
            delta -= instance.dist[d][e];
            delta -= instance.dist[b][c];
        }
    } else {
        delta += instance.dist[b][f];
        delta += instance.dist[a][e];
        delta += instance.dist[e][c];
        delta += instance.dist[d][b];
        delta -= instance.dist[a][b];
        delta -= instance.dist[b][c];
        delta -= instance.dist[d][e];
        delta -= instance.dist[e][f];
    }

    return delta;
}


void swapVertex(int* c, int i1, int i2) {
    int t = c[i1];
    c[i1] = c[i2];
    c[i2] = t;
}

void baseline_alg(const Instance& instance, SolutionOut* outSolution) {
    Solution sol(outSolution);
    sol.eval(instance);
    
    bool improved;
    do {
        improved = false;
        bool edgeMove = true;
        std::pair<int,int> bestMove;
        int bestDelta = 0;

        // Przegląd wszystkich możliwych ruchów (zamian)
        for (int c = 0; c < 2; c++) {
            for (int i = 0; i < CYCLE_SIZE; i++) {
                for (int j = 0; j < CYCLE_SIZE; j++) {
                    if(i >= j ) continue;
                    if(absMod(i+1, CYCLE_SIZE) == j || absMod(j+1, CYCLE_SIZE) == i) continue;

                    int delta = evalEdgeSwap(instance, sol.cycles, i+c*CYCLE_SIZE, j+c*CYCLE_SIZE);

                    if (delta < 0 && delta < bestDelta) {
                        bestDelta = delta;
                        bestMove = {i+c*CYCLE_SIZE, j+c*CYCLE_SIZE};
                        improved = true;
                    }
                }
            }
        }

        for (int i = 0; i < CYCLE_SIZE; i++) {
            for (int j = 0; j < CYCLE_SIZE; j++) {
                int delta = evalVertexSwap(instance, sol.cycles, i, j+CYCLE_SIZE);

                if (delta < 0 && delta < bestDelta) {
                    edgeMove = false;
                    bestDelta = delta;
                    bestMove = {i, j+CYCLE_SIZE};
                    improved = true;
                }
            }
        }


        // Jeśli znaleziono lepsze rozwiązanie, przejdź do niego
        if (improved) {
            auto[f, s] = bestMove;

            if(edgeMove) {
                swapEdge(sol.cycles, f, s);
            } else {
                swapVertex(sol.cycles, f, s);
            }

            sol.value += bestDelta;
        }

    } while (improved);

    sol.backToOut(outSolution);
}

SolutionOut random_solution() {
    SolutionOut solOut;
    int nodes[NODE_COUNT];
    for (int i = 0; i < NODE_COUNT; i++) nodes[i] = i;
    
    std::shuffle(nodes, nodes+NODE_COUNT, std::default_random_engine(std::rand()));
    for (int i = 0; i < CYCLE_SIZE; i++) {
        solOut.cycle1[i] = nodes[i];
        solOut.cycle2[i] = nodes[i + CYCLE_SIZE];
    }
    solOut.value = -1;
    return solOut;
}

inline long time_elapsed(std::chrono::high_resolution_clock::time_point startTime) {
    return  std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime).count();
}

long multiple_start(const Instance& instance, SolutionOut* best, int iterCount=200) {
    best->value = std::numeric_limits<int>::max();

    auto start = std::chrono::high_resolution_clock::now();
    for(int i=0; i<iterCount; i++) {
        SolutionOut newSol = random_solution();

        baseline_alg(instance, &newSol);
        if(newSol.value < best->value) *best = newSol;
    }
    return time_elapsed(start);
}

SolutionOut peturbate(SolutionOut* sol, int count=5) {
    SolutionOut newSol = *sol;

    for(int i=0;i<count; i++) {
        int a = std::rand() % CYCLE_SIZE;
        int b = std::rand() % CYCLE_SIZE;
        if(a == b) a = (a+1) % CYCLE_SIZE;
        
        if(std::rand()%2) { //edge swap
            int cycle = std::rand()%2;
            if(cycle) std::reverse(newSol.cycle1 + a + 1, newSol.cycle1 + b+ 1);
            else std::reverse(newSol.cycle2 + a + 1, newSol.cycle2 + b + 1);
        } else { //vertex swap
            std::swap(newSol.cycle1[a], newSol.cycle1[b]);
        }
    }
    return newSol;
}

int iterated_alg(const Instance& instance, SolutionOut* sol, long timeLimit) {
    baseline_alg(instance, sol);
    int c = 0;
    auto start = std::chrono::high_resolution_clock::now();
    while(time_elapsed(start) < timeLimit) {
        SolutionOut newSol = peturbate(sol);
        baseline_alg(instance, &newSol);

        if(newSol.value < sol->value) *sol = newSol;
        c++;
    }
    return c;
}

//#MARK: destroy_repair
SolutionOut destroy_repair(const Instance& instance, SolutionOut* sol, float destructionFactor = 0.6) {
    SolutionOut newSol = *sol;
    int nodes = CYCLE_SIZE * destructionFactor;

    int bestRegions[2];
    for(int c=0; c<2; c++) {
        int* cycle1 = c ? newSol.cycle1 : newSol.cycle2;
        int* cycle2 = c ? newSol.cycle2 : newSol.cycle1;


        int regions[CYCLE_SIZE];
        int bestRegion = 0;
        int bestSize = 0;
        int currentSize = 0;
        for(int i=0; i<nodes; i++) {
            int best = std::numeric_limits<int>::max();
            for(int j=0; j<CYCLE_SIZE; j++) {
                if(best > instance.dist[cycle1[i]][cycle2[j]]) {
                    best = instance.dist[cycle1[i]][cycle2[j]];
                }
            }
            regions[i] = best;
            currentSize += best;
        }
        bestSize = currentSize;
        for(int i=nodes; i<CYCLE_SIZE; i++) {
            int best = std::numeric_limits<int>::max();
            for(int j=0; j<CYCLE_SIZE; j++) {
                if(best > instance.dist[cycle1[i]][cycle2[j]]) {
                    best = instance.dist[cycle1[i]][cycle2[j]];
                }
            }
            regions[i] = best;
            currentSize += best - regions[i-nodes];
            if(currentSize < bestSize) {
                bestSize = currentSize;
                bestRegion = i - nodes + 1;
            }
        }
        bestRegions[c] = bestRegion;

    }
    

    int c = 0;
    for (const auto& bestregion : bestRegions){
        int* cycle = c == 0 ? newSol.cycle1 : newSol.cycle2;
        for(int i=bestregion; i<bestregion+nodes; i++) {
            cycle[i] = -1;
        }
        c++;
    }


    weighted_regret_heuristic_4(instance, &newSol);

    return newSol;
}

int destroy_repair_alg(const Instance& instance, SolutionOut* sol, long timeLimit, bool doLocalSearch=false, float destructionFactor = 0.4) {
    baseline_alg(instance, sol);

    auto start = std::chrono::high_resolution_clock::now();
    int c = 0;
    while(time_elapsed(start) < timeLimit) {
        SolutionOut newSol = destroy_repair(instance, sol, destructionFactor);
        if(doLocalSearch) baseline_alg(instance, &newSol);
        if(newSol.value < sol->value) *sol = newSol;
        c++;
    }
    return c;
}

bool operator<(const SolutionOut& a, const SolutionOut& b) {
    return a.value < b.value;
}

struct Edge {
    int v1;
    int v2;

    bool operator< (const Edge& other) const {
        return std::tie(v1, v2) < std::tie(other.v1, other.v2);
    };
};

SolutionOut recombine_parents(const Instance& instance, const std::vector<SolutionOut>& parents) {
    SolutionOut newSol = parents[0];
    std::set<Edge> parent0_egdes;
    std::set<Edge> parent1_egdes;
    for (int i = 0; i < CYCLE_SIZE; i++) {
        parent0_egdes.insert({parents[0].cycle1[i], parents[0].cycle1[(i+1)%CYCLE_SIZE]}); 
        parent0_egdes.insert({parents[0].cycle2[i], parents[0].cycle2[(i+1)%CYCLE_SIZE]});
        parent1_egdes.insert({parents[1].cycle1[i], parents[1].cycle1[(i+1)%CYCLE_SIZE]}); 
        parent1_egdes.insert({parents[1].cycle2[i], parents[1].cycle2[(i+1)%CYCLE_SIZE]}); 
    }

    std::set<Edge> parents_intersection;
    std::set_intersection(parent0_egdes.begin(), parent0_egdes.end(),
        parent1_egdes.begin(), parent1_egdes.end(), std::inserter(parents_intersection, parents_intersection.begin()));

    std::set<Edge> parent0_unique;
    std::set_difference(parent0_egdes.begin(), parent0_egdes.end(),
        parents_intersection.begin(), parents_intersection.end(), std::inserter(parent0_unique, parent0_unique.begin()));

    std::array<int, NODE_COUNT> how_many = {0};
    for (const auto& it : parent0_unique) {
        how_many[it.v1]++;
        how_many[it.v2]++;
    }

    std::set<int> vertices_to_remove;
    for (int i = 0; i < NODE_COUNT; i++) {
        if(how_many[i]>0) vertices_to_remove.insert(i);
    }

#define MUTATION_CHANCE_PER_VERTEX 0.05
    std::uniform_real_distribution<> dis(0.0, 1.0);
    static std::random_device rd;
    static std::mt19937 g(rd());

    for (int i = 0; i < CYCLE_SIZE; i++) {
        for (const auto& v : vertices_to_remove) {
            //if(i==0) std::cout << v << " ";
            if(newSol.cycle1[i] == v) newSol.cycle1[i] = -1;            
            if(newSol.cycle2[i] == v) newSol.cycle2[i] = -1;
        }
        if (dis(g) <= MUTATION_CHANCE_PER_VERTEX) {
            newSol.cycle1[i] = -1;
        }
        if (dis(g) <= MUTATION_CHANCE_PER_VERTEX) {
            newSol.cycle2[i] = -1;
        }


    }

    weighted_regret_heuristic_5(instance, &newSol);
    return newSol;
}

#define POPULATION_SIZE 20

int hybrid_alg(const Instance& instance, SolutionOut* sol, int timeLimit, bool doLocalSearch=false) {
    auto start = std::chrono::high_resolution_clock::now();
    static std::random_device rd;
    static std::mt19937 g(rd());
    std::vector<SolutionOut> parents;
    parents.reserve(2);

    std::set<SolutionOut> population;
    while (population.size() < POPULATION_SIZE/2) {
        SolutionOut newSol = random_solution();
        baseline_alg(instance, &newSol);
        population.insert(newSol);
    }
    while (population.size() < POPULATION_SIZE) {
        SolutionOut newSol = random_solution();
        weighted_regret_heuristic(instance, &newSol);
        baseline_alg(instance, &newSol);
        population.insert(newSol);
    }
    int c = 0;
    while(time_elapsed(start) < timeLimit) {
        while (population.size() < POPULATION_SIZE+1) {
            if (time_elapsed(start) >= timeLimit) break;
            parents.clear();
            std::sample(population.begin(), population.end(), std::back_inserter(parents), 2, g);
            SolutionOut child = recombine_parents(instance, parents);
            if(doLocalSearch) baseline_alg(instance, &child);
            population.insert(child);
        }
        if(population.size()>0) population.erase(--population.end());
        c++;
    }

    *sol = *population.begin();
    return c;
}

std::mutex m;

void thread_alg(const Instance& instance, std::set<SolutionOut>& solutions, int timeLimit, bool doLocalSearch=false) {
    SolutionOut s;
    hybrid_alg(instance, &s, timeLimit, doLocalSearch);
    //std::lock_guard<std::mutex> lock(m);
    solutions.insert(s);
}

void mt_alg(const Instance& instance, SolutionOut* sol, int timeLimit, bool doLocalSearch=false) {
    std::set<SolutionOut> allSolutions;
    std::vector<std::thread> threads;
    for (int i = 0; i < MAX_THREADS; i++) {
        threads.push_back(std::thread(thread_alg, std::ref(instance), std::ref(allSolutions), timeLimit, doLocalSearch));
    }
    for (auto& th : threads) {
            th.join();
    }

    *sol = *allSolutions.begin();

}

int main() {
    auto kroA = loadData("kroA200.tsp");
    auto kroB = loadData("kroB200.tsp");
    Instance instance = genInst(kroA);

    int lowest = 40000;

    for (int i = 0; i < 100; i ++) {
        SolutionOut solOut = random_solution();
        mt_alg(instance, &solOut, 5500, true);
        std::cout<<i <<": Score: "<<solOut.value<<std::endl;

        if (solOut.value < lowest) {
            lowest = solOut.value;
        }
    }

    std::cout<<"Best: "<<lowest<<std::endl;

}

}