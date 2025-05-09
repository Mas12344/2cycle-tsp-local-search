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

#define NEIGHBOURHOOD_SIZE 10
#define NODE_COUNT 200
#define CYCLE_SIZE (NODE_COUNT / 2)

extern "C" {


    typedef std::array<std::array<int, NEIGHBOURHOOD_SIZE>, NODE_COUNT> Neighbourhood;

    struct Instance { 
        int dist[NODE_COUNT][NODE_COUNT];
    };

    struct SolutionOut {
        int cycle1[CYCLE_SIZE];
        int cycle2[CYCLE_SIZE];
        int value;
    };

    enum Applicability {
        APPLICABLE,
        FUTURE_MAYBE,
        NEEDS_REVERSAL,
        NOT_APPLICABLE
    };

    struct Move {
        int delta;
        bool isEdgeMove;
        int a, b, c, d, e, f; // zapamiętane krawędzie

        void canon() {
            if(isEdgeMove) {
                int t = std::min({a,b,c,d});
                int i=0;
                for (const auto& v : {a,b,c,d}) {
                    if (v==t) break;
                    i++;
                }

                switch (i) {
                    case 0:
                        break;
                    case 1:
                        std::swap(a,b);
                        std::swap(c,d);
                        break;
                    case 2:
                        std::swap(a,c);
                        std::swap(b,d);
                        break;
                    case 3:
                        std::swap(a,d);
                        std::swap(b,c);
                        break;
                }
            } else {
                if (a > c) {
                    std::swap(a,c);
                }
                if (d > f) {
                    std::swap(d,f);
                }
                if (a > d) {
                    std::swap(a,d);
                    std::swap(b,e);
                    std::swap(c,f);
                }
            }
        }

    //#define CONTAINER_LIST

#ifdef CONTAINER_LIST
        bool operator<(const Move& other) const {
            return delta < other.delta;
        }
#else
        bool operator<(const Move& other) const {
            return std::tie(delta, isEdgeMove, a, b, c, d, e, f) <
                   std::tie(other.delta, other.isEdgeMove, other.a, other.b, other.c, other.d, other.e, other.f);

        }
#endif
    };
    


#ifdef CONTAINER_LIST
    typedef std::list<Move> MoveContainer;
    void add_move(MoveContainer& LM, Move a) {
        LM.push_back(a);
    }
#else
    typedef std::set<Move> MoveContainer;
    void add_move(MoveContainer& LM, Move a) {
        a.canon();
        LM.insert(a);
    }

#endif

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

    enum edge_type {
        direct,
        reversed,
        no_edge
    };

    edge_type check_edge(const int* cycles, int index_u, int v) {
        int v1 = cycles[(index_u + 1) % CYCLE_SIZE + (index_u >= CYCLE_SIZE ? CYCLE_SIZE : 0)];
        int v2 = cycles[(index_u - 1 + CYCLE_SIZE) % CYCLE_SIZE + (index_u >= CYCLE_SIZE ? CYCLE_SIZE : 0)];
        if (v==v1) return direct;
        if (v==v2) return reversed;
        return no_edge;
    }

    Applicability isEdgeMoveApplicable(const Move& m, const int* cycles) {
        int index_a = indexOf(cycles, m.a);
        int index_c = indexOf(cycles, m.c);

        if(!sameCycle(index_a,index_c)) return NOT_APPLICABLE;

        edge_type ab_edge = check_edge(cycles, index_a, m.b);
        edge_type cd_edge = check_edge(cycles, index_c, m.d);

        if (ab_edge==no_edge || cd_edge==no_edge) return NOT_APPLICABLE;
        if (ab_edge==direct&&cd_edge==direct) return APPLICABLE;
        if (ab_edge==reversed&&cd_edge==reversed) return NEEDS_REVERSAL;
        return FUTURE_MAYBE;
    }

    Applicability isVertexMoveApplicable(const Move& m, const int* cycles) {
        int index_b = indexOf(cycles, m.b);
        int index_e = indexOf(cycles, m.e);

        if(sameCycle(index_b,index_e)) return NOT_APPLICABLE;

        int index_a = indexOf(cycles, m.a);
        int index_d = indexOf(cycles, m.d);

        edge_type ab_exists = check_edge(cycles, index_a, m.b);
        edge_type bc_exists = check_edge(cycles, index_b, m.c);

        edge_type de_exists = check_edge(cycles, index_d, m.e);
        edge_type ef_exists = check_edge(cycles, index_e, m.f);

        if (ab_exists==no_edge || bc_exists==no_edge || de_exists==no_edge || ef_exists==no_edge) return NOT_APPLICABLE;
       // TODO
        if (((ab_exists==direct && bc_exists==direct)||(ab_exists==reversed && bc_exists==reversed)) && ((de_exists==direct && ef_exists==direct)||(de_exists==reversed && ef_exists==reversed))) return APPLICABLE;

        return FUTURE_MAYBE;
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
                delta += instance.dist[this->cycles[i-1]][this->cycles[i]];
                delta += instance.dist[this->cycles[i-1+CYCLE_SIZE]][this->cycles[i+CYCLE_SIZE]];
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

    int evalEdgeSwapValues(const Instance& instance, const int cycles[NODE_COUNT], int a, int b, int c, int d) {
        (void)cycles;
        int delta = 0;
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

    // tylko dla międzycyklowych
    int evalVertexSwapValue(const Instance& instance, const int cycles[NODE_COUNT], int a, int b, int c, int d, int e, int f) {
        int delta = 0;

        delta += instance.dist[b][f];
        delta += instance.dist[a][e];
        delta += instance.dist[e][c];
        delta += instance.dist[d][b];
        delta -= instance.dist[a][b];
        delta -= instance.dist[b][c];
        delta -= instance.dist[d][e];
        delta -= instance.dist[e][f];
        

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

    Neighbourhood generate_candidates(const Instance& instance) {
        Neighbourhood c = {};
        for (int i = 0; i < NODE_COUNT; i++) {
            for(int j = 0; j < NEIGHBOURHOOD_SIZE; j++) {
                c[i][j] = -1;
            }
        }

        auto put_in_correct_spot = [&c, &instance](int first, int second){
            for(int i = 0; i < NEIGHBOURHOOD_SIZE; i++) {
                if (c[first][i] == -1) {
                    c[first][i] = second;
                    return;
                }
                int new_dist = instance.dist[first][second];
                int biggest_of_smallest = instance.dist[first][c[first][i]];
                if(new_dist < biggest_of_smallest) {
                    std::shift_right(c[first].begin()+i, c[first].end(), 1);
                    c[first][i] = second;
                    return;
                }
            }
        };
        for (int i = 0; i < NODE_COUNT; i++) {
            for (int j = 0; j < NODE_COUNT; j++) {
                if (i==j) continue;
                if (c[i][NEIGHBOURHOOD_SIZE-1] == -1) {
                    put_in_correct_spot(i, j);
                    continue;
                }
                if (instance.dist[i][j] < instance.dist[i][c[i][NEIGHBOURHOOD_SIZE-1]]) {
                    put_in_correct_spot(i, j);
                }
            }
        }

        return c;
    }

    void candidates_alg(const Instance& instance, SolutionOut* outSolution) {
        Solution sol(outSolution);
        sol.eval(instance);

        Neighbourhood candidates = generate_candidates(instance);
        static std::random_device rd;
        static std::mt19937 g(rd());
        bool improved;
        do {
            improved = false;
            Move bestMove = {0, false, -1, -1, -1, -1, -1, -1};
            std::array<int, NODE_COUNT> indices;
            std::iota(indices.begin(), indices.end(), 0);
            std::shuffle(indices.begin(), indices.end(), g);
           
            // raczej powinno być, że jeśli są w jednym cyklu to zamieniamy krawędzie, a nie wierzchołki
            for (int i = 0; i < NODE_COUNT; i++) {
                for (int j = 0; j < NEIGHBOURHOOD_SIZE; j++) {
                    int f = indices[i];
                    int s = indexOf(sol.cycles, candidates[f][j]);
                    int fp = absMod(f-1, CYCLE_SIZE) + (f>=CYCLE_SIZE?CYCLE_SIZE:0);
                    int fn = absMod(f+1, CYCLE_SIZE) + (f>=CYCLE_SIZE?CYCLE_SIZE:0);

                    if (!sameCycle(fp, s)) {
                        int delta = evalVertexSwap(instance, sol.cycles, fp, s);
                        if (delta < 0 && delta < bestMove.delta) {
                            bestMove = {delta, false, fp, s, -1, -1, -1, -1};
                            improved = true;
                        }
                    }

                    if (!sameCycle(fn, s)) {
                        int delta = evalVertexSwap(instance, sol.cycles, fn, s);
                        if (delta < 0 && delta < bestMove.delta) {
                            bestMove = {delta, false, fn, s, -1, -1, -1, -1};
                            improved = true;
                        }
                    }

                    if(f==s)continue;
                    if (sameCycle(f, s)) {
                        int delta = evalEdgeSwap(instance, sol.cycles, f, s);
                        if (delta < 0 && delta < bestMove.delta) {
                            bestMove = {delta, true, f, s, -1, -1, -1, -1};
                            improved = true;
                        }
                    }
                    
                }
            }
            
    
            // Jeśli znaleziono lepsze rozwiązanie, przejdź do niego
            if (improved) {
                if(bestMove.isEdgeMove) {
                    swapEdge(sol.cycles, bestMove.a, bestMove.b);
                } else {
                    swapVertex(sol.cycles, bestMove.a, bestMove.b);
                }

                sol.value += bestMove.delta;
            }

        } while (improved);
    
        sol.backToOut(outSolution);
    }



    void cache_alg(const Instance& instance, SolutionOut* outSolution) {
        Solution sol(outSolution);
        sol.eval(instance);
        MoveContainer LM;

        auto initializeMoveList = [&](MoveContainer& LM) {

            for (int c = 0; c < 2; ++c) {
                for (int i = 0; i < CYCLE_SIZE; ++i) {
                    for (int j = 0; j < CYCLE_SIZE; ++j) {
                        if(i == j ) continue;
                        if (absMod(i + 1, CYCLE_SIZE) == j || absMod(j + 1, CYCLE_SIZE) == i) continue;
                        int idx1 = i + c * CYCLE_SIZE;
                        int idx2 = j + c * CYCLE_SIZE;
                        int idx2succ = absMod(idx2 + 1, CYCLE_SIZE) + (c * CYCLE_SIZE);
                        int a = sol.cycles[idx1];
                        int b = sol.cycles[absMod(idx1 + 1, CYCLE_SIZE) + (c * CYCLE_SIZE)];
                        int c1 = sol.cycles[idx2];
                        int d = sol.cycles[idx2succ];
                        
                        int delta = evalEdgeSwap(instance, sol.cycles, idx1, idx2);
                        if (delta < 0) {
                            add_move(LM, {delta, true, a, b, c1, d, -1, -1});
                        }
                        delta = evalEdgeSwapValues(instance, sol.cycles, a, b, d, c1);
                        if (delta < 0) {
                            add_move(LM, {delta, true, a, b, d, c1, -1, -1});
                        }
                        delta = evalEdgeSwapValues(instance, sol.cycles, b, a, c1, d);
                        if (delta < 0) {
                            add_move(LM, {delta, true, b, a, c1, d, -1, -1});
                        }
                    }
                }
            }

            for (int i = 0; i < CYCLE_SIZE; ++i) {
                for (int j = 0; j < CYCLE_SIZE; ++j) {
                    int i1 = i;
                    int i2 = j + CYCLE_SIZE;

                    int temp1n = absMod(i1+1, CYCLE_SIZE);
                    int temp2n = absMod(i2+1, CYCLE_SIZE);
                    int temp1p = absMod(i1-1, CYCLE_SIZE);
                    int temp2p = absMod(i2-1, CYCLE_SIZE);
                    int i1n = i1 < CYCLE_SIZE ? temp1n : temp1n + CYCLE_SIZE;
                    int i2n = i2 < CYCLE_SIZE ? temp2n : temp2n + CYCLE_SIZE;
                    int i1p = i1 < CYCLE_SIZE ? temp1p : temp1p + CYCLE_SIZE;
                    int i2p = i2 < CYCLE_SIZE ? temp2p : temp2p + CYCLE_SIZE;

                    int a = sol.cycles[i1p];
                    int b = sol.cycles[i1];
                    int c = sol.cycles[i1n];
                    int d = sol.cycles[i2p];
                    int e = sol.cycles[i2];
                    int f = sol.cycles[i2n];
                    int delta = evalVertexSwap(instance, sol.cycles, i, j + CYCLE_SIZE);
                    if (delta < 0) {
                        add_move(LM, {delta, false, a, b, c, d, e, f});
                    }
                }
            }
            #ifdef CONTAINER_LIST
            LM.sort();
            #endif
        };

        auto updateMoveList = [&](MoveContainer& LM, const Move& performedMove) {
            // po zamianie krawędzi
            if (performedMove.isEdgeMove) {

                for (const auto& val : {performedMove.a, performedMove.b}) {
                    int cycle = indexOf(sol.cycles, val) > CYCLE_SIZE ? 1 : 0;
                    for (int j = 0; j < CYCLE_SIZE; ++j) {
                        int idx1 = indexOf(sol.cycles, val);
                        int idx1succ = absMod(idx1 + 1, CYCLE_SIZE) + (cycle * CYCLE_SIZE);
                        int idx2 = j + cycle * CYCLE_SIZE;
                        if (idx2==idx1 || idx2==idx1succ) continue;
                        int idx2succ = absMod(idx2 + 1, CYCLE_SIZE) + (cycle * CYCLE_SIZE);
                        int a = sol.cycles[idx1];
                        int b = sol.cycles[idx1succ];
                        int c1 = sol.cycles[idx2];
                        int d = sol.cycles[idx2succ];
                        int delta = evalEdgeSwap(instance, sol.cycles, idx1, idx2);
                        if (delta < 0) {
                            add_move(LM, {delta, true, a, b, c1, d, -1, -1});
                        }
                        delta = evalEdgeSwapValues(instance, sol.cycles, a, b, d, c1);
                        if (delta < 0) {
                            add_move(LM, {delta, true, a, b, d, c1, -1, -1});
                        }
                        delta = evalEdgeSwapValues(instance, sol.cycles, b, a, c1, d);
                        if (delta < 0) {
                            add_move(LM, {delta, true, b, a, c1, d, -1, -1});
                        }
                    }
                }

                for (const auto& val : {performedMove.a, performedMove.b, performedMove.c, performedMove.d}) {
                    int cycle;
                    if (indexOf(sol.cycles, val) < CYCLE_SIZE) {
                        cycle = 1;
                    } else {
                        cycle = 0;
                    }
                    for (int i = 0; i < CYCLE_SIZE; ++i) {
                        int idx1 = indexOf(sol.cycles, val);
                        int idx2 = i+CYCLE_SIZE*cycle;
                        int i1 = idx1;
                        int i2 = idx2;

                        int temp1n = absMod(i1+1, CYCLE_SIZE);
                        int temp2n = absMod(i2+1, CYCLE_SIZE);
                        int temp1p = absMod(i1-1, CYCLE_SIZE);
                        int temp2p = absMod(i2-1, CYCLE_SIZE);
                        int i1n = i1 < CYCLE_SIZE ? temp1n : temp1n + CYCLE_SIZE;
                        int i2n = i2 < CYCLE_SIZE ? temp2n : temp2n + CYCLE_SIZE;
                        int i1p = i1 < CYCLE_SIZE ? temp1p : temp1p + CYCLE_SIZE;
                        int i2p = i2 < CYCLE_SIZE ? temp2p : temp2p + CYCLE_SIZE;

                        int a = sol.cycles[i1p];
                        int b = sol.cycles[i1];
                        int c = sol.cycles[i1n];
                        int d = sol.cycles[i2p];
                        int e = sol.cycles[i2];
                        int f = sol.cycles[i2n];

                        int delta = evalVertexSwap(instance, sol.cycles, idx1, idx2);
                        if (delta < 0) {
                            add_move(LM, {delta, false, a, b, c, d, e, f});
                        }
                    }
                }
            } else {
            // po zamianie wierzchołków

                for (const auto& val : {performedMove.a, performedMove.e, performedMove.b, performedMove.f}) {
                    int cycle = indexOf(sol.cycles, val) > CYCLE_SIZE ? 1 : 0;
                    for (int j = 0; j < CYCLE_SIZE; ++j) {
                        int idx1 = indexOf(sol.cycles, val);
                        int idx1succ = absMod(idx1 + 1, CYCLE_SIZE) + (cycle * CYCLE_SIZE);
                        int idx2 = j + cycle * CYCLE_SIZE;
                        if (idx2==idx1 || idx2==idx1succ) continue;
                        int idx2succ = absMod(idx2 + 1, CYCLE_SIZE) + (cycle * CYCLE_SIZE);
                        int a = sol.cycles[idx1];
                        int b = sol.cycles[idx1succ];
                        int c1 = sol.cycles[idx2];
                        int d = sol.cycles[idx2succ];
                        int delta = evalEdgeSwap(instance, sol.cycles, idx1, idx2);
                        if (delta < 0) {
                            add_move(LM, {delta, true, a, b, c1, d, -1, -1});
                        }
                        delta = evalEdgeSwapValues(instance, sol.cycles, a, b, d, c1);
                        if (delta < 0) {
                            add_move(LM, {delta, true, a, b, d, c1, -1, -1});
                        }
                        delta = evalEdgeSwapValues(instance, sol.cycles, b, a, c1, d);
                        if (delta < 0) {
                            add_move(LM, {delta, true, b, a, c1, d, -1, -1});
                        }
                    }
                }

                for (const auto& val : {performedMove.a, performedMove.b, performedMove.c, performedMove.d, performedMove.e, performedMove.f}) {
                    int cycle;
                    if (indexOf(sol.cycles, val) < CYCLE_SIZE) {
                        cycle = 1;
                    } else {
                        cycle = 0;
                    }
                    for (int i = 0; i < CYCLE_SIZE; ++i) {
                        int idx1 = indexOf(sol.cycles, val);
                        int idx2 = i+CYCLE_SIZE*cycle;
                        int i1 = idx1;
                        int i2 = idx2;

                        int temp1n = absMod(i1+1, CYCLE_SIZE);
                        int temp2n = absMod(i2+1, CYCLE_SIZE);
                        int temp1p = absMod(i1-1, CYCLE_SIZE);
                        int temp2p = absMod(i2-1, CYCLE_SIZE);
                        int i1n = i1 < CYCLE_SIZE ? temp1n : temp1n + CYCLE_SIZE;
                        int i2n = i2 < CYCLE_SIZE ? temp2n : temp2n + CYCLE_SIZE;
                        int i1p = i1 < CYCLE_SIZE ? temp1p : temp1p + CYCLE_SIZE;
                        int i2p = i2 < CYCLE_SIZE ? temp2p : temp2p + CYCLE_SIZE;

                        int a = sol.cycles[i1p];
                        int b = sol.cycles[i1];
                        int c = sol.cycles[i1n];
                        int d = sol.cycles[i2p];
                        int e = sol.cycles[i2];
                        int f = sol.cycles[i2n];

                        int delta = evalVertexSwap(instance, sol.cycles, idx1, idx2);
                        if (delta < 0) {
                            add_move(LM, {delta, false, a, b, c, d, e, f});
                        }
                    }
                }
            }

            #ifdef CONTAINER_LIST
            LM.sort();
            #endif
            
        };


        initializeMoveList(LM);


        bool improved;
        do {
EOI:
            improved = false;

            for (auto it = LM.begin(); it != LM.end(); ) {
                Applicability app = APPLICABLE;

                if (it->isEdgeMove) {
                    app = isEdgeMoveApplicable(*it, sol.cycles);
                } else {
                    app = isVertexMoveApplicable(*it, sol.cycles);
                }

                int a = it->a;
                int b = it->b;
                int c = it->c;
                int d = it->d;
                int e = it->e;
                int f = it->f;
                int delta = it->delta;

                if(app==NEEDS_REVERSAL){
                    if (it->isEdgeMove){
                        a = it->b;
                        b = it->a;
                        c = it->d;
                        d = it->c;
                    } else {
                        a = it->c;
                        c = it->a;
                        d = it->f;
                        f = it->d;
                    }
                    app = APPLICABLE;
                }

                switch (app) {
                    case APPLICABLE:

                        if (it->isEdgeMove) {
                            swapEdge(sol.cycles, indexOf(sol.cycles, a), indexOf(sol.cycles, c));
                        } else {
                            swapVertex(sol.cycles, indexOf(sol.cycles, b), indexOf(sol.cycles, e));
                        }
                        sol.value += delta;

                        {
                            Move performedMove = {delta,it->isEdgeMove,a,b,c,d,e,f};
                            LM.erase(it);
                            updateMoveList(LM, performedMove);
                        }
                        improved = true;

                        goto EOI;

                    case FUTURE_MAYBE:
                        // Jeszcze nie teraz - nie usuwamy ruchu
                        ++it;
                        break;

                    case NOT_APPLICABLE:
                        it = LM.erase(it);
                        break;
                }
            }
        } while (improved);
        sol.backToOut(outSolution);
    }

    int main() {

        SolutionOut solOut = {{2,4,0,1,3,5},{6,7,8,9,10,11},0};
        Instance instance =  {{
            {0, 2, 4, 7, 6, 4, 12, 22, 5, 32, 54, 6},
            {2, 0, 2, 5, 5, 2, 10, 20, 4, 30, 51, 4},
            {4, 2, 0, 4, 7, 3, 9, 18, 5, 28, 50, 3},
            {7, 5, 4, 0, 7, 4, 5, 15, 6, 25, 46, 2},
            {6, 5, 7, 7, 0, 4, 10, 20, 1, 30, 52, 5},
            {4, 2, 3, 4, 4, 0, 8, 18, 2, 28, 50, 2},
            {12, 10, 9, 5, 10, 8, 0, 10, 9, 20, 42, 6},
            {22, 20, 18, 15, 20, 18, 10, 0, 19, 10, 32, 16},
            {5, 4, 5, 6, 1, 2, 9, 19, 0, 29, 51, 4},
            {32, 30, 28, 25, 30, 28, 20, 10, 29, 0, 22, 26},
            {54, 51, 50, 46, 52, 50, 42, 32, 51, 22, 0, 48},
            {6, 4, 3, 2, 5, 2, 6, 16, 4, 26, 48, 0}
        }};

        cache_alg(instance, &solOut);
    }

}