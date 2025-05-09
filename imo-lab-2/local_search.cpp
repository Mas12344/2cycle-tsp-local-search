#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <time.h>

// #define NODE_COUNT = 4

extern "C" {
    const int NODE_COUNT = 200;
    const int CYCLE_SIZE = 100;

    struct Instance { 
        int dist[NODE_COUNT][NODE_COUNT];
    };

    struct SolutionOut {
        int cycle1[CYCLE_SIZE];
        int cycle2[CYCLE_SIZE];
        int value;
    }; 

    struct Move {
        int cA;
        int cB;
        int nA;
        int nB;
    };

    struct EdgeMove {
        int c;
        int nA;
        int nB;
    };
    
    int absMod(int val, int mod) {
        if(val<0) val+=mod;
        return val%mod;
    }

    class Solution {
    public:
        std::vector<std::vector<int>> cycles;
        int value;

        void evaluate(const Instance& instance) {
            int val = 0;
            for(auto cycle : this->cycles) {
                for(int i=1; i<cycle.size(); i++) 
                    val += instance.dist[cycle[i-1]][cycle[i]];
                val += instance.dist[cycle[CYCLE_SIZE-1]][cycle[0]];
            }
            this->value =val;
        }

        int evaluateMove(const Instance& instance, Move move) {
            int val = 0;
            int nodeA = this->cycles[move.cA][move.nA];
            int nodeB = this->cycles[move.cB][move.nB];
            int prevA = this->cycles[move.cA][absMod(move.nA-1, CYCLE_SIZE)];
            int prevB = this->cycles[move.cB][absMod(move.nB-1, CYCLE_SIZE)];
            int nextA = this->cycles[move.cA][absMod(move.nA+1, CYCLE_SIZE)];
            int nextB = this->cycles[move.cB][absMod(move.nB+1, CYCLE_SIZE)];
        
            if(move.cA==move.cB && (nodeA==prevB || nodeB==prevA)) {
                if (nodeA==prevB) {
                    val -= instance.dist[prevA][nodeA];
                    val -= instance.dist[nodeB][nextB];
        
                    val += instance.dist[prevA][nodeB];
                    val += instance.dist[nodeA][nextB];
                } else {
                    val -= instance.dist[prevB][nodeB];
                    val -= instance.dist[nodeA][nextA];
        
                    val += instance.dist[prevB][nodeA];
                    val += instance.dist[nodeB][nextA];
                }
            } else {
                val -= instance.dist[prevA][nodeA];
                val -= instance.dist[prevB][nodeB];
                val -= instance.dist[nodeA][nextA];
                val -= instance.dist[nodeB][nextB];
                
                val += instance.dist[prevB][nodeA];
                val += instance.dist[prevA][nodeB];
                val += instance.dist[nodeB][nextA];
                val += instance.dist[nodeA][nextB];
            }
            return val;
        }

        int evaluateEdgeMove(const Instance& instance, EdgeMove move) {
            int val = 0;
            int nodeA = this->cycles[move.c][move.nA];
            int nodeB = this->cycles[move.c][move.nB];

            int nextA = this->cycles[move.c][absMod(move.nA+1, CYCLE_SIZE)];
            int nextB = this->cycles[move.c][absMod(move.nB+1, CYCLE_SIZE)];

            val -= instance.dist[nodeA][nextA];
            val -= instance.dist[nodeB][nextB];

            val += instance.dist[nodeA][nodeB];
            val += instance.dist[nextA][nextB];

            return val;
        }



        Solution() {};
        Solution(Solution& other) {
            for(int i=0; i<2; i++) {
                std::vector<int> copy = other.cycles[i];
                this->cycles.push_back(copy);
            }
            this->value = other.value;
        }
        Solution(SolutionOut* other) {
            this->cycles.resize(2);
            for(int i=0; i<2; i++) this->cycles[i].reserve(CYCLE_SIZE);
            for(int i=0; i<CYCLE_SIZE; i++) {
                this->cycles[0].push_back(other->cycle1[i]);
            }
            for(int i=0; i<CYCLE_SIZE; i++) {
                this->cycles[1].push_back(other->cycle2[i]);
            }
        }

        void convertToOutput(SolutionOut* out) {
            for(int i=0; i<cycles[0].size(); i++) {
                out->cycle1[i] = this->cycles[0][i];
            }
            for(int i=0; i<cycles[1].size(); i++) {
                out->cycle2[i] = this->cycles[1][i];
            }
            out->value = this->value;
        }

        void printCycle() {
            std::cout<<cycles[0][0];
            for(int i=1; i<cycles[0].size();i++) {
                std::cout<<"->"<<cycles[0][i];
            }
            std::cout<<std::endl<<cycles[1][0];;
            for(int i=1; i<cycles[1].size();i++) {
                std::cout<<"->"<<cycles[1][i];
            }
            std::cout<<std::endl;
        }
    };

    void randomSearch(Instance instance, SolutionOut* outSolution, int iterCount) {
        std::srand((unsigned)time(nullptr));
        Solution sol(outSolution);
        sol.evaluate(instance);
        Solution bestSol(sol);

        for(int i=0; i<1000000; i++) {
            int cA = std::rand() % 2;
            int cB = std::rand() % 2;

            int nodeA = std::rand() % CYCLE_SIZE;
            int nodeB = std::rand() % CYCLE_SIZE;
            if(cA == cB && nodeA == nodeB) nodeB = (nodeB+1) % CYCLE_SIZE;
            std::swap(sol.cycles[cA][nodeA], sol.cycles[cB][nodeB]);
            sol.evaluate(instance);

            if(sol.value < bestSol.value) {
                bestSol = Solution(sol);
            }
        }

        bestSol.convertToOutput(outSolution);
    }

    void greedySearch(Instance instance, SolutionOut* outSolution) {
        std::srand((unsigned)time(nullptr));
        Solution sol(outSolution);
        sol.evaluate(instance);
    
        // Generowanie wszystkich możliwych zamian
        std::vector<Move> moves;
        for (int cA = 0; cA < 2; cA++) {
            for (int cB = 0; cB < 2; cB++) {
                for (int i = 0; i < CYCLE_SIZE; i++) {
                    for (int j = 0; j < CYCLE_SIZE; j++) {
                        if (cA == cB && i >= j) continue;
                        moves.push_back({cA, cB, i, j});
                    }
                }
            }
        }
    
        bool improved;
        do {
            improved = false;
    
            int dir = std::rand() % 2;
    
            if (dir == 0) {
                for (size_t i = 0; i < moves.size(); i++) {
                    Move m = moves[i];
                    int delta = sol.evaluateMove(instance, m);
                    if (delta < 0) {
                        std::swap(sol.cycles[m.cA][m.nA], sol.cycles[m.cB][m.nB]);
                        sol.value += delta;
                        improved = true;
                        break;
                    }
                }
            } else {
                for (size_t i = moves.size(); i-- > 0;) {
                    Move m = moves[i];
                    int delta = sol.evaluateMove(instance, m);
                    if (delta < 0) {
                        std::swap(sol.cycles[m.cA][m.nA], sol.cycles[m.cB][m.nB]);
                        sol.value += delta;
                        improved = true;
                        break;
                    }
                }
            }
    
        } while (improved);
    
        sol.convertToOutput(outSolution);
    }

    void greedySearchRnd(Instance instance, SolutionOut* outSolution) {
        std::srand((unsigned)time(nullptr));
        Solution sol(outSolution);
        sol.evaluate(instance);

        // generowanie wszystkich możliwych zamian
        std::vector<Move> moves;
        for (int cA = 0; cA < 2; cA++) {
            for (int cB = 0; cB < 2; cB++) {
                for (int i = 0; i < CYCLE_SIZE; i++) {
                    for (int j = 0; j < CYCLE_SIZE; j++) {
                        if (cA == cB && i >= j) continue; //chyba?
                        moves.push_back({cA, cB, i, j});
                    }
                }
            }
        }

        bool improved;
        do {
            improved = false;
            std::shuffle(moves.begin(), moves.end(), std::default_random_engine(std::rand()));

            for (auto m : moves) {
                int delta = sol.evaluateMove(instance, m);
                
                if (delta < 0) {
                    std::swap(sol.cycles[m.cA][m.nA], sol.cycles[m.cB][m.nB]);

                    sol.value += delta;
                    improved = true;
                    break;
                }
            }
            //break;

        } while (improved);

        sol.convertToOutput(outSolution);
    }

    void steepestSearch(Instance instance, SolutionOut* outSolution) {
        std::srand((unsigned)time(nullptr));
        Solution sol(outSolution);
        sol.evaluate(instance);
        
        bool improved;
        do {
            improved = false;
            Move bestMove;
            int bestDelta = 0;
    
            // Przegląd wszystkich możliwych ruchów (zamian)
            for (int cA = 0; cA < 2; cA++) {
                for (int cB = 0; cB < 2; cB++) {
                    for (int i = 0; i < CYCLE_SIZE; i++) {
                        for (int j = 0; j < CYCLE_SIZE; j++) {
                            if (cA == cB && i >= j) continue;
                            Move move = {cA,cB,i,j};
                            int delta = sol.evaluateMove(instance, move);
    
                            if (delta < 0 && delta < bestDelta) {
                                bestDelta = delta;
                                bestMove = move;
                                improved = true;
                            }
                        }
                    }
                }
            }
    
            // Jeśli znaleziono lepsze rozwiązanie, przejdź do niego
            if (improved) {
                std::swap(sol.cycles[bestMove.cA][bestMove.nA], sol.cycles[bestMove.cB][bestMove.nB]);
                sol.value += bestDelta;
            }
            
    
        } while (improved);
    
        sol.convertToOutput(outSolution);
    }

    void greedySearchEdgesRnd(Instance instance, SolutionOut* outSolution) {
        std::srand((unsigned)time(nullptr));
        Solution sol(outSolution);
        sol.evaluate(instance);

        // generowanie wszystkich możliwych zamian
        std::vector<EdgeMove> moves;
        for (int c = 0; c < 2; c++) {
            for (int i = 0; i < CYCLE_SIZE; i++) {
                for (int j = 0; j < CYCLE_SIZE; j++) {
                    if(i >= j ) continue;
                    if((i+1)%CYCLE_SIZE == j || (j+1)%CYCLE_SIZE == i) continue;
                    moves.push_back({c, i, j});
                }
            }
        }

        bool improved;
        do {
            improved = false;
            std::shuffle(moves.begin(), moves.end(), std::default_random_engine(std::rand()));

            for (auto m : moves) {
                int delta = sol.evaluateEdgeMove(instance, m);
                
                if (delta < 0) {
                    std::reverse(sol.cycles[m.c].begin() + m.nA + 1, sol.cycles[m.c].begin() + m.nB + 1);
                    sol.value += delta;
                    improved = true;
                    break;
                }
            }
            //break;

        } while (improved);

        sol.convertToOutput(outSolution);
    }

    void greedySearchEdges(Instance instance, SolutionOut* outSolution) {
        std::srand((unsigned)time(nullptr));
        Solution sol(outSolution);
        sol.evaluate(instance);
    
        // Generowanie możliwych zamian krawędzi
        std::vector<EdgeMove> edgeMoves;
        for (int c = 0; c < 2; c++) {
            for (int i = 0; i < CYCLE_SIZE; i++) {
                for (int j = 0; j < CYCLE_SIZE; j++) {
                    if (i >= j) continue;
                    if ((i + 1) % CYCLE_SIZE == j || (j + 1) % CYCLE_SIZE == i) continue;
                    edgeMoves.push_back({c, i, j});
                }
            }
        }
        std::vector<int> edgeIndices(edgeMoves.size());
        std::iota(edgeIndices.begin(), edgeIndices.end(), 0);
        std::vector<int> iIndices(CYCLE_SIZE), jIndices(CYCLE_SIZE);
        std::iota(iIndices.begin(), iIndices.end(), 0);
        std::iota(jIndices.begin(), jIndices.end(), 0);
    
        bool improved;
        do {
            improved = false;
    
            // Losowa kolejność przeglądania edgeMoves
            if (std::rand() % 2) std::reverse(edgeIndices.begin(), edgeIndices.end());
    
            for (int idx : edgeIndices) {
                EdgeMove m = edgeMoves[idx];
                int delta = sol.evaluateEdgeMove(instance, m);
                if (delta < 0) {
                    std::reverse(sol.cycles[m.c].begin() + m.nA + 1, sol.cycles[m.c].begin() + m.nB + 1);
                    sol.value += delta;
                    improved = true;
                    break;
                }
            }
    
            // Losowa kolejność przeglądania Move
            if (std::rand() % 2) std::reverse(iIndices.begin(), iIndices.end());
            if (std::rand() % 2) std::reverse(jIndices.begin(), jIndices.end());
    
            for (int i : iIndices) {
                for (int j : jIndices) {
                    Move m = {0, 1, i, j};  // swap vertex i from cycle 0 with j from cycle 1
                    int delta = sol.evaluateMove(instance, m);
                    if (delta < 0) {
                        std::swap(sol.cycles[m.cA][m.nA], sol.cycles[m.cB][m.nB]);
                        sol.value += delta;
                        improved = true;
                        goto nextIteration; // break both loops
                    }
                }
            }
    
            nextIteration:
            continue;
    
        } while (improved);
    
        sol.convertToOutput(outSolution);
    }

    void steepestSearchEdges(Instance instance, SolutionOut* outSolution) {
        std::srand((unsigned)time(nullptr));
        Solution sol(outSolution);
        sol.evaluate(instance);
        
        bool improved;
        do {
            improved = false;
            bool edgeMove = true;
            EdgeMove bestMoveEdge;
            Move bestMove;
            int bestDelta = 0;
    
            // Przegląd wszystkich możliwych ruchów (zamian)
            for (int c = 0; c < 2; c++) {
                for (int i = 0; i < CYCLE_SIZE; i++) {
                    for (int j = 0; j < CYCLE_SIZE; j++) {
                        if(i >= j ) continue;
                        if((i+1)%CYCLE_SIZE == j || (j+1)%CYCLE_SIZE == i) continue;
                        EdgeMove move = {c,i,j};
                        int delta = sol.evaluateEdgeMove(instance, move);

                        if (delta < 0 && delta < bestDelta) {
                            bestDelta = delta;
                            bestMoveEdge = move;
                            improved = true;
                        }
                    }
                }
            }

            for (int i = 0; i < CYCLE_SIZE; i++) {
                for (int j = 0; j < CYCLE_SIZE; j++) {
                    Move move = {0, 1, i, j};
                    int delta = sol.evaluateMove(instance, move);

                    if (delta < 0 && delta < bestDelta) {
                        edgeMove = false;
                        bestDelta = delta;
                        bestMove = move;
                        improved = true;
                    }
                }
            }

    
            // Jeśli znaleziono lepsze rozwiązanie, przejdź do niego
            if (improved) {
                if(edgeMove) {
                    std::reverse(sol.cycles[bestMoveEdge.c].begin() + bestMoveEdge.nA + 1, 
                                 sol.cycles[bestMoveEdge.c].begin() + bestMoveEdge.nB + 1);
                } else {
                    std::swap(sol.cycles[bestMove.cA][bestMove.nA], sol.cycles[bestMove.cB][bestMove.nB]);
                }
                sol.value += bestDelta;
            }

        } while (improved);
    
        sol.convertToOutput(outSolution);
    }

    // int main() {
    //     SolutionOut solOut = {{2,4,0,1,3,5},{6,7,8,9,10,11},0};
    //     Instance instance =  {{{0, 2, 4, 7, 6, 4, 12, 22, 5, 32, 54, 6},
    //     {2, 0, 2, 5, 5, 2, 10, 20, 4, 30, 51, 4},
    //     {4, 2, 0, 4, 7, 3, 9, 18, 5, 28, 50, 3},
    //     {7, 5, 4, 0, 7, 4, 5, 15, 6, 25, 46, 2},
    //     {6, 5, 7, 7, 0, 4, 10, 20, 1, 30, 52, 5},
    //     {4, 2, 3, 4, 4, 0, 8, 18, 2, 28, 50, 2},
    //     {12, 10, 9, 5, 10, 8, 0, 10, 9, 20, 42, 6},
    //     {22, 20, 18, 15, 20, 18, 10, 0, 19, 10, 32, 16},
    //     {5, 4, 5, 6, 1, 2, 9, 19, 0, 29, 51, 4},
    //     {32, 30, 28, 25, 30, 28, 20, 10, 29, 0, 22, 26},
    //     {54, 51, 50, 46, 52, 50, 42, 32, 51, 22, 0, 48},
    //     {6, 4, 3, 2, 5, 2, 6, 16, 4, 26, 48, 0}}};

    //     steepestSearchEdges(instance, &solOut);
    //     std::cout<<solOut.value<<std::endl;
    // }

}