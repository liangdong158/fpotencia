/****************************************************************
 * Written by Bibek Subedi
 * http://www.programming-techniques.com/2012/07/breadth-first-search-in-c-algorithm-and.html
 * 
 * Modifyed by Santiago Peñate Vera to be correct and detect loops
 ****************************************************************/

#include "fpotencia_libs.h"

using namespace std;

/****************************************************************
 
Class Queue represent a Queue data structure which is First In
First Out [FIFO] structured. It has operations like Enqueue which
adds an element at the rear side and Dequeue which removes the
element from front.
 
 *****************************************************************/

struct node {
    int info;
    node *next;
};

class Queue {
private:

    node *front;

    node *rear;

public:
    Queue();

    ~Queue();

    bool isEmpty();

    void enqueue(int);

    int dequeue();

    void display();

};

/************************************************************
 
Class Graph represents a Graph [V,E] having vertices V and
edges E.
 
 ************************************************************/
class Edge {
public:
    Edge(int node_1, int node_2, int label_);
    int label;
    int node1;
    int node2;
};

/************************************************************
 
Class Graph represents a Graph [V,E] having vertices V and
edges E.
 
 ************************************************************/
class Graph {
private:

    // n is the number of vertices in the graph
    int n;

    //number of edges in the graph
    int m;

    // A stores the edges between two vertices
    int **A;

    vector<Edge> edges;

public:
    Graph(int size = 2);

    ~Graph();

    bool isConnected(int u, int v, int &val);

    void addEdge(int u, int v);

    void BFS(int);

    int nodes_size();

    int edges_size();

    int *Visited;

    int *Node_labels;

    int *Edge_labels;
};