#include "Graph.h"

/*******************************************************************************
 * Queue
 ******************************************************************************/
void Queue::display() {
    node *p = new node;
    p = front;
    if (front == NULL) {
        cout << "\nNothing to Display\n";
    } else {
        while (p != NULL) {
            cout << endl << p->info;
            p = p->next;
        }
    }
}

Queue::Queue() {
    front = NULL;
    rear = NULL;
}

Queue::~Queue() {
    delete front;
}

void Queue::enqueue(int data) {
    node *temp = new node();
    temp->info = data;
    temp->next = NULL;
    if (front == NULL) {
        front = temp;
    } else {
        rear->next = temp;
    }
    rear = temp;
}

int Queue::dequeue() {
    node *temp = new node();
    int value;
    if (front == NULL) {
        cout << "\nQueue is Emtpty\n";
    } else {
        temp = front;
        value = temp->info;
        front = front->next;
        delete temp;
    }
    return value;
}

bool Queue::isEmpty() {
    return (front == NULL);
}

/*******************************************************************************
 * Edge
 ******************************************************************************/
Edge::Edge(int node_1, int node_2, int label_) {
    node1 = node_1;
    node2 = node_2;
    label = label_;
}

/*******************************************************************************
 * Graph
 ******************************************************************************/

Graph::Graph(int size) {
    int i, j;

    m = 0;

    if (size < 2)
        n = 2;
    else
        n = size;

    A = new int*[n];

    for (i = 0; i < n; ++i)
        A[i] = new int[n];

    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            A[i][j] = -1;
}

Graph::~Graph() {
    for (int i = 0; i < n; ++i)
        delete [] A[i];
    delete [] A;
    delete [] Visited;
    delete [] Node_labels;
    delete [] Edge_labels;
}

int Graph::nodes_size() {
    return n;
}

int Graph::edges_size() {
    return m;
}

/******************************************************
Checks if two given vertices are connected by an edge
@param u vertex
@param v vertex
@return true if connected false if not connected
 ******************************************************/
bool Graph::isConnected(int u, int v, int &val) {
    val = A[u][v];
    return (val > -1);
}

/*****************************************************
adds an edge E to the graph G.
@param u vertex
@param v vertex
 *****************************************************/
void Graph::addEdge(int u, int v) {
    A[u][v] = A[v][u] = m; //Always bigger than -1 and stores the edge number   
    edges.push_back(Edge(u, v, m));
    m++;
}

/*****************************************************
performs Breadth First Search
@param s initial vertex
 *****************************************************/
void Graph::BFS(int s) {
    int d = 0; //node discovery index
    int h = 0; //edge discovery index
    int edge;
    int w;
    Queue Q;

    /** Keeps track of explored vertices */
    bool *explored = new bool[n];
    Visited = new int[n];
    Node_labels = new int[n];
    Edge_labels = new int[m];

    /** Initailized all vertices as unexplored */
    for (int i = 1; i < n; ++i) {
        explored[i] = false;
        Visited[i] = 0;
    }

    /** Push initial vertex to the queue */
    Q.enqueue(s);
    explored[s] = true; /** mark it as explored */
    Visited[s]++;
    cout << "Breadth first Search starting from vertex ";
    cout << s << " : " << endl;

    /** Unless the queue is empty */
    while (!Q.isEmpty()) {
        /** Pop the vertex from the queue */
        int v = Q.dequeue();

        /** display the explored vertices */
        cout << v << " ";
        Node_labels[d] = v;


        /** From the explored vertex v try to explore all the
        connected vertices */
        //for (int w = 0; w < n; w++) {
        for (Edge e : edges) {
            if (e.node1 == v || e.node2 == v){
                
                if (e.node1 ==v)
                    w = e.node2;
                else
                    w = e.node1;
                
                /** Explores the vertex w if it is connected to v
                and and if it is unexplored */
                if (Visited[w] == 0) {
                    /** adds the vertex w to the queue */
                    Q.enqueue(w);
                    /** marks the vertex w as visited */
                    explored[w] = true;

                    Edge_labels[h] = e.label;
                    h++;
                    Visited[w]++;
                }
            }
        }
        d++;
    }
    cout << endl;
    delete [] explored;
}