#include <iostream>
#include "util.h"
#include "node.h"

using namespace std;

void printV(const int *a, int n, const char *comment) {
    cout << comment << ": n= " << n << " : ";
    for (int i = 0; i < Min(n, 20); i++) {
        cout << a[i] << " ";
    }
    cout << endl;
}

void print_node(int i, const Node& t) {
    printf("%2d: %d [%2d,%2d,%2d];", i, t._val, 
          t._parent, t._left, t._right);

}

void print_tree(const Node *tree, int n, const char *comment) {
    cout << comment << ": n= " << n << " : ";
    for (int i = 0; i < Min(n, 20); i++) {
        print_node(i, tree[i]);
        cout  << endl;
    }
}

void __assert(const char *file, int line) {
    cerr << "\nAssertion violation " << __FILE__ << ":" << __LINE__ << endl;
}

int **new2(int n, int m) {
    int **v = new int *[n];
    for (int i = 0; i < n; i++) {
        v[i] = new int[m];
    }
    return v;
}

void delete2(int **v, int n) {
    for (int i = 0; i < n; i++) {
        int *x = v[i];
        delete[] x;
    }
    delete[] v;
}
