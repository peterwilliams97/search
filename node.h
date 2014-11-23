#ifndef NODE_H
#define NODE_H

struct Node {
    int _val;
    int _left;
    int _right;
    int _parent;
    Node() : _val(0), _left(-1), _right(-1), _parent(-1) {}
};

void print_node(int i, const Node& node); 
void print_tree(const Node *tree, int n, const char *comment); 

#endif // #ifndef NODE_H
