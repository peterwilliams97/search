/*
 * http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=lowestCommonAncestor#Range_Minimum_Query_%28RMQ%29
 *
 */

#include <iostream>
#include "util.h"
using namespace std;

/*
 * Build naive RMQ look-up-table rmq_lut for sequence s
 * s[n]
 * rmq_lut[n][n]  rmq_lut[i][j] = rmq_lut(i, j)
 *
 * O(n^2) time and space construction
 * O(1) operation
 */
void preprocess_rmq_1(const int *s, int n, int **rmq_lut) {
    int i, j;
    for (i = 0; i < n; i++) {
        rmq_lut[i][i] = i;
    }
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            int j1 = rmq_lut[i][j - 1];
            rmq_lut[i][j] = s[j1] < s[j] ? j1 : j;
        }
    }
}

/*
 * Build RMQ look-up-table rmq_lut for sequence s
 * s[n]
 * rmq_lut[sqrt(n)]
 *
 * O(n) time, O(sqrt(n)) space construction
 * O(sqrt(n)) operation
 */
void preprocess_rmq_2(const int *s, int n, int *rmq_lut) {
    int i, j;
    int r = int(floor(sqrt(n)));
    int nr = (n + r - 1) / r;
    int min_v;

    for (i = 0, j = -1; i < n; i++) {
        int v = s[i];
        if (i % r == 0) {
            j++;
            rmq_lut[j] = i;
            min_v = v;
        } else if (v < min_v) {
            rmq_lut[j] = i;
            min_v = v;
        }
    }
}

int get_rmq_2(const int *s, int n, const int *rmq_lut, int i, int j) {
    int r = int(floor(sqrt(n)));
    int nr = (n + r - 1) / r;

    int min_k = -1;
    int min_v;
    int k;
    int v;

    //int ti = (i + r - 1) / r; // 1 offset
    int ti = (i + r) / r - 1; // 0 offset
    if (ti * r != i) {
        for (k = i; k < ti * r; k++) {
            v = s[k];
            if (min_k < 0 || v < min_v ) { min_k = k; min_v = v; }
        }
    }

    //int tj = j / r;
    int tj = (j + 1)/ r - 1;
    if (tj * r != j) {
        for (k = tj * r + 1; k <= j; k++) {
            v = s[k];
            if (min_k < 0 || v < min_v ) { min_k = k; min_v = v; }
        }
    }

    for (int t = ti; t <= tj ; t++) {
        k = rmq_lut[t];
        v = s[k];
        if (min_k < 0 || v < min_v ) { min_k = k; min_v = v; }
    }

    return min_k;
}

/*
 * Sparse Table (ST) algorithm
 * Build RMQ look-up-table rmq_lut for sequence s
 * s[n]
 * rmq_lut[n][log(n)]
 *
 * O(n log(n)) time and space construction
 * O(1) operation
 */
void preprocess_rmq_3(const int *s, int n, int **rmq_lut) {
    int i, j;

    // initialize rmq_lut for the intervals with length 1
    for (i = 0; i < n; i++) {
        rmq_lut[i][0] = i;
    }
    // compute values from smaller to bigger intervals
    for (j = 1; 1 << j <= n; j++) {
        //cout << j << endl; // !@#$
        for (i = 0; i + (1 << j) - 1 < n; i++) {
            int rmq1 = rmq_lut[i][j - 1];
            int rmq2 = rmq_lut[i + (1 << (j - 1))][j - 1];
            rmq_lut[i][j] = s[rmq1] < s[rmq2] ? rmq1 : rmq2;
        }
    }
}

int get_rmq_3(const int *s, int n, int * const*rmq_lut, int i, int j) {
    int k = int(log2(j - i + 1));
    int j1 = j - ipow(2, k) + 1;

    int rmq1 = rmq_lut[i][k];
    int rmq2 = rmq_lut[j1][k];
    return s[rmq1] < s[rmq2] ? rmq1 : rmq2;
}

#if 0
void preprocess_rmq_3(const int *s, int n, int **rmq_lut) {
    if (b == e)
        M[node] = b;
    else
    {
        //compute the values in the left and right subtrees
        initialize(2 * node, b, (b + e) / 2, M, A, N);
        initialize(2 * node + 1, (b + e) / 2 + 1, e, M, A, N);
        //search for the minimum value in the first and
        //second half of the interval
        if (A[M[2 * node]] <= A[M[2 * node + 1]])
            M[node] = M[2 * node];
        else
            M[node] = M[2 * node + 1];
    }
}
#endif

/*
* Build Cartesian tree
*  http://wcipeg.com/wiki/Cartesian_tree
* s[n]
* tree[n]
*
*/
void compute_tree(const int *s, int n, int *tree) {
    int *st = new int[n];
    int i, k, top = -1;

    //we start with an empty stack
    //at step i we insert A[i] in the stack
    for (i = 0; i < n; i++)  {
        //compute the position of the first element that is
        //equal or smaller than s[i]
        k = top;
        while (k >= 0 && s[st[k]] > s[i]) {
            k--;
        }
        //we modify the tree as explained above
        if (k != -1) {
            tree[i] = st[k];
            cout << "   tree[" << i << "]=" << tree[i] << ",v=" << s[tree[i]] << endl;
        }
        if (k < top) {
            tree[st[k + 1]] = i;
            cout << "   tree[" << st[k + 1] << "]=" << i << ",v=" << s[i] << endl;;
        }
        //we insert A[i] in the stack and remove
        //any bigger elements
        st[++k] = i;
        cout << "i=" << i << ", ";
        printV(st, k, "stack");
        top = k;

    }
    //the first element in the stack is the root of
    //the tree, so it has no father
    tree[st[0]] = -1;
    cout << "   tree[" << st[0] << "]=" << -1 << ",v=*" << endl;;


    delete[] st;
}

int
test_rmq_3(const int *s, int n, int i, int j, int rmq_ij) {
    int lg_n = int(ceil(log2(n)));
    //cout << "test_rmq_3: n = " << n << ", lg_n = " << lg_n << endl;
    int **rmq_lut = new2(n, lg_n + 1);

    int rmq_k = -1;
    preprocess_rmq_3(s, n, rmq_lut);
#if 1
    rmq_k = get_rmq_3(s, n, rmq_lut, i, j);
#endif
    delete2(rmq_lut, n);
    return rmq_k;
}

void
test_compute_tree(const int *s, int n) {
    int *tree = new int[n];

    compute_tree(s, n, tree);

    printV(tree, n, "tree");

    delete [] tree;
}

bool test_rmq() {
#if 1
    for (int k = 1; k < 1000; k++) {
        int **rmq_lut = new2(k, 10);
        delete2(rmq_lut, k);
    }
#endif

    bool ok = false;
    int n = 4000;
#if 1
    int *s = new int[n];
    for (int k = 0; k < n; k++) {
        s[k] = 100 + (k + n / 2) % n;
        //s[k] = k;
    }
#else
    int s[] = { 2, 4, 3, 1, 6, 7, 8, 9, 1, 7 };
#endif
    printV(s, n, "s");
    //test_compute_tree(s, n);

    for (int i = 0; i <= 3 * n / 2; i++) {
        for (int j = i; j < n; j++) {

            int rmq_ij = i;
            int rmq_v = s[i];
            for (int k = i + 1; k <= j; k++) {
                if (s[k] < rmq_v) {
                    rmq_ij = k;
                    rmq_v = s[k];
                }
            }

            if (j <= 10) {
                cout << "i=" << i << ",j=" << j << ",rmq=" << rmq_ij << ",v=" << rmq_v << endl;
            }
            int rmq_k = test_rmq_3(s, n, i, j, rmq_ij);
            if (j <= 10) {
                cout << "rmq=" << rmq_k << endl;
            }
            if (rmq_k != rmq_ij) {
                cerr << "BAD!\n";
                goto Done;
            }
        }
    }

    ok = true;
Done:
    delete[] s;
    return ok;
}
