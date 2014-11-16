/*
 * http://web.stanford.edu/class/cs97si/suffix-array.pdf
 *
 */

#include <iostream>
#include "util.h"
using namespace std;
#if 0
/*
 * Build naive RMQ look-up-table rmq_lut for sequence s
 * s[n]
 * rmq_lut[n][n]  rmq_lut[i][j] = rmq_lut(i, j)
 *
 * O(n^2) time and space construction
 * O(1) operation
 */
int lcp(int x, int y) {
    int k, ret = 0;
    if (x == y)
        return N - x;
    for (k = stp - 1; k >= 0 && x < N && y < N; k--) {
        if (P[k][x] == P[k][y]) {
            x += 1 << k;
            y += 1 << k;
            ret += 1 << k;
        }
    }
    return ret;
}
#endif
