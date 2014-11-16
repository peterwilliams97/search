// http://people.mpi-inf.mpg.de/~sanders/programs/suffix/drittel.C

#include <vector>
#include <algorithm>
using namespace std;
#include "util.h"

// lexic. order for pairs
inline bool leq(int a1, int a2, int b1, int b2) {
    return a1 < b1 || (a1 == b1 && a2 <= b2);
}

// and triples
inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3) {
    return a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3));
}

// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K-1 in r[]
static void radixPass(const int *a, int *b, const int *r, int n, int K) {
    // count occurrences
    int *c = new int[K + 1];                        // counter array
    for (int i = 0; i <= K; i++) c[i] = 0;          // reset counters
    for (int i = 0; i < n; i++) c[r[a[i]]]++;       // c[k] <- #occurrences of r[k]  in a[]
    for (int i = 0, sum = 0; i <= K; i++) {         // exclusive prefix sums
        int t = c[i];  c[i] = sum;  sum += t;       // c[k] <- #occurrences of r[0]..r[k-1]  in a[]
    }                                               // ie offset of 1st occurrence of k in sorted array
    for (int i = 0; i < n; i++) {
        b[c[r[a[i]]]++] = a[i];                     // sort b[i]
    }
    delete[] c;
}

// find the suffix array SA of s[0..n-1] in {1..K}^n
// require s[n]=s[n+1]=s[n+2]=0, n>=2
void suffixArray(const int *s, int n, int K, int *SA) {

    int n0 = (n + 2) / 3, n1 = (n + 1) / 3, n2 = n / 3, n02 = n0 + n2;
    int* s12 = new int[n02 + 3];  s12[n02] = s12[n02 + 1] = s12[n02 + 2] = 0;
    int* SA12 = new int[n02 + 3]; SA12[n02] = SA12[n02 + 1] = SA12[n02 + 2] = 0;
    int* s0 = new int[n0];
    int* SA0 = new int[n0];

    // generate positions of mod 1 and mod  2 suffixes
    // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
    for (int i = 0, j = 0; i < n + (n0 - n1); i++) {
        if (i % 3 != 0)
            s12[j++] = i;
    }

    // lsb radix sort the mod 1 and mod 2 triples
    radixPass(s12, SA12, s + 2, n02, K);
    radixPass(SA12, s12, s + 1, n02, K);
    radixPass(s12, SA12, s + 0, n02, K);

    // find lexicographic names of triples
    int name = 0, c0 = -1, c1 = -1, c2 = -1;
    for (int i = 0; i < n02; i++) {
        if (s[SA12[i]] != c0 || s[SA12[i] + 1] != c1 || s[SA12[i] + 2] != c2) {
            name++;
            c0 = s[SA12[i]];
            c1 = s[SA12[i] + 1];
            c2 = s[SA12[i] + 2];
        }
        if (SA12[i] % 3 == 1) {
            s12[SA12[i] / 3] = name;  // left half
         } else {
            s12[SA12[i] / 3 + n0] = name;  // right half
        }
    }

    // recurse if names are not yet unique
    if (name < n02) {
        suffixArray(s12, n02, name, SA12);
        // store unique names in s12 using the suffix array
        for (int i = 0; i < n02; i++) {
            s12[SA12[i]] = i + 1;
        }
    }
    else {  // generate the suffix array of s12 directly
        for (int i = 0; i < n02; i++) {
            SA12[s12[i] - 1] = i;
        }
    }

    // stably sort the mod 0 suffixes from SA12 by their first character
    for (int i = 0, j = 0; i < n02; i++) {
        if (SA12[i] < n0) {
            s0[j++] = 3 * SA12[i];
        }
    }
    radixPass(s0, SA0, s, n0, K);

    // merge sorted SA0 suffixes and sorted SA12 suffixes
    for (int p = 0, t = n0 - n1, k = 0; k < n; k++) {

#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)

        int i = GetI(); // pos of current offset 12 suffix
        int j = SA0[p]; // pos of current offset 0  suffix
        if (SA12[t] < n0 ?
            leq(s[i], s12[SA12[t] + n0], s[j], s12[j / 3]) :
            leq(s[i], s[i + 1], s12[SA12[t] - n0 + 1], s[j], s[j + 1], s12[j / 3 + n0]))
        { // suffix from SA12 is smaller
            SA[k] = i;  t++;
            if (t == n02) { // done --- only SA0 suffixes left
                for (k++; p < n0; p++, k++) SA[k] = SA0[p];
            }
        }
        else {
            SA[k] = j;  p++;
            if (p == n0)  { // done --- only SA12 suffixes left
                for (k++; t < n02; t++, k++) SA[k] = GetI();
            }
        }
    }

    delete[] s12;
    delete[] SA12;
    delete[] SA0;
    delete[] s0;
}

void calc_lcp(const int *s, const int *sa, int n, int *lcp) {
    int *rank = new int[n];

    for (int i = 0; i < n; i++) {
        rank[sa[i]] = i;
    }
    for (int i = 0, h = 0; i < n; i++) {
        if (rank[i] < n - 1) {
            for (int j = sa[rank[i] + 1]; s[i + h] == s[j + h]; h++) {
                ;;
            }
            lcp[rank[i]] = h;
            if (h > 0) {
                --h;
            }
        }
    }
    delete[] rank;
}

/*
 * http://webglimpse.net/pubs/suffix.pdf

  height(i) = lcp(A[Pos[i-1]], A[Pos[i])], 
               1 <= i <= N-1,
               A is input string
               Pos[i] = position of ith lexicallu sorted suffix in A

 for i [1, N 1] such that BH[i] and Hgt[i] > N do
 { a <- Prm[Pos[i 1] + H]
   b <- Prm[Pos[i] + H]
   Set( i, H + Min_Height(min(a, b) + 1, max(a, b))) { these routines are defined below }
 }

 */

static int 
calclcp_recurse(const int *hgt, int l, int r, int *llcp, int *rlcp) {
    if (r - l > 1) {
        int m = (r + l) / 2;
        llcp[m] = calclcp_recurse(hgt, l, m, llcp, rlcp);
        rlcp[m] = calclcp_recurse(hgt, m, r, llcp, rlcp);
        return Min(llcp[m], rlcp[m]);
    }
    else {
        return hgt[r];
    }
}

void calc_lcp_lr(const int *lcp, int n, int *llcp, int *rlcp) {
    calclcp_recurse(lcp, 0, n - 1, llcp, rlcp);
}

void suffixArrayLcp(const int *s, int n, int K, int *sa, int *lcp) {
    suffixArray(s, n, K, sa);
    calc_lcp(s, sa, n, lcp);
}

void suffixArrayLcpLR(const int *s, int n, int K, int *sa, int *lcp, int *llcp, int *rlcp) {
    suffixArray(s, n, K, sa);
    calc_lcp(s, sa, n, lcp);
    calc_lcp_lr(lcp, n, llcp, rlcp);
}

static inline int
comp(const int *s, const int *sa, int n, const int *pattern, int p, int i) {
    int pos = sa[i];
    //int m = Min(n - pos, p);
    for (int k = 0; k < p; k++) {
        if (pos + k >= n) {
            return +1;
        }
        int d = pattern[k] - s[pos + k];
        if (d < 0) return -1;
        if (d > 0) return +1;
    }
    return 0;
}

/*
 * Find offset of string `pattern` of length `p` in string `s` of length n
 *  using suffix array sa
 */
vector<int>
sa_search(const int *s, const int *sa, const int *lcp, const int *llcp, const int *rlcp,
             int n,
             const int *pattern, int p) {
    int l = 0;
    int r = n - 1;
    int m;
    vector<int> matches;

    if (comp(s, sa, n, pattern, p, 0) < 0) {
        return matches;
    }
    if (comp(s, sa, n, pattern, p, n - 1) > 0) {
        return matches;
    }
    
    while (l < r) {
        m = (l + r) / 2;
        int d = comp(s, sa, n, pattern, p, m);
        if (d < 0) {
            r = m;
        } else if (d > 0) {
            l = m;
        } else {
            matches.push_back(sa[m]);
            for (int k = m + 1; k <= r; k++) {
                if (comp(s, sa, n, pattern, p, k) != 0) break;
                matches.push_back(sa[k]);
            }
            for (int k = m - 1; k >= 1; k--) {
                if (comp(s, sa, n, pattern, p, k) != 0) break;
                matches.push_back(sa[k]);
            }
            sort(matches.begin(), matches.end());
            break;
        }
    }

    return matches;
}
