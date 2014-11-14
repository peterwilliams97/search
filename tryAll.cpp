#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
using namespace std;
#define DEBUGLEVEL 1
#include "util.h"

void suffixArray(const int *s, int *SA, int n, int K);
void suffixArrayLcp(const int *s, int *SA, int *lcp, int n, int K);

int min(int a, int b) {
    return a < b ? a : b;
}

void printV(const int *a, int n, const char *comment) {
    cout << comment << ": n= " << n << " : ";
    for (int i = 0; i < min(n, 20); i++) {
        cout << a[i] << " ";
    }
    cout << endl;
}

bool isPermutation(const int *SA, int n) {
    bool *seen = new bool[n];
    for (int i = 0; i < n; i++) seen[i] = 0;
    for (int i = 0; i < n; i++) seen[SA[i]] = 1;
    for (int i = 0; i < n; i++) if (!seen[i]) return 0;
    return 1;
}

bool sleq(const int *s1, const int *s2) {
#if 0
    if (s1[0] < s2[0]) return 1;
    if (s1[0] > s2[0]) return 0;
    return sleq(s1 + 1, s2 + 1);
#else
    for (int i = 0;; i++) {
        if (s1[i] < s2[i]) return 1;
        if (s1[i] > s2[i]) return 0;
    }
#endif
}

// is SA a sorted suffix array for s?
bool isSorted(const int *SA, const int *s, int n) {
    for (int i = 0; i < n - 1; i++) {
        if (!sleq(s + SA[i], s + SA[i + 1])) return 0;
    }
    return 1;
}

double 
test_one(const int *s, int *SA, int *lcp, int n, int b, bool do_check=false) {

    cout << "-------------------\n";
    Debug1(printV(s, n, "s"));

    clock_t start = clock();
    suffixArrayLcp(s, SA, lcp, n, b);
    double duration = (double)(clock() - start) / CLOCKS_PER_SEC;


    if (do_check) {
        Assert0(s[n] == 0);
        Assert0(s[n + 1] == 0);
        Assert0(SA[n] == 0);
        Assert0(SA[n + 1] == 0);
        Assert0(isPermutation(SA, n));
        Assert0(isSorted(SA, s, n));
        Debug1(printV(SA, n, "SA"));
    }

    return duration;
}

double 
test_n_b(int n, int b) {

    int *s = new int[n + 3];
    int *SA = new int[n + 3];
    int *lcp = new int[n + 3];
    for (int i = 0; i < n; i++) {
        lcp[i] = SA[i] = -1;
    }
    for (int i = 0; i < n; i++) {
        s[i] = i % b;
    }
    s[n] = s[n + 1] = s[n + 2] = SA[n] = SA[n + 1] = SA[n + 2] = 0;


    double duration =  test_one(s, SA, lcp, n, b);

    delete[] s;
    delete[] SA;
    delete[] lcp;

    return duration;
}

void 
test_n_all_b(int n, int b) {

    int N = int(pow(double(b), n) + 0.5);
    int *s = new int[n + 3];
    int *SA = new int[n + 3];
    int *lcp = new int[n + 3];
    for (int i = 0; i < n; i++) {
        lcp[i] = s[i] = SA[i] = -1;
    }
    s[n] = s[n + 1] = s[n + 2] = SA[n] = SA[n + 1] = SA[n + 2] = 0;
    for (int j = 0; j < N; j++) {
        test_one(s, SA, lcp, n, b);

        // generate next s
        int i;
        for (i = 0; s[i] == b; i++) {
            s[i] = 1;
        }
        s[i]++;
    }
    delete[] s;
    delete[] SA;
    delete[] lcp;
}

// try all inputs from {1,..,b}^n for 1 <= n <= nmax
int main(int argc, char **argv) {

    //int n = 13;
    //int s1[] = {2,1,4,4,1,4,4,1,3,3,1,0,0,0}; // mississippi
    //int s2[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //int n = 8;
    //int s1[] = {2,1,3,1,3,1,0,0,0}; // banana
    //int s2[] = {0,0,0,0,0,0,0,0,0};

    int n, b;
    b = 256;
    n = 1000; 

    for (n = 16 * 1024; n <= 1024 * 1024 * 1024; n *= 2) {
        cout << "=====================" << endl;
        cout << "n=" << n << ",b=" << b << endl;
        for (int i = 0; i < 1; i++) {
            double duration = test_n_b(n, b);
            int rate = (int)((double)n / duration);
            cout << "Time  = " << duration << endl;
            cout << "n/sec = " << rate << endl;
        }
    }

#if 0
    int nmax = atoi(argv[1]);
    int b = atoi(argv[2]);
    // try all strings from (1..b)^n
    for (int n = 2; n <= nmax; n++) {
        cout << "n=" << n << ",b=" << b << endl;
        test_n_all_b(n, b);
    }
#endif
    cout << "Done\n";
}