#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
using namespace std;
#define DEBUGLEVEL 1
#include "util.h"

void suffixArray(const int *s, int n, int K, int *SA);
void suffixArrayLcp(const int *s, int n, int K, int *sa, int *lcp);
void suffixArrayLcpLR(const int *s, int n, int K, int *sa, int *lcp, int *llcp, int *rlcp);
vector<int> sa_search(const int *s, const int *sa, const int *llcp, const int *rlcp,
              int n,
              const int *pattern, int p);
bool test_rmq();


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
    for (int i = 0; ; i++) {
        if (s1[i] < s2[i]) return 1;
        if (s1[i] > s2[i]) return 0;
    }
#endif
}

// is SA a sorted suffix array for s?
bool
isSorted(const int *SA, const int *s, int n) {
    for (int i = 0; i < n - 1; i++) {
        if (!sleq(s + SA[i], s + SA[i + 1])) {
            return false;
        }
    }
    return true;
}

static double
test_one(const int *s, int n, int b, bool do_check, int *SA, int *lcp, int *llcp, int *rlcp) {

    cout << "-------------------\n";
    Debug1(printV(s, n, "s"));

    clock_t start = clock();
    suffixArrayLcpLR(s, n, b, SA, lcp, llcp, rlcp);
    double duration = (double)(clock() - start) / CLOCKS_PER_SEC;

    int p = n / 3;
    int actual_offset = p / 2  +1;
    p = Min(p, 20);
    int *pattern = new int[p];
    for (int i = 0; i < p; i++) {
        pattern[i] = s[actual_offset + i];
    }
    Debug1(printV(pattern, p, "p"));

    vector<int> offsets = sa_search(s, SA, llcp, rlcp, n, pattern, p);
    cout << "offset actual = " << actual_offset << ", detected = " << offsets.size() << ": " ;
    for (int i = 0; i < Min(offsets.size(), 10); i++) {
        cout << offsets[i] << ", ";
    }
    cout << endl;
    Assert0(find(offsets.begin(), offsets.end(), actual_offset) != offsets.end());


    if (do_check) {
        Assert0(s[n] == 0);
        Assert0(s[n + 1] == 0);
        Assert0(SA[n] == 0);
        Assert0(SA[n + 1] == 0);
        Assert0(isPermutation(SA, n));
        Assert0(isSorted(SA, s, n));
        Debug1(printV(SA, n, "SA"));
        Debug1(printV(lcp, n, "lcp"));
    }

    return duration;
}


double
test_n_b(int n, int b) {

    bool do_check = false;

    int *s = new int[n + 3];
    int *SA = new int[n + 3];
    int *lcp = new int[n + 3];
    int *llcp = new int[n + 3];
    int *rlcp = new int[n + 3];
    for (int i = 0; i < n; i++) {
        llcp[i] = rlcp[i] = lcp[i] = SA[i] = -1;
    }
    for (int i = 0; i < n; i++) {
        s[i] = i % b;
    }
    s[n] = s[n + 1] = s[n + 2] = SA[n] = SA[n + 1] = SA[n + 2] = 0;

    double duration = test_one(s, n, b, false, SA, lcp, llcp, rlcp);

    delete[] s;
    delete[] SA;
    delete[] lcp;
    delete[] llcp;
    delete[] rlcp;

    return duration;
}

void
test_sa(const char *text) {
    int n = (int)strlen(text);
    int k = 0;
    map<char, int> char_int;
    vector<char> chars;

    for (const char *t = text; *t; t++) {
        if (find(chars.begin(), chars.end(), *t) != chars.end()) {
            continue;
        }
        chars.push_back(*t);
    }

    k = (int)chars.size();
    sort(chars.begin(), chars.end());

    for (int i = 0; i < k; i++) {
        char_int[chars[i]] = i;
    }

    int *s = new int[n + 3];
    int *sa = new int[n + 3];
    int *lcp = new int[n + 3];
    int *llcp = new int[n + 3];
    int *rlcp = new int[n + 3];

    for (int i = 0; i < n + 3; i++) {
        s[i] = sa[i] = lcp[i] = 0;
    }

    for (int i = 0; i < n; i++) {
        s[i] = char_int[text[i]];
    }

    test_one(s, n, k, true, sa, lcp, llcp, rlcp);
    cout << "Done: " << text << endl;
}

void
banana() {
    test_sa("banana");
}

void
test_n_all_b(int n, int b) {

    int N = int(pow(double(b), n) + 0.5);
    int *s = new int[n + 3];
    int *SA = new int[n + 3];
    int *lcp = new int[n + 3];
    int *llcp = new int[n + 3];
    int *rlcp = new int[n + 3];
    for (int i = 0; i < n; i++) {
        llcp[i] = rlcp[i] = lcp[i] = SA[i] = -1;
    }
    for (int i = n; i < n + 3 ; i++) {
        lcp[i] = s[i] = SA[i] = 0;
    }

    for (int j = 0; j < N; j++) {
        test_one(s, n, b, false, SA, lcp, llcp, rlcp);

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
    delete[] llcp;
    delete[] rlcp;
}

// try all inputs from {1,..,b}^n for 1 <= n <= nmax
int main(int argc, char **argv) {

    //int n = 13;
    //int s1[] = {2,1,4,4,1,4,4,1,3,3,1,0,0,0}; // mississippi
    //int s2[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    //int n = 8;
    //int s1[] = {2,1,3,1,3,1,0,0,0}; // banana
    //int s2[] = {0,0,0,0,0,0,0,0,0};

#if 0
    test_rmq();
#endif

#if 1
    banana();

    int n, b;
    b = 256;
    n = 1000;

    for (n = 2; n <= 256; n *= 2) {
        cout << "=====================" << endl;
        cout << "n=" << n << ",b=" << b << endl;
        for (int i = 0; i < 1; i++) {
            double duration = test_n_b(n, 2);
            int rate = (int)((double)n / duration);
            cout << "Time  = " << duration << endl;
            cout << "n/sec = " << rate << endl;
        }
    }

    for (n = 4; n <= 1024 * 1024 * 1024; n *= 2) {
        cout << "=====================" << endl;
        cout << "n=" << n << ",b=" << b << endl;
        for (int i = 0; i < 1; i++) {
            double duration = test_n_b(n, b);
            int rate = (int)((double)n / duration);
            cout << "Time  = " << duration << endl;
            cout << "n/sec = " << rate << endl;
        }
    }
#endif

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
