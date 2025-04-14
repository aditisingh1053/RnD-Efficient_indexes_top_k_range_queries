#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <cassert>

using namespace std;

struct Point {
    int id;                  
    vector<double> coords;  
    double weight;
};

double dotProduct(const Point &p, const vector<double>& u) {
    double sum = 0.0;
    for (size_t i = 0; i < p.coords.size(); i++)
        sum += p.coords[i] * u[i];
    return sum;
}

double score(const Point &p, const vector<double>& u, double lambda, int k) {
    return dotProduct(p, u) + (lambda / (k - 1)) * p.weight;
}

// Priority queue entry for a given direction u
struct PQEntry {
    double scoreVal;
    int pid;  // index of the point in global vector
};

struct PQComp {
    bool operator()(const PQEntry &a, const PQEntry &b) {
        return a.scoreVal < b.scoreVal;
    }
};

void cleanPQ(priority_queue<PQEntry, vector<PQEntry>, PQComp>& pq, const vector<bool>& removed) {
    while (!pq.empty() && removed[pq.top().pid])
        pq.pop();
}

int main() {
    int n, d, k;
    double lambda, epsilon=0.1;
    cin >> n >> d >> k >> lambda;

    vector<Point> points(n);
    for (int i = 0; i < n; i++) {
        points[i].id = i;
        points[i].coords.resize(d);
        for (int j = 0; j < d; j++) {
            cin >> points[i].coords[j];
        }
        cin >> points[i].weight;
    }

    // Generate an epsilon-net N on the sphere.
    vector<vector<double>> net;
    for (int i = 0; i < d; i++) {
        vector<double> u(d, 0.0);
        u[i] = 1.0;
        net.push_back(u);
        vector<double> neg_u = u;
        neg_u[i] = -1.0;
        net.push_back(neg_u);
    }

    int numDirections = net.size();
    vector<priority_queue<PQEntry, vector<PQEntry>, PQComp>> pqList(numDirections);
    for (int i = 0; i < numDirections; i++) {
        for (int j = 0; j < n; j++) {
            double s = score(points[j], net[i], lambda, k);
            PQEntry entry { s, j };
            pqList[i].push(entry);
        }
    }

    vector<bool> removed(n, false);
    vector<int> solution; 

    int rounds = k / 2;
    for (int round = 0; round < rounds; round++) {
        double bestSum = -numeric_limits<double>::infinity();
        int bestDir = -1;

        // For each symmetric pair: we consider indices in net in pairs.
        // Since we constructed net as [u1, -u1, u2, -u2, ...] we iterate with step 2.
        pair<int, int> bestPairCandidate(-1, -1); 
        for (int i = 0; i < numDirections; i += 2) {
            int idx_u = i;
            int idx_neg = i + 1;
            cleanPQ(pqList[idx_u], removed);
            cleanPQ(pqList[idx_neg], removed);
            if (pqList[idx_u].empty() || pqList[idx_neg].empty()) continue;

            PQEntry top_u = pqList[idx_u].top();
            PQEntry top_neg = pqList[idx_neg].top();

            int cand1 = top_u.pid, cand2 = top_neg.pid;
            double sumScore = score(points[cand1], net[idx_u], lambda, k) +
                              score(points[cand2], net[idx_neg], lambda, k);

            if (cand1 == cand2) {
                // Temporarily remove the top to check the next best candidate.
                PQEntry second_u, second_neg;
                bool foundSecondU = false, foundSecondNeg = false;

                // For queue u try to get the next best candidate
                {
                    PQEntry cur = top_u;
                    pqList[idx_u].pop();
                    cleanPQ(pqList[idx_u], removed);
                    if (!pqList[idx_u].empty()) {
                        second_u = pqList[idx_u].top();
                        foundSecondU = true;
                    }
                    // Push back the top element to restore the queue
                    pqList[idx_u].push(cur);
                }
                // Similarly for queue -u
                {
                    PQEntry cur = top_neg;
                    pqList[idx_neg].pop();
                    cleanPQ(pqList[idx_neg], removed);
                    if (!pqList[idx_neg].empty()) {
                        second_neg = pqList[idx_neg].top();
                        foundSecondNeg = true;
                    }
                    pqList[idx_neg].push(cur);
                }
                // Now try both alternatives if available.
                double alt1 = -numeric_limits<double>::infinity(), alt2 = -numeric_limits<double>::infinity();
                int cand1_alt1 = -1, cand2_alt1 = -1;
                int cand1_alt2 = -1, cand2_alt2 = -1;
                if (foundSecondNeg) {
                    // Use top_u and second_neg
                    alt1 = score(points[top_u.pid], net[idx_u], lambda, k) +
                           score(points[second_neg.pid], net[idx_neg], lambda, k);
                    cand1_alt1 = top_u.pid;
                    cand2_alt1 = second_neg.pid;
                }
                if (foundSecondU) {
                    // Use second_u and top_neg
                    alt2 = score(points[second_u.pid], net[idx_u], lambda, k) +
                           score(points[top_neg.pid], net[idx_neg], lambda, k);
                    cand1_alt2 = second_u.pid;
                    cand2_alt2 = top_neg.pid;
                }
                if (alt1 >= alt2) {
                    sumScore = alt1;
                    cand1 = cand1_alt1;
                    cand2 = cand2_alt1;
                } else {
                    sumScore = alt2;
                    cand1 = cand1_alt2;
                    cand2 = cand2_alt2;
                }
            }
            // Update the best symmetric pair if necessary.
            if (sumScore > bestSum) {
                bestSum = sumScore;
                bestDir = idx_u; 
                bestPairCandidate = {cand1, cand2};
            }
        }

        if (bestDir == -1) break;

        int p1 = bestPairCandidate.first, p2 = bestPairCandidate.second;
        solution.push_back(p1);
        solution.push_back(p2);
        removed[p1] = true;
        removed[p2] = true;

        for (int i = 0; i < numDirections; i++) {
            cleanPQ(pqList[i], removed);
        }
    }

    // If k is odd we can add one more highest scoring point among all remaining points.
    if (solution.size() < (size_t)k) {
        double best = -numeric_limits<double>::infinity();
        int bestID = -1;
        // For simplicity, consider all points
        for (int i = 0; i < n; i++) {
            if (!removed[i]) {
                double s = score(points[i], net[0], lambda, k);
                if (s > best) {
                    best = s;
                    bestID = i;
                }
            }
        }
        if (bestID != -1) {
            solution.push_back(bestID);
            removed[bestID] = true;
        }
    }

    for (int id : solution) {
        for (double coord : points[id].coords)
            cout << coord << " ";
        cout << points[id].weight << endl;
    }

    return 0;
}
