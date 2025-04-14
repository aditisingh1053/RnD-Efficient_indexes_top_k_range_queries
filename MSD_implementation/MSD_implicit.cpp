#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <cassert>

using namespace std;

// ------------------------------
// Data structures for implicit representation

struct Point {
    int id;                   // unique identifier for the point
    vector<double> coords;    // d-dimensional coordinate
    double weight;
};

struct Box {
    int id;
    vector<Point> points;     // Points in this box sorted in decreasing order of weight
    int currentIdx;           // Pointer to the next candidate in the box (initially 0)
};

// ------------------------------
// Global parameters (used in score computation)
double g_lambda;
int g_k;

// ------------------------------
// Dot product function.
double dotProduct(const Point &p, const vector<double>& u) {
    double sum = 0.0;
    for (size_t i = 0; i < p.coords.size(); i++)
        sum += p.coords[i] * u[i];
    return sum;
}

// Score function for a point in a given direction u.
double scoreCandidate(const Point &p, const vector<double>& u) {
    return dotProduct(p, u) + (g_lambda / (g_k - 1)) * p.weight;
}

// ------------------------------
// Priority queue entry for the implicit version; now each entry comes
// from a particular box, and we record the candidate's index in that box.
struct PQEntry {
    double scoreVal;  // score computed using a specific direction u
    int boxId;        // from which box this candidate comes
    int candIdx;      // the index in box.points at the time of insertion
};

struct PQComp {
    bool operator()(const PQEntry &a, const PQEntry &b) {
        return a.scoreVal < b.scoreVal;  // max-heap: highest score first
    }
};

// Clean the PQ by checking if the candidate from a box is still current.
void cleanPQ(priority_queue<PQEntry, vector<PQEntry>, PQComp>& pq,
             const vector<Box>& boxes) {
    while (!pq.empty()) {
        const PQEntry &top = pq.top();
        // If the candidate is stale (i.e. the box has advanced its pointer) then pop.
        if (top.candIdx != boxes[top.boxId].currentIdx)
            pq.pop();
        else
            break;
    }
}

// Helper function: for a given direction u and a given box, if the box still has candidates,
// compute the candidate's score and push a new PQEntry into the given priority queue.
void pushCandidateForBox(const vector<double>& u, 
                         const Box &box, 
                         priority_queue<PQEntry, vector<PQEntry>, PQComp>& pq) {
    if (box.currentIdx < (int)box.points.size()) {
        const Point &p = box.points[box.currentIdx];
        double s = scoreCandidate(p, u);
        PQEntry entry { s, box.id, box.currentIdx };
        pq.push(entry);
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int s, d, k;
    double lambda, epsilon = 0.1;
    // Read number of boxes, dimensionality, k and lambda.
    // (In practice, the implicit input would be provided via the index Φ.)
    cin >> s >> d >> k >> lambda;
    g_lambda = lambda;
    g_k = k;

    // Read boxes. For each box, first read the number of points in that box.
    vector<Box> boxes(s);
    for (int i = 0; i < s; i++) {
        boxes[i].id = i;
        boxes[i].currentIdx = 0;
        int numPoints;
        cin >> numPoints;
        boxes[i].points.resize(numPoints);
        for (int j = 0; j < numPoints; j++) {
            boxes[i].points[j].id = i * 1000 + j; // Example: unique id constructed from box and local index.
            boxes[i].points[j].coords.resize(d);
            for (int l = 0; l < d; l++) {
                cin >> boxes[i].points[j].coords[l];
            }
            cin >> boxes[i].points[j].weight;
        }
        // It is assumed that the points in each box are given in decreasing order of weight.
    }

    // ------------------------------
    // Construct an epsilon-net on the sphere.
    // Here we use the d unit vectors and their opposites.
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

    // ------------------------------
    // Instead of having one PQ per explicit point, we now have one PQ per direction.
    // Each PQ initially holds one candidate per box (its current top point).
    vector<priority_queue<PQEntry, vector<PQEntry>, PQComp>> pqList(numDirections);
    for (int i = 0; i < numDirections; i++) {
        for (int b = 0; b < s; b++) {
            pushCandidateForBox(net[i], boxes[b], pqList[i]);
        }
    }

    // We no longer have a global "removed" vector.
    // Instead, once a candidate from a box is selected, that box’s pointer is advanced.
    vector<Point> solution;

    // Each round selects a symmetric pair from the net directions.
    // We run for k/2 rounds (if k is even; otherwise we add one more candidate at the end).
    int rounds = k / 2;
    for (int round = 0; round < rounds; round++) {
        double bestSum = -numeric_limits<double>::infinity();
        int bestDir = -1;
        pair<pair<int, int>, pair<int, int>> bestPairCandidate; 
        // bestPairCandidate.first = (boxId, candidate index) for candidate1
        // bestPairCandidate.second = (boxId, candidate index) for candidate2

        // For each symmetric pair of directions. (net was constructed as [u, -u, u, -u, ...])
        for (int i = 0; i < numDirections; i += 2) {
            int idx_u = i;
            int idx_neg = i + 1;

            cleanPQ(pqList[idx_u], boxes);
            cleanPQ(pqList[idx_neg], boxes);
            if (pqList[idx_u].empty() || pqList[idx_neg].empty()) continue;

            PQEntry top_u = pqList[idx_u].top();
            PQEntry top_neg = pqList[idx_neg].top();

            int candBox1 = top_u.boxId, candIdx1 = top_u.candIdx;
            int candBox2 = top_neg.boxId, candIdx2 = top_neg.candIdx;

            double score1 = scoreCandidate(boxes[candBox1].points[candIdx1], net[idx_u]);
            double score2 = scoreCandidate(boxes[candBox2].points[candIdx2], net[idx_neg]);
            double sumScore = score1 + score2;

            // If the two candidates come from the same box (and same candidate) then we try alternatives.
            if (candBox1 == candBox2 && candIdx1 == candIdx2) {
                // Temporarily remove top from each PQ to peek at the second best candidate.
                PQEntry second_u, second_neg;
                bool foundSecondU = false, foundSecondNeg = false;
                
                {
                    PQEntry cur = top_u;
                    pqList[idx_u].pop();
                    cleanPQ(pqList[idx_u], boxes);
                    if (!pqList[idx_u].empty()) {
                        second_u = pqList[idx_u].top();
                        foundSecondU = true;
                    }
                    pqList[idx_u].push(cur);
                }
                {
                    PQEntry cur = top_neg;
                    pqList[idx_neg].pop();
                    cleanPQ(pqList[idx_neg], boxes);
                    if (!pqList[idx_neg].empty()) {
                        second_neg = pqList[idx_neg].top();
                        foundSecondNeg = true;
                    }
                    pqList[idx_neg].push(cur);
                }
                double alt1 = -numeric_limits<double>::infinity(), alt2 = -numeric_limits<double>::infinity();
                int bestBoxA = -1, bestIdxA = -1, bestBoxB = -1, bestIdxB = -1;
                if (foundSecondNeg) {
                    alt1 = scoreCandidate(boxes[top_u.boxId].points[top_u.candIdx], net[idx_u]) +
                           scoreCandidate(boxes[second_neg.boxId].points[second_neg.candIdx], net[idx_neg]);
                    bestBoxA = top_u.boxId; bestIdxA = top_u.candIdx;
                    bestBoxB = second_neg.boxId; bestIdxB = second_neg.candIdx;
                }
                if (foundSecondU) {
                    alt2 = scoreCandidate(boxes[second_u.boxId].points[second_u.candIdx], net[idx_u]) +
                           scoreCandidate(boxes[top_neg.boxId].points[top_neg.candIdx], net[idx_neg]);
                    if (alt2 > alt1) {
                        bestBoxA = second_u.boxId; bestIdxA = second_u.candIdx;
                        bestBoxB = top_neg.boxId; bestIdxB = top_neg.candIdx;
                        alt1 = alt2;
                    }
                }
                sumScore = alt1;
                candBox1 = bestBoxA; candIdx1 = bestIdxA;
                candBox2 = bestBoxB; candIdx2 = bestIdxB;
            }
            // Update overall best symmetric pair if necessary.
            if (sumScore > bestSum) {
                bestSum = sumScore;
                bestDir = idx_u;
                bestPairCandidate = { {candBox1, candIdx1}, {candBox2, candIdx2} };
            }
        }
        if (bestDir == -1) break; // no candidate found

        // Retrieve the selected candidates.
        Point p1 = boxes[bestPairCandidate.first.first].points[bestPairCandidate.first.second];
        Point p2 = boxes[bestPairCandidate.second.first].points[bestPairCandidate.second.second];
        solution.push_back(p1);
        solution.push_back(p2);

        // When a candidate is chosen, advance the corresponding box iterator.
        {
            int b = bestPairCandidate.first.first;
            Box &box = boxes[b];
            if (box.currentIdx == bestPairCandidate.first.second)
                box.currentIdx++;
        }
        {
            int b = bestPairCandidate.second.first;
            Box &box = boxes[b];
            if (box.currentIdx == bestPairCandidate.second.second)
                box.currentIdx++;
        }
        
        // For every direction, push a new candidate from the boxes that just advanced (if available)
        for (int i = 0; i < numDirections; i++) {
            // For first candidate's box:
            pushCandidateForBox(net[i], boxes[bestPairCandidate.first.first], pqList[i]);
            // For second candidate's box:
            pushCandidateForBox(net[i], boxes[bestPairCandidate.second.first], pqList[i]);
            // Clean outdated candidates from these boxes will be done by cleanPQ later.
        }
    } // end rounds

    // If k is odd, select one more candidate: choose the best among all boxes using an arbitrary direction (say, net[0])
    if ((int)solution.size() < k) {
        double bestScore = -numeric_limits<double>::infinity();
        int bestBox = -1;
        int bestCand = -1;
        for (int b = 0; b < s; b++) {
            if (boxes[b].currentIdx < (int)boxes[b].points.size()) {
                double sVal = scoreCandidate(boxes[b].points[boxes[b].currentIdx], net[0]);
                if (sVal > bestScore) {
                    bestScore = sVal;
                    bestBox = b;
                    bestCand = boxes[b].currentIdx;
                }
            }
        }
        if (bestBox != -1) {
            solution.push_back(boxes[bestBox].points[bestCand]);
            boxes[bestBox].currentIdx++;
        }
    }

    // Output the solution: print the coordinates and weight for every selected point.
    for (const Point &pt : solution) {
        for (double coord : pt.coords)
            cout << coord << " ";
        cout << pt.weight << "\n";
    }
    return 0;
}
