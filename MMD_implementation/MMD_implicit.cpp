#include <bits/stdc++.h>
using namespace std;

struct Point {
    vector<double> coord;
    double weight;
    int id;
};

double euclideanDistance(const Point &a, const Point &b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.coord.size(); i++)
        sum += (a.coord[i] - b.coord[i]) * (a.coord[i] - b.coord[i]);
    return sqrt(sum);
}

double closestPairDistance(const vector<Point> &pts) {
    double best = 1e18;
    for (size_t i = 0; i < pts.size(); ++i)
        for (size_t j = i + 1; j < pts.size(); ++j)
            best = min(best, euclideanDistance(pts[i], pts[j]));
    return (pts.size() < 2 ? 0.0 : best);
}

double minWeight(const vector<Point> &pts) {
    double m = 1e18;
    for (auto &p : pts)
        m = min(m, p.weight);
    return m;
}

// --- Simulated Box & Phi interface ---
struct Box {
    int id;
};

struct Phi {
    unordered_map<int, vector<Point>> boxData;
    unordered_map<int, int> nextIndex;

    void addBoxData(int boxId, vector<Point> points) {
        sort(points.begin(), points.end(), [](auto &a, auto &b) {
            return a.weight > b.weight;
        });
        boxData[boxId] = points;
        nextIndex[boxId] = 0;
    }

    Point next(int boxId) {
        int &idx = nextIndex[boxId];
        if (idx >= boxData[boxId].size())
            return Point{{}, -1e9, -1};
        return boxData[boxId][idx++];
    }

    void reset(int boxId) {
        nextIndex[boxId] = 0;
    }

    vector<Point> getAllInBoundingBox() {
        vector<Point> all;
        for (auto &[id, pts] : boxData)
            all.insert(all.end(), pts.begin(), pts.end());
        return all;
    }
};

// --- Implicit MMD Implementation ---
int main() {

    int d, k, s; // dimension, subset size, number of boxes
    double lambda;
    cin >> s >> d >> k >> lambda;

    vector<Box> boxes(s);
    Phi phi;

    // Input: each box gets some points
    for (int i = 0; i < s; i++) {
        int m;
        cin >> m;
        vector<Point> pts(m);
        for (int j = 0; j < m; j++) {
            pts[j].coord.resize(d);
            for (int t = 0; t < d; t++) cin >> pts[j].coord[t];
            cin >> pts[j].weight;
            pts[j].id = i * INT16_MAX + j;  // unique id
        }
        phi.addBoxData(i, pts);
        boxes[i].id = i;
    }

    // --- Build H1: one heaviest point per box ---
    vector<Point> H1;
    for (auto &box : boxes) {
        Point p = phi.next(box.id);
        if (p.id != -1) H1.push_back(p);
    }

    // --- Build H2: k heaviest from bounding box ---
    vector<Point> allPoints = phi.getAllInBoundingBox();
    sort(allPoints.begin(), allPoints.end(), [](auto &a, auto &b) {
        return a.weight > b.weight;
    });
    vector<Point> H2(allPoints.begin(), allPoints.begin() + min(k, (int)allPoints.size()));

    // --- Combine and remove duplicates based on unique IDs ---
    vector<Point> H = H1;
    H.insert(H.end(), H2.begin(), H2.end());

    unordered_set<int> seen_ids;
    vector<Point> uniqueH;
    for (auto &pt : H) {
        if (seen_ids.insert(pt.id).second)
            uniqueH.push_back(pt);
    }
    H = uniqueH;

    // --- Brute-force all subsets of size k from H ---
    int nH = H.size();
    vector<Point> bestSolution;
    double bestUtility = -1e18;

    for (int mask = 0; mask < (1 << nH); ++mask) {
        if (__builtin_popcount(mask) != k) continue;
        vector<Point> subset;
        for (int j = 0; j < nH; ++j)
            if (mask & (1 << j)) subset.push_back(H[j]);

        double util = closestPairDistance(subset) + lambda * minWeight(subset);
        if (util > bestUtility) {
            bestUtility = util;
            bestSolution = subset;
        }
    }

    // --- Output ---
    for (auto &pt : bestSolution) {
        for (auto x : pt.coord) cout << x << " ";
        cout << pt.weight << "\n";
    }

    return 0;
}
