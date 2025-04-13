#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <functional>

using namespace std;

// -------------------- Data Structures --------------------

struct Point {
    double weight;
    vector<double> coords;

    bool operator==(const Point &other) const {
        return weight == other.weight && coords == other.coords;
    }
};

bool compareByWeight(const Point &a, const Point &b) {
    return a.weight > b.weight;
}

// Checks if point p lies inside query rectangle R
bool isInside(const Point &p, const vector<pair<double, double>> &R) {
    for (int i = 0; i < R.size(); i++) {
        if (p.coords[i] < R[i].first || p.coords[i] > R[i].second)
            return false;
    }
    return true;
}

// -------------------- BBD Tree --------------------

// BBDNode represents a node in the BBD-tree.
struct BBDNode {
    vector<Point> points;
    BBDNode *left = nullptr;
    BBDNode *right = nullptr;
    bool isLeaf = true;

    BBDNode(const vector<Point> &pts) : points(pts) {}
};

// BBDTree partitions the points hierarchically.
struct BBDTree {
    BBDNode *root;

    BBDTree(const vector<Point> &pts) {
        root = build(pts);
    }

    // Recursively build the BBD-tree.
    BBDNode* build(const vector<Point> &pts) {
        if (pts.size() <= 2) {
            return new BBDNode(pts);
        }
        // Sort by the first coordinate
        vector<Point> sorted = pts;
        sort(sorted.begin(), sorted.end(), [](const Point &a, const Point &b) {
            return a.coords[0] < b.coords[0];
        });
        int mid = sorted.size() / 2;
        vector<Point> leftPoints(sorted.begin(), sorted.begin() + mid);
        vector<Point> rightPoints(sorted.begin() + mid, sorted.end());
        BBDNode *node = new BBDNode({});
        node->isLeaf = false;
        node->left = build(leftPoints);
        node->right = build(rightPoints);
        return node;
    }

    // Collect all nodes in the tree (for debugging or filtering purposes)
    vector<BBDNode*> allNodes() {
        vector<BBDNode*> nodes;
        function<void(BBDNode*)> collect = [&](BBDNode* node) {
            if (!node) return;
            nodes.push_back(node);
            collect(node->left);
            collect(node->right);
        };
        collect(root);
        return nodes;
    }
};

// -------------------- Cover Region --------------------

// coverRegion returns all BBDNodes (leaf nodes) whose points are within radius of point q.
// Radius = (1 + epsilon)*delta.
vector<BBDNode*> coverRegion(BBDTree &tree, const Point &q, double delta, double epsilon) {
    vector<BBDNode*> covered;
    double radius = (1 + epsilon) * delta;
    cout << "[coverRegion] Using radius: " << radius << endl;

    function<void(BBDNode*)> dfs = [&](BBDNode* node) {
        if (!node) return;
        if (node->isLeaf) {
            for (const Point &p : node->points) {
                double dist = 0;
                for (int i = 0; i < p.coords.size(); ++i) {
                    dist += pow(p.coords[i] - q.coords[i], 2);
                }
                double euclidean = sqrt(dist);
                cout << "  Checking point (weight " << p.weight << ") at distance " << euclidean << endl;
                if (euclidean <= radius) {
                    cout << "   Adding leaf node (contains weight " << p.weight << ") to cover region" << endl;
                    covered.push_back(node);
                    return; // Add node once; assume its entire set is covered.
                }
            }
        } else {
            dfs(node->left);
            dfs(node->right);
        }
    };

    dfs(tree.root);
    if (covered.empty()) {
        cout << "  ⚠️ No nodes found within radius in coverRegion" << endl;
    }
    return covered;
}

// -------------------- Range Tree --------------------

// RangeTree stores all points sorted by weight.
struct RangeTree {
    vector<Point> points;

    void build(const vector<Point> &input) {
        points = input;
        sort(points.begin(), points.end(), compareByWeight);
        cout << "[RangeTree] Built with " << points.size() << " points (sorted by weight)" << endl;
    }

    // rangeMaxQuery: return the highest-weight point that is inside R and not contained in any node of Z.
    Point rangeMaxQuery(const vector<pair<double, double>> &R, const vector<BBDNode*> &Z) {
        for (const Point &p : points) {
            if (!isInside(p, R)) continue;
            bool skip = false;
            for (auto* node : Z) {
                // Here, we check if p is one of the points in the node.
                if (find(node->points.begin(), node->points.end(), p) != node->points.end()) {
                    skip = true;
                    break;
                }
            }
            if (skip) {
                cout << "   Skipping point (weight " << p.weight << ") as it is in a covered region (Z)" << endl;
                continue;
            }
            cout << "   Candidate from RangeTree: weight " << p.weight << " at (";
            for (double c : p.coords) cout << c << " ";
            cout << ")" << endl;
            return p;
        }
        return Point{-1, {}}; // No valid point found.
    }

    // Remove point p from the tree's point list.
    void erasePoint(const Point &p) {
        auto it = remove_if(points.begin(), points.end(), [&](const Point &q) {
            return q == p;
        });
        points.erase(it, points.end());
    }
};

// -------------------- CDT Query Procedure --------------------

// queryKHeaviestDiverse follows the paper's algorithm:
//   Q ⊆ P∩R is maintained; initially Q = ∅.
//   Z is the set of nodes covering previously selected points' δ-balls.
//   In each iteration, we select the highest-weight point not in ∪Z.
vector<Point> queryKHeaviestDiverse(BBDTree &tree, RangeTree &rangeTree, int k, double delta, double epsilon, const vector<pair<double, double>> &R) {
    vector<Point> result;
    vector<BBDNode*> Z; // Covered nodes set

    // INITIALIZATION: Get initial candidate from all points (Z is initially empty)
    Point initial = rangeTree.rangeMaxQuery(R, Z);
    if (initial.weight == -1) {
        cout << " No initial point found in range." << endl;
        return result;
    }
    cout << "Initial point: " << initial.weight << " at (";
    for (double c : initial.coords) cout << c << " ";
    cout << ")" << endl;

    result.push_back(initial);
    // Compute cover region for the initial point and add to Z.
    vector<BBDNode*> initCover = coverRegion(tree, initial, delta, epsilon);
    for (auto* node : initCover) {
        if (find(Z.begin(), Z.end(), node) == Z.end()) {
            Z.push_back(node);
        }
    }
    rangeTree.erasePoint(initial);

    cout << "After initialization, Z (cover region nodes) contains:" << endl;
    for (auto* z : Z) {
        for (const auto &p : z->points) {
            cout << "  " << p.weight << " at (";
            for (double c : p.coords) cout << c << " ";
            cout << ")" << endl;
        }
    }
    cout << "Result so far:" << endl;
    for (const auto &p : result) {
        cout << "  " << p.weight << " at (";
        for (double c : p.coords) cout << c << " ";
        cout << ")" << endl;
    }

    // ITERATIVE SELECTION
    int iter = 1;
    while (result.size() < k) {
        cout << "\n--- Iteration " << iter++ << " ---" << endl;
        Point candidate = rangeTree.rangeMaxQuery(R, Z);
        if (candidate.weight == -1) {
            cout << " No more valid candidates found outside covered regions." << endl;
            break;
        }
        cout << "Selected candidate: " << candidate.weight << " at (";
        for (double c : candidate.coords) cout << c << " ";
        cout << ")" << endl;

        result.push_back(candidate);
        rangeTree.erasePoint(candidate);

        // Update cover region: compute cover region for candidate and add to Z.
        vector<BBDNode*> newCover = coverRegion(tree, candidate, delta, epsilon);
        cout << "New cover region nodes for candidate " << candidate.weight << ":" << endl;
        for (auto* node : newCover) {
            // If node is not already in Z, add it.
            if (find(Z.begin(), Z.end(), node) == Z.end()) {
                Z.push_back(node);
                cout << "  Added node containing: ";
                for (auto &p : node->points) {
                    cout << p.weight << " ";
                }
                cout << endl;
            }
        }

        cout << "Updated Z (cover region nodes):" << endl;
        for (auto* z : Z) {
            for (const auto &p : z->points) {
                cout << "  " << p.weight << " at (";
                for (double c : p.coords) cout << c << " ";
                cout << ")" << endl;
            }
        }

        cout << "Result set so far:" << endl;
        for (const auto &p : result) {
            cout << "  " << p.weight << " at (";
            for (double c : p.coords) cout << c << " ";
            cout << ")" << endl;
        }
    }

    return result;
}

// -------------------- Main Function --------------------

int main() {
    int n, k, d;
    double delta, epsilon = 0.1;
    cin >> n >> k >> delta >> d;

    vector<Point> points(n);
    for (int i = 0; i < n; i++) {
        cin >> points[i].weight;
        points[i].coords.resize(d);
        for (int j = 0; j < d; j++) {
            cin >> points[i].coords[j];
        }
    }

    vector<pair<double, double>> R(d);
    for (int i = 0; i < d; i++) {
        cin >> R[i].first >> R[i].second;
    }

    cout << "\n Building BBD Tree..." << endl;
    BBDTree bbdTree(points);
    cout << " BBD Tree built successfully." << endl;

    cout << "\n Building Range Tree..." << endl;
    RangeTree rangeTree;
    rangeTree.build(points);
    cout << " Range Tree built successfully." << endl;

    cout << "\n Executing Query Procedure..." << endl;
    vector<Point> result = queryKHeaviestDiverse(bbdTree, rangeTree, k, delta, epsilon, R);

    cout << "\n Final Result:" << endl;
    for (const auto &p : result) {
        cout << p.weight << " ";
        for (double coord : p.coords) {
            cout << coord << " ";
        }
        cout << endl;
    }

    return 0;
}
