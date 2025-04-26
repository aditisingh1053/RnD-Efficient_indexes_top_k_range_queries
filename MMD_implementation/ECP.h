// ---------------------------------------------------------------------------
// Dynamic ECP Data Structure
//
// This class maintains a set of points with a cached closest pair distance.
// Upon every insertion or deletion the “dirty” flag is marked. On query, if dirty,
// the closest pair is recomputed in O(n^2) time (for simplicity), so that subsequent
// queries can be answered in O(1) time.
class DynamicECP {
private:
    vector<Point> points;
    double cachedClosest;
    bool dirty;
    
    double computeClosest() {
        return bruteForceClosestPair(points);
    }
public:
    DynamicECP() : cachedClosest(numeric_limits<double>::max()), dirty(true) {}
    
    void insert(const Point &p) {
        points.push_back(p);
        dirty = true;
    }
    
    void remove(const Point &p) {
        for (auto it = points.begin(); it != points.end(); ++it) {
            if(it->id == p.id) {
                points.erase(it);
                dirty = true;
                break;
            }
        }
    }
    
    // Returns the cached closest pair distance (O(1) if not dirty).
    double getClosestPair() {
        if (dirty) {
            cachedClosest = computeClosest();
            dirty = false;
        }
        return cachedClosest;
    }
    
    // Get the current set of points (for possible iteration in pruning).
    const vector<Point>& getPoints() const { return points; }
};