#include "clipper2/clipper.core.h"
#include "clipper2/clipper.h"
#include "clipper2/clipper.offset.h"
#include <emscripten/bind.h>
#include <cmath>
#include <algorithm>

using namespace emscripten;
using namespace Clipper2Lib;

template <typename T>
void ReversePath(Path<T>& path) {
    std::reverse(path.begin(), path.end());
}


// Helper to convert PathD to Path64 with given precision
static Path64 PathDToPath64(const PathD& pathD, double scale) {
    Path64 result;
    result.reserve(pathD.size());
    for (const auto& pt : pathD) {
#ifdef USINGZ
        result.push_back(Point64(
            static_cast<int64_t>(pt.x * scale),
            static_cast<int64_t>(pt.y * scale),
            pt.z));
#else
        result.push_back(Point64(
            static_cast<int64_t>(pt.x * scale),
            static_cast<int64_t>(pt.y * scale)));
#endif
    }
    return result;
}

// Helper to convert Path64 to PathD with given inverse scale
static PathD Path64ToPathD(const Path64& path64, double invScale) {
    PathD result;
    result.reserve(path64.size());
    for (const auto& pt : path64) {
#ifdef USINGZ
        result.push_back(PointD(
            pt.x * invScale,
            pt.y * invScale,
            pt.z));
#else
        result.push_back(PointD(
            pt.x * invScale,
            pt.y * invScale));
#endif
    }
    return result;
}

// Helper to convert PathsD to Paths64
static Paths64 PathsDToPaths64(const PathsD& pathsD, double scale) {
    Paths64 result;
    result.reserve(pathsD.size());
    for (const auto& path : pathsD) {
        result.push_back(PathDToPath64(path, scale));
    }
    return result;
}

// Helper to convert Paths64 to PathsD
static PathsD Paths64ToPathsD(const Paths64& paths64, double invScale) {
    PathsD result;
    result.reserve(paths64.size());
    for (const auto& path : paths64) {
        result.push_back(Path64ToPathD(path, invScale));
    }
    return result;
}

// Auto-scaling constants
// Use 10^15 as safe max to leave headroom for Clipper2's internal calculations
// (Clipper2's actual limit is INT64_MAX >> 2 ≈ 4.6 * 10^18)
static const double SAFE_MAX_COORD = 1.0e15;
static const double MAX_SCALE = 1.0e15;  // Cap scale to avoid issues with tiny coords
static const double MIN_SCALE = 1.0;     // Floor to avoid issues with huge coords

// Find the maximum absolute coordinate value in a PathD
static double GetMaxAbsCoord(const PathD& path) {
    double maxVal = 0.0;
    for (const auto& pt : path) {
        maxVal = std::max(maxVal, std::abs(pt.x));
        maxVal = std::max(maxVal, std::abs(pt.y));
    }
    return maxVal;
}

// Find the maximum absolute coordinate value across all PathsD
static double GetMaxAbsCoord(const PathsD& paths) {
    double maxVal = 0.0;
    for (const auto& path : paths) {
        maxVal = std::max(maxVal, GetMaxAbsCoord(path));
    }
    return maxVal;
}

// Compute optimal scale factor based on coordinate range and delta
// This maximizes precision while staying within int64 safe bounds
static double ComputeOptimalScale(double maxAbsCoord, double absDelta) {
    // Account for the offset expanding the path
    double maxExpectedCoord = maxAbsCoord + absDelta;

    // Handle edge case of zero/tiny coordinates
    if (maxExpectedCoord < 1.0e-10) {
        return MAX_SCALE;
    }

    // Compute scale: we want maxExpectedCoord * scale <= SAFE_MAX_COORD
    double scale = SAFE_MAX_COORD / maxExpectedCoord;

    // Clamp to reasonable bounds
    scale = std::min(scale, MAX_SCALE);
    scale = std::max(scale, MIN_SCALE);

    return scale;
}

// Single-call offset function with automatic precision scaling
PathsD InflatePathsDAutoScale(const PathsD& paths, double delta,
    JoinType jt, EndType et, double miterLimit = 2.0, double arcTolerance = 0.0)
{
    if (paths.empty()) return PathsD();

    double maxAbsCoord = GetMaxAbsCoord(paths);
    double scale = ComputeOptimalScale(maxAbsCoord, std::abs(delta));
    double invScale = 1.0 / scale;

    ClipperOffset co(miterLimit, arcTolerance * scale);
    co.AddPaths(PathsDToPaths64(paths, scale), jt, et);

    Paths64 result64;
    co.Execute(delta * scale, result64);

    return Paths64ToPathsD(result64, invScale);
}

// Compute optimal scale for boolean operations (no delta involved)
static double ComputeOptimalScaleForBooleanOp(double maxAbsCoord) {
    // Handle edge case of zero/tiny coordinates
    if (maxAbsCoord < 1.0e-10) {
        return MAX_SCALE;
    }

    // Compute scale: we want maxAbsCoord * scale <= SAFE_MAX_COORD
    double scale = SAFE_MAX_COORD / maxAbsCoord;

    // Clamp to reasonable bounds
    scale = std::min(scale, MAX_SCALE);
    scale = std::max(scale, MIN_SCALE);

    return scale;
}

// Auto-scaling boolean operation: Union
PathsD UnionDAutoScale(const PathsD& subjects, const PathsD& clips, FillRule fillRule)
{
    double maxAbsCoord = std::max(GetMaxAbsCoord(subjects), GetMaxAbsCoord(clips));
    if (maxAbsCoord == 0.0 && subjects.empty() && clips.empty()) return PathsD();

    double scale = ComputeOptimalScaleForBooleanOp(maxAbsCoord);
    double invScale = 1.0 / scale;

    Paths64 subjects64 = PathsDToPaths64(subjects, scale);
    Paths64 clips64 = PathsDToPaths64(clips, scale);

    Paths64 result64 = Union(subjects64, clips64, fillRule);

    return Paths64ToPathsD(result64, invScale);
}

// Auto-scaling boolean operation: Union (self - single path set)
PathsD UnionSelfDAutoScale(const PathsD& subjects, FillRule fillRule)
{
    if (subjects.empty()) return PathsD();

    double maxAbsCoord = GetMaxAbsCoord(subjects);
    double scale = ComputeOptimalScaleForBooleanOp(maxAbsCoord);
    double invScale = 1.0 / scale;

    Paths64 subjects64 = PathsDToPaths64(subjects, scale);

    Paths64 result64 = Union(subjects64, fillRule);

    return Paths64ToPathsD(result64, invScale);
}

// Auto-scaling boolean operation: Intersect
PathsD IntersectDAutoScale(const PathsD& subjects, const PathsD& clips, FillRule fillRule)
{
    if (subjects.empty() || clips.empty()) return PathsD();

    double maxAbsCoord = std::max(GetMaxAbsCoord(subjects), GetMaxAbsCoord(clips));
    double scale = ComputeOptimalScaleForBooleanOp(maxAbsCoord);
    double invScale = 1.0 / scale;

    Paths64 subjects64 = PathsDToPaths64(subjects, scale);
    Paths64 clips64 = PathsDToPaths64(clips, scale);

    Paths64 result64 = Intersect(subjects64, clips64, fillRule);

    return Paths64ToPathsD(result64, invScale);
}

// Auto-scaling boolean operation: Difference
PathsD DifferenceDAutoScale(const PathsD& subjects, const PathsD& clips, FillRule fillRule)
{
    if (subjects.empty()) return PathsD();

    double maxAbsCoord = std::max(GetMaxAbsCoord(subjects), GetMaxAbsCoord(clips));
    double scale = ComputeOptimalScaleForBooleanOp(maxAbsCoord);
    double invScale = 1.0 / scale;

    Paths64 subjects64 = PathsDToPaths64(subjects, scale);
    Paths64 clips64 = PathsDToPaths64(clips, scale);

    Paths64 result64 = Difference(subjects64, clips64, fillRule);

    return Paths64ToPathsD(result64, invScale);
}

// Auto-scaling boolean operation: Xor
PathsD XorDAutoScale(const PathsD& subjects, const PathsD& clips, FillRule fillRule)
{
    double maxAbsCoord = std::max(GetMaxAbsCoord(subjects), GetMaxAbsCoord(clips));
    if (maxAbsCoord == 0.0 && subjects.empty() && clips.empty()) return PathsD();

    double scale = ComputeOptimalScaleForBooleanOp(maxAbsCoord);
    double invScale = 1.0 / scale;

    Paths64 subjects64 = PathsDToPaths64(subjects, scale);
    Paths64 clips64 = PathsDToPaths64(clips, scale);

    Paths64 result64 = Xor(subjects64, clips64, fillRule);

    return Paths64ToPathsD(result64, invScale);
}

// Normalize a path's start point by rotating it so the vertex with the
// most negative Y comes first (lowest point). Ties are broken by most
// negative X (leftmost). This matches the convention Clipper2 uses
// internally for orientation detection, and the extreme point is always
// on the outer surface of the path by definition.
static void NormalizePathStart(PathD& path) {
    if (path.size() < 3) {
        return;
    }

    size_t bestIdx = 0;
    for (size_t i = 1; i < path.size(); ++i) {
        if (path[i].y < path[bestIdx].y ||
            (path[i].y == path[bestIdx].y && path[i].x < path[bestIdx].x)) {
            bestIdx = i;
        }
    }

    if (bestIdx != 0) {
        std::rotate(path.begin(), path.begin() + bestIdx, path.end());
    }
}

static void NormalizePathsStart(PathsD& paths) {
    for (auto& path : paths) {
        NormalizePathStart(path);
    }
}

// Stored path entry for deferred scaling
struct PathEntry {
    PathD path;
    JoinType joinType;
    EndType endType;
};

// Wrapper class for ClipperOffset that handles auto-scaling transparently
// Paths are stored and scaled at Execute time when the delta is known
class ClipperOffsetD {
private:
    std::vector<PathEntry> storedPaths_;
    double maxAbsCoord_ = 0.0;

    double miterLimit_ = 2.0;
    double arcTolerance_ = 0.0;
    bool preserveCollinear_ = false;
    bool reverseSolution_ = false;
    int stepCount_ = 0;
    double superellipseExp_ = 2.5;
    double endExtensionScale_ = 1.0;
    double arrowBackSweep_ = 0.0;
    double teardropPinch_ = 0.5;
    double joinAngleThreshold_ = 0.0;
    bool chebyshevSpacing_ = false;
    int errorCode_ = 0;
    double criticalAngleThreshold_ = 0.3; // radians
    // Segment fraction: 0 = default (0.15), 1 = disabled
    double criticalSegmentFraction_ = 0.0;
    // Simplification epsilon: 0 = default (auto-scaled)
    double simplifyEpsilon_ = 0.0;
    std::vector<double> critical_t_values_;
    // Start-point normalization: rotate each output path to start at
    // the lowest (most negative Y) point, with leftmost as tiebreaker.
    bool normalizeStart_ = false;

public:
    ClipperOffsetD(double miterLimit = 2.0, double arcTolerance = 0.0)
        : miterLimit_(miterLimit), arcTolerance_(arcTolerance) {}

    void Clear() {
        storedPaths_.clear();
        maxAbsCoord_ = 0.0;
        errorCode_ = 0;
        critical_t_values_.clear();
    }

    void AddPath(const PathD& path, JoinType jt, EndType et) {
        maxAbsCoord_ = std::max(maxAbsCoord_, GetMaxAbsCoord(path));
        storedPaths_.push_back({path, jt, et});
    }

    void AddPaths(const PathsD& paths, JoinType jt, EndType et) {
        for (const auto& path : paths) {
            AddPath(path, jt, et);
        }
    }

    PathsD Execute(double delta) {
        if (storedPaths_.empty()) return PathsD();

        double scale = ComputeOptimalScale(maxAbsCoord_, std::abs(delta));
        double invScale = 1.0 / scale;

        ClipperOffset co(miterLimit_, arcTolerance_ * scale);
        co.PreserveCollinear(preserveCollinear_);
        co.ReverseSolution(reverseSolution_);
        co.SetStepCount(stepCount_);
        co.SetSuperellipseExponent(superellipseExp_);
        co.SetEndCapParams(endExtensionScale_, arrowBackSweep_);
        co.SetTeardropPinch(teardropPinch_);
        co.SetJoinAngleThreshold(joinAngleThreshold_);
        co.SetChebyshevSpacing(chebyshevSpacing_);
        co.SetAngleThreshold(criticalAngleThreshold_);
        co.SetCriticalSegmentFraction(criticalSegmentFraction_);
        // Scale simplify epsilon if provided (it's in input coordinate units)
        const double simplifyEpsScaled = (simplifyEpsilon_ > 0.0) ?
            (simplifyEpsilon_ * scale) : 0.0;
        co.SetSimplifyEpsilon(simplifyEpsScaled);

        for (const auto& entry : storedPaths_) {
            co.AddPath(PathDToPath64(entry.path, scale), entry.joinType, entry.endType);
        }

        Paths64 result64;
        co.Execute(delta * scale, result64);
        errorCode_ = co.ErrorCode();
        critical_t_values_ = co.GetCriticalTValues();

        PathsD result = Paths64ToPathsD(result64, invScale);

        if (normalizeStart_) {
            NormalizePathsStart(result);
        }

        return result;
    }

    double MiterLimit() const { return miterLimit_; }
    void SetMiterLimit(double val) { miterLimit_ = val; }
    double ArcTolerance() const { return arcTolerance_; }
    void SetArcTolerance(double val) { arcTolerance_ = val; }
    bool PreserveCollinear() const { return preserveCollinear_; }
    void SetPreserveCollinear(bool val) { preserveCollinear_ = val; }
    bool ReverseSolution() const { return reverseSolution_; }
    void SetReverseSolution(bool val) { reverseSolution_ = val; }
    void SetStepCount(int val) { stepCount_ = val; }
    int StepCount() const { return stepCount_; }
    void SetSuperellipseExponent(double val) { superellipseExp_ = val; }
    double SuperellipseExponent() const { return superellipseExp_; }
    void SetEndCapParams(double ext, double sweep) { endExtensionScale_ = ext; arrowBackSweep_ = sweep; }
    double EndExtensionScale() const { return endExtensionScale_; }
    double ArrowBackSweep() const { return arrowBackSweep_; }
    void SetTeardropPinch(double val) { teardropPinch_ = val; }
    double TeardropPinch() const { return teardropPinch_; }
    void SetJoinAngleThreshold(double val) { joinAngleThreshold_ = val; }
    double JoinAngleThreshold() const { return joinAngleThreshold_; }
    void SetChebyshevSpacing(bool val) { chebyshevSpacing_ = val; }
    bool ChebyshevSpacing() const { return chebyshevSpacing_; }
    int ErrorCode() const { return errorCode_; }
    void SetAngleThreshold(double val) { criticalAngleThreshold_ = val; } // radians
    double AngleThreshold() const { return criticalAngleThreshold_; }
    // Segment fraction: 0 = default (0.15), 1 = disabled
    void SetCriticalSegmentFraction(double val) { criticalSegmentFraction_ = val; }
    double CriticalSegmentFraction() const { return criticalSegmentFraction_; }
    // Simplification epsilon: 0 = default (auto-scaled based on delta)
    void SetSimplifyEpsilon(double val) { simplifyEpsilon_ = val; }
    double SimplifyEpsilon() const { return simplifyEpsilon_; }
    std::vector<double> GetCriticalTValues() const { return critical_t_values_; }
    // Start-point normalization: rotate each output path to start at
    // the lowest (most negative Y) vertex, leftmost as tiebreaker.
    void SetNormalizeStart(bool val) { normalizeStart_ = val; }
    bool NormalizeStart() const { return normalizeStart_; }
};

ClipperOffsetD* CreateClipperOffsetD(double miterLimit, double arcTolerance) {
    return new ClipperOffsetD(miterLimit, arcTolerance);
}

/*
Clipper64* CreateClipper64(bool preserveCollinear) {
    Clipper64* clipper = new Clipper64();
    clipper->PreserveCollinear(preserveCollinear);
    return clipper;
}

ClipperD* CreateClipperD(bool preserveCollinear) {
    ClipperD* clipper = new ClipperD();
    clipper->PreserveCollinear(preserveCollinear);
    return clipper;
}
*/

template <typename T> intptr_t getVecDataPtr(std::vector<T> &vec) {
  return reinterpret_cast<intptr_t>(vec.data());
}

template <typename T>
emscripten::class_<std::vector<T>> register_vector_custom(const char *name) {
  typedef std::vector<T> VecType;

  void (VecType::*resize)(const size_t, const T &) = &VecType::resize;
  size_t (VecType::*size)() const = &VecType::size;
  return emscripten::class_<std::vector<T>>(name)
      .template constructor<>()
      .function("resize", resize)
      .function("size", size)
      .function("data", &getVecDataPtr<T>, emscripten::allow_raw_pointers());
}

EMSCRIPTEN_BINDINGS(clipper_module) {
        register_vector_custom<double>("VectorDouble");
        register_vector_custom<int>("VectorInt");
        /*
        class_<ClipperBase>("ClipperBase")
        .function("Clear", &ClipperBase::Clear)
        .function("SetPreserveCollinear", select_overload<void(bool)>(&ClipperBase::PreserveCollinear))
        .function("GetPreserveCollinear", select_overload<bool() const>(&ClipperBase::PreserveCollinear));
        */

        enum_<FillRule>("FillRule")
        .value("EvenOdd", FillRule::EvenOdd)
        .value("NonZero", FillRule::NonZero)
        .value("Positive", FillRule::Positive)
        .value("Negative", FillRule::Negative);

        /*
        enum_<ClipType>("ClipType")
        .value("Intersection", ClipType::Intersection)
        .value("Union", ClipType::Union)
        .value("Difference", ClipType::Difference)
        .value("Xor", ClipType::Xor);

        enum_<PathType>("PathType")
        .value("Subject", PathType::Subject)
        .value("Clip", PathType::Clip);
        */

        enum_<JoinType>("JoinType")
        .value("Square", JoinType::Square)
        .value("Bevel", JoinType::Bevel)
        .value("Round", JoinType::Round)
        .value("Miter", JoinType::Miter)
        .value("Superellipse", JoinType::Superellipse)
        .value("Knob", JoinType::Knob)
        .value("Step", JoinType::Step)
        .value("Spike", JoinType::Spike);

        enum_<EndType>("EndType")
        .value("Polygon", EndType::Polygon)
        .value("Joined", EndType::Joined)
        .value("Butt", EndType::Butt)
        .value("Square", EndType::Square)
        .value("Round", EndType::Round)
        .value("Superellipse", EndType::Superellipse)
        .value("Triangle", EndType::Triangle)
        .value("Arrow", EndType::Arrow)
        .value("Teardrop", EndType::Teardrop);

        /*
        enum_<PointInPolygonResult>("PointInPolygonResult")
        .value("IsOn", PointInPolygonResult::IsOn)
        .value("IsInside", PointInPolygonResult::IsInside)
        .value("IsOutside", PointInPolygonResult::IsOutside);
        */

        /*
        // #############################
        // ###### 64 bit bindings ######
        // #############################
        // NOTE: 64-bit integer path bindings are commented out.
        // Use the double-precision (PathD) bindings below instead.
        // Clipper2 handles precision scaling internally.

	// Point64 bindings (for now only support USINGZ=ON)
	#ifdef USINGZ
	class_<Point64>("Point64")
        .constructor<int64_t, int64_t, int64_t>()
        .property("x", &Point64::x)
        .property("y", &Point64::y)
        .property("z", &Point64::z)
        .function("SetZ", &Point64::SetZ);
	#endif

	// Path64
	class_<Path64>("Path64")
        .constructor<>()
        .function("size", &Path64::size)
        .function("clear", &Path64::clear)
        .function("push_back", select_overload<void(const Point64&)>(&Path64::push_back))
        .function("get", select_overload<Point64&(size_t)>(&Path64::operator[]), allow_raw_pointers());

	// Paths64
	class_<Paths64>("Paths64")
        .constructor<>()
        .function("size", &Paths64::size)
        .function("clear", &Paths64::clear)
        .function("push_back", select_overload<void(const Path64&)>(&Paths64::push_back))
        .function("get", select_overload<Path64&(size_t)>(&Paths64::operator[]), allow_raw_pointers());

        // Misc64
        function("AreaPath64", select_overload<double(const Path64&)>(&Area), allow_raw_pointers());
        function("AreaPaths64", select_overload<double(const Paths64&)>(&Area), allow_raw_pointers());
        function("IsPositive64", select_overload<bool(const Path64&)>(&IsPositive), allow_raw_pointers());
        // function("PointInPolygon64", select_overload<PointInPolygsonResult(const Point64&, const Path64&)>(&PointInPolygon), allow_raw_pointers());
        // function("ReversePath64", &ReversePath<int64_t>, allow_raw_pointers());

        // Geometry64
        class_<Rect64>("Rect64")
        .constructor<>()
        .constructor<int64_t, int64_t, int64_t, int64_t>()
        .property("left", &Rect64::left)
        .property("top", &Rect64::top)
        .property("right", &Rect64::right)
        .property("bottom", &Rect64::bottom)
        .function("IsValid", &Rect64::IsValid)
        .function("Width", select_overload<int64_t() const>(&Rect<int64_t>::Width))
        .function("Height", select_overload<int64_t() const>(&Rect<int64_t>::Height))
        .function("MidPoint", &Rect64::MidPoint)
        .function("AsPath", &Rect64::AsPath)
        .function("ContainsPoint", select_overload<bool(const Point<int64_t>&) const>(&Rect<int64_t>::Contains))
        .function("ContainsRect", select_overload<bool(const Rect<int64_t>&) const>(&Rect<int64_t>::Contains))
        .function("Scale", &Rect64::Scale)
        .function("IsEmpty", &Rect64::IsEmpty)
        .function("Intersects", &Rect64::Intersects)
        .function("Equals", &Rect64::operator==);

        // function("Ellipse64", select_overload<Path64(const Point64&, double, double, size_t)>(&Ellipse), allow_raw_pointers());
        // function("EllipseFromRect64", select_overload<Path64(const Rect64&, size_t)>(&Ellipse), allow_raw_pointers());

        // Translate64
        // function("TranslatePath64", select_overload<Path64(const Path64&, int64_t, int64_t)>(&TranslatePath), allow_raw_pointers());
        // function("TranslatePaths64", select_overload<Paths64(const Paths64&, int64_t, int64_t)>(&TranslatePaths), allow_raw_pointers());

        // RectClip64
        function("RectClipPaths64", select_overload<Paths64(const Rect64&, const Paths64&)>(&RectClip), allow_raw_pointers());
        function("RectClipPath64", select_overload<Paths64(const Rect64&, const Path64&)>(&RectClip), allow_raw_pointers());
        function("RectClipLinesPaths64", select_overload<Paths64(const Rect64&, const Paths64&)>(&RectClipLines), allow_raw_pointers());
        function("RectClipLinesPath64", select_overload<Paths64(const Rect64&, const Path64&)>(&RectClipLines), allow_raw_pointers());

        // Minkowski
        function("MinkowskiSum64", select_overload<Paths64(const Path64&, const Path64&, bool)>(&MinkowskiSum), allow_raw_pointers());
        function("MinkowskiDiff64", select_overload<Paths64(const Path64&, const Path64&, bool)>(&MinkowskiDiff), allow_raw_pointers());

        // BooleanOps
        function("BooleanOp64", select_overload<Paths64(ClipType, FillRule, const Paths64&, const Paths64&)>(&BooleanOp), allow_raw_pointers());
        function("BooleanOpOut64", select_overload<void(ClipType, FillRule, const Paths64&, const Paths64&, PolyTree64&)>(&BooleanOp), allow_raw_pointers());
        function("Intersect64", select_overload<Paths64(const Paths64&, const Paths64&, FillRule)>(&Intersect), allow_raw_pointers());
        function("Union64", select_overload<Paths64(const Paths64&, const Paths64&, FillRule)>(&Union), allow_raw_pointers());
        function("UnionSelf64", select_overload<Paths64(const Paths64&, FillRule)>(&Union), allow_raw_pointers());
        function("Difference64", select_overload<Paths64(const Paths64&, const Paths64&, FillRule)>(&Difference), allow_raw_pointers());
        function("Xor64", select_overload<Paths64(const Paths64&, const Paths64&, FillRule)>(&Xor), allow_raw_pointers());

        // Offset
        function("InflatePaths64", select_overload<Paths64(const Paths64&, double, JoinType, EndType, double, double)>(&InflatePaths), allow_raw_pointers());
        */

        // ClipperOffsetD - Double-precision offset with automatic scaling
        // This stores paths and computes optimal scale at Execute time
        class_<ClipperOffsetD>("ClipperOffsetD")
            .constructor<double, double>()
            .function("Clear", &ClipperOffsetD::Clear)
            .function("AddPath", &ClipperOffsetD::AddPath)
            .function("AddPaths", &ClipperOffsetD::AddPaths)
            .function("Execute", &ClipperOffsetD::Execute)
            .function("MiterLimit", &ClipperOffsetD::MiterLimit)
            .function("SetMiterLimit", &ClipperOffsetD::SetMiterLimit)
            .function("ArcTolerance", &ClipperOffsetD::ArcTolerance)
            .function("SetArcTolerance", &ClipperOffsetD::SetArcTolerance)
            .function("PreserveCollinear", &ClipperOffsetD::PreserveCollinear)
            .function("SetPreserveCollinear", &ClipperOffsetD::SetPreserveCollinear)
            .function("ReverseSolution", &ClipperOffsetD::ReverseSolution)
            .function("SetReverseSolution", &ClipperOffsetD::SetReverseSolution)
            .function("SetStepCount", &ClipperOffsetD::SetStepCount)
            .function("StepCount", &ClipperOffsetD::StepCount)
            .function("SetSuperellipseExponent", &ClipperOffsetD::SetSuperellipseExponent)
            .function("SuperellipseExponent", &ClipperOffsetD::SuperellipseExponent)
            .function("SetEndCapParams", &ClipperOffsetD::SetEndCapParams)
            .function("EndExtensionScale", &ClipperOffsetD::EndExtensionScale)
            .function("ArrowBackSweep", &ClipperOffsetD::ArrowBackSweep)
            .function("SetTeardropPinch", &ClipperOffsetD::SetTeardropPinch)
            .function("TeardropPinch", &ClipperOffsetD::TeardropPinch)
            .function("SetJoinAngleThreshold", &ClipperOffsetD::SetJoinAngleThreshold)
            .function("JoinAngleThreshold", &ClipperOffsetD::JoinAngleThreshold)
            .function("SetChebyshevSpacing", &ClipperOffsetD::SetChebyshevSpacing)
            .function("ChebyshevSpacing", &ClipperOffsetD::ChebyshevSpacing)
            .function("ErrorCode", &ClipperOffsetD::ErrorCode)
            .function("SetAngleThreshold", &ClipperOffsetD::SetAngleThreshold)
            .function("AngleThreshold", &ClipperOffsetD::AngleThreshold)
            .function("SetCriticalSegmentFraction", &ClipperOffsetD::SetCriticalSegmentFraction)
            .function("CriticalSegmentFraction", &ClipperOffsetD::CriticalSegmentFraction)
            .function("SetSimplifyEpsilon", &ClipperOffsetD::SetSimplifyEpsilon)
            .function("SimplifyEpsilon", &ClipperOffsetD::SimplifyEpsilon)
            .function("GetCriticalTValues", &ClipperOffsetD::GetCriticalTValues)
            .function("SetNormalizeStart", &ClipperOffsetD::SetNormalizeStart)
            .function("NormalizeStart", &ClipperOffsetD::NormalizeStart);

        function("CreateClipperOffsetD", &CreateClipperOffsetD, allow_raw_pointers());

        // Simple auto-scaling offset function (recommended for most use cases)
        function("InflatePathsD", &InflatePathsDAutoScale, allow_raw_pointers());

        // Auto-scaling boolean operations for PathsD
        function("UnionD", &UnionDAutoScale, allow_raw_pointers());
        function("UnionSelfD", &UnionSelfDAutoScale, allow_raw_pointers());
        function("IntersectD", &IntersectDAutoScale, allow_raw_pointers());
        function("DifferenceD", &DifferenceDAutoScale, allow_raw_pointers());
        function("XorD", &XorDAutoScale, allow_raw_pointers());

        /*
        // Simplify64 - commented out, use SimplifyPathD/SimplifyPathsD instead
        function("SimplifyPath64", select_overload<Path64(const Path64&, double, bool)>(&SimplifyPath), allow_raw_pointers());
        function("SimplifyPaths64", select_overload<Paths64(const Paths64&, double, bool)>(&SimplifyPaths), allow_raw_pointers());
        function("TrimCollinear64", select_overload<Path64(const Path64&, bool)>(&TrimCollinear), allow_raw_pointers());
        */

        // PolyPath
        /*
        class_<PolyPath>("PolyPath")
        .function("isHole", &PolyPath::IsHole);

        // PolyPath64
        class_<PolyPath64, base<PolyPath>>("PolyPath64")
        .constructor<>()
        .function("addChild", &PolyPath64::AddChild, allow_raw_pointers())
        .function("clear", &PolyPath64::Clear)
        .function("count", &PolyPath64::Count)
        .function("polygon", &PolyPath64::Polygon)
        .function("area", &PolyPath64::Area)
        .function("child", &PolyPath64::Child, allow_raw_pointers());
        */

        // Clipper64
        /*
        class_<Clipper64, base<ClipperBase>>("Clipper64")
        .constructor<>()
        .function("AddSubject", &Clipper64::AddSubject, allow_raw_pointers())
        .function("AddOpenSubject", &Clipper64::AddOpenSubject, allow_raw_pointers())
        .function("AddClip", &Clipper64::AddClip, allow_raw_pointers())
        .function("Clear", &Clipper64::Clear)
        .function("ExecutePath", select_overload<bool(ClipType, FillRule, Paths64&)>(&Clipper64::Execute), allow_raw_pointers())
        .function("ExecutePath", select_overload<bool(ClipType, FillRule, Paths64&, Paths64&)>(&Clipper64::Execute), allow_raw_pointers())
        .function("ExecutePoly", select_overload<bool(ClipType, FillRule, PolyTree64&)>(&Clipper64::Execute), allow_raw_pointers())
        .function("ExecutePoly", select_overload<bool(ClipType, FillRule, PolyTree64&, Paths64&)>(&Clipper64::Execute), allow_raw_pointers());

        function("CreateClipper64", &CreateClipper64, allow_raw_pointers());
        */

        // #############################
        // ###### Double precision bindings ######
        // #############################

        // PointD bindings (for now only support USINGZ=ON)
	#ifdef USINGZ
	class_<PointD>("PointD")
        .constructor<double, double, double>()
        .property("x", &PointD::x)
        .property("y", &PointD::y)
        .property("z", &PointD::z)
        .function("SetZ", &PointD::SetZ);
	#endif

	// PathD
	class_<PathD>("PathD")
        .constructor<>()
        .function("size", &PathD::size)
        .function("clear", &PathD::clear)
        .function("push_back", select_overload<void(const PointD&)>(&PathD::push_back))
        .function("get", select_overload<PointD&(size_t)>(&PathD::operator[]), allow_raw_pointers());

	// PathsD
	class_<PathsD>("PathsD")
		.constructor<>()
		.function("size", &PathsD::size)
		.function("clear", &PathsD::clear)
		.function("push_back", select_overload<void(const PathD&)>(&PathsD::push_back))
		.function("get", select_overload<PathD&(size_t)>(&PathsD::operator[]), allow_raw_pointers());

        // MiscD
        // function("AreaPathD", select_overload<double(const PathD&)>(&Area), allow_raw_pointers());
        // function("AreaPathsD", select_overload<double(const PathsD&)>(&Area), allow_raw_pointers());
        // function("IsPositiveD", select_overload<bool(const PathD&)>(&IsPositive), allow_raw_pointers());
        // function("PointInPolygonD", select_overload<PointInPolygonResult(const PointD&, const PathD&)>(&PointInPolygon), allow_raw_pointers());
        // function("ReversePathD", &ReversePath<double>, allow_raw_pointers());

        // GeometryD
        class_<RectD>("RectD")
        .constructor<>()
        .constructor<double, double, double, double>()
        .property("left", &RectD::left)
        .property("top", &RectD::top)
        .property("right", &RectD::right)
        .property("bottom", &RectD::bottom)
        .function("IsValid", &RectD::IsValid)
        .function("Width", select_overload<double() const>(&Rect<double>::Width))
        .function("Height", select_overload<double() const>(&Rect<double>::Height))
        .function("MidPoint", &RectD::MidPoint)
        .function("AsPath", &RectD::AsPath)
        .function("ContainsPoint", select_overload<bool(const Point<double>&) const>(&Rect<double>::Contains))
        .function("ContainsRect", select_overload<bool(const Rect<double>&) const>(&Rect<double>::Contains))
        .function("Scale", &RectD::Scale)
        .function("IsEmpty", &RectD::IsEmpty)
        .function("Intersects", &RectD::Intersects)
        .function("Equals", &RectD::operator==);

        // function("EllipseD", select_overload<PathD(const PointD&, double, double, size_t)>(&Ellipse), allow_raw_pointers());
        // function("EllipseFromRectD", select_overload<PathD(const RectD&, size_t)>(&Ellipse), allow_raw_pointers());

        // TranslateD
        // function("TranslatePathD", select_overload<PathD(const PathD&, double, double)>(&TranslatePath), allow_raw_pointers());
        // function("TranslatePathsD", select_overload<PathsD(const PathsD&, double, double)>(&TranslatePaths), allow_raw_pointers());

        // RectClipD
        /*
        function("RectClipPathsD", select_overload<PathsD(const RectD&, const PathsD&, int)>(&RectClip), allow_raw_pointers());
        function("RectClipPathD", select_overload<PathsD(const RectD&, const PathD&, int)>(&RectClip), allow_raw_pointers());
        function("RectClipLinesPathsD", select_overload<PathsD(const RectD&, const PathsD&, int)>(&RectClipLines), allow_raw_pointers());
        function("RectClipLinesPathD", select_overload<PathsD(const RectD&, const PathD&, int)>(&RectClipLines), allow_raw_pointers());
        */

        // Minkowski
        /*
        function("MinkowskiSumD", select_overload<PathsD(const PathD&, const PathD&, bool, int)>(&MinkowskiSum), allow_raw_pointers());
        function("MinkowskiDiffD", select_overload<PathsD(const PathD&, const PathD&, bool, int)>(&MinkowskiDiff), allow_raw_pointers());
        */

        // BooleanOps
        /*
        function("BooleanOpD", select_overload<PathsD(ClipType, FillRule, const PathsD&, const PathsD&, int)>(&BooleanOp), allow_raw_pointers());
        function("BooleanOpOutD", select_overload<void(ClipType, FillRule, const PathsD&, const PathsD&, PolyTreeD&, int)>(&BooleanOp), allow_raw_pointers());
        function("IntersectD", select_overload<PathsD(const PathsD&, const PathsD&, FillRule, int)>(&Intersect), allow_raw_pointers());
        function("UnionD", select_overload<PathsD(const PathsD&, const PathsD&, FillRule, int)>(&Union), allow_raw_pointers());
        function("UnionSelfD", select_overload<PathsD(const PathsD&, FillRule, int)>(&Union), allow_raw_pointers());
        function("DifferenceD", select_overload<PathsD(const PathsD&, const PathsD&, FillRule, int)>(&Difference), allow_raw_pointers());
        function("XorD", select_overload<PathsD(const PathsD&, const PathsD&, FillRule, int)>(&Xor), allow_raw_pointers());
        */

        // Offset - using auto-scaling InflatePathsD defined above
        // (The built-in InflatePaths with manual precision is not exposed)

        // Simplify
        function("SimplifyPathD", select_overload<PathD(const PathD&, double, bool)>(&SimplifyPath), allow_raw_pointers());
        function("SimplifyPathsD", select_overload<PathsD(const PathsD&, double, bool)>(&SimplifyPaths), allow_raw_pointers());
        function("TrimCollinearD", select_overload<PathD(const PathD&, int, bool)>(&TrimCollinear), allow_raw_pointers());

        // PolyPathD
        /*
        class_<PolyPathD, base<PolyPath>>("PolyPathD")
        .constructor<>()
        .function("addChild", select_overload<PolyPathD*(const PathD&)>(&PolyPathD::AddChild), allow_raw_pointers())
        .function("clear", &PolyPathD::Clear)
        .function("count", &PolyPathD::Count)
        .function("polygon", &PolyPathD::Polygon)
        .function("area", &PolyPathD::Area)
        .function("child", &PolyPathD::Child, allow_raw_pointers());
        */

        // ClipperD
        /*
        class_<ClipperD, base<ClipperBase>>("ClipperD")
        .constructor<int>()
        .function("AddSubject", &ClipperD::AddSubject, allow_raw_pointers())
        .function("AddOpenSubject", &ClipperD::AddOpenSubject, allow_raw_pointers())
        .function("AddClip", &ClipperD::AddClip, allow_raw_pointers())
        .function("Clear", &ClipperD::Clear)
        .function("ExecutePath", select_overload<bool(ClipType, FillRule, PathsD&)>(&ClipperD::Execute), allow_raw_pointers())
        .function("ExecutePath", select_overload<bool(ClipType, FillRule, PathsD&, PathsD&)>(&ClipperD::Execute), allow_raw_pointers())
        .function("ExecutePoly", select_overload<bool(ClipType, FillRule, PolyTreeD&)>(&ClipperD::Execute), allow_raw_pointers())
        .function("ExecutePoly", select_overload<bool(ClipType, FillRule, PolyTreeD&, PathsD&)>(&ClipperD::Execute), allow_raw_pointers());

        function("CreateClipperD", &CreateClipperD, allow_raw_pointers());
        */
}
