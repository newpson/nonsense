#include <ostream>
#include <cmath>

double floor_eps(const double value, const double eps = 1e-12)
{
    const double floored = std::floor(value);
    const double ceiled = std::ceil(value);
    if (std::abs(ceiled - value) < eps)
        return ceiled;
    return floored;
}

bool are_equal(const double a, const double b, const double eps = 1e-12)
{
    return std::abs(a - b) < eps;
}

struct Vector2D
{
    double x;
    double y;

    Vector2D(const double x = 0.0, const double y = 0.0): x(x), y(y) {}
    Vector2D(const Vector2D &v): x(v.x), y(v.y) {}
    Vector2D operator+(const Vector2D &v) const { return Vector2D(x + v.x, y + v.y); }
    Vector2D operator-(const Vector2D &v) const { return Vector2D(x - v.x, y - v.y); }

    // fancy vector operators
    Vector2D operator*(const Vector2D &v) const { return Vector2D(x * v.x, y * v.y); }
    Vector2D operator/(const Vector2D &v) const { return Vector2D(x / v.x, y / v.y); }
    Vector2D &operator=(const Vector2D &v) { x = v.x; y = v.y; return *this; }
    Vector2D floor() const { return Vector2D(floor_eps(x), floor_eps(y)); }
    Vector2D floor(const Vector2D &grid) const { return Vector2D(floor_eps(x/grid.x), floor_eps(y/grid.y)); }
    double length() const { return std::sqrt(x*x + y*y); }

    bool operator==(const Vector2D &v) const { return are_equal(x, v.x) && are_equal(y, v.y); }
    bool operator!=(const Vector2D &v) const { return !(*this == v); }

    friend std::ostream &operator<<(std::ostream &out, const Vector2D &v)
    {
        out << v.x << " " << v.y;
        return out;
    }
};

inline double min(const double a, const double b)
{
    return a < b ? a : b;
}

inline double min(const double a, const double b, const double c)
{
    return min(min(a, b), c);
}

inline double min(const double a, const double b, const double c, const double d)
{
    return min(min(min(a, b), c), d);
}

inline double is_inside(const double x, const double a, const double b)
{
    return (x >= a && x <= b);
}

inline const Vector2D &min(const Vector2D &a, const Vector2D &b)
{
    return a.y < b.y ? a : b;
}

inline const Vector2D &min(const Vector2D &a, const Vector2D &b, const Vector2D &c)
{
    return min(min(a, b), c);
}

inline const Vector2D &min(const Vector2D &a, const Vector2D &b, const Vector2D &c, const Vector2D &d)
{
    return min(min(min(a, b), c), d);
}

bool is_convex(const Vector2D &a, const Vector2D &center, const Vector2D &b)
{
    const double delta_ac = a.y - center.y;
    const double delta_bc = b.y - center.y;
    return (delta_ac >= 0) && (delta_bc >= 0) && (delta_ac + delta_bc > 0);
}
