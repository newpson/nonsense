#include <ostream>

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

    bool operator==(const Vector2D &v) const { return x == v.x && y == v.y; }

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
