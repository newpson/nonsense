#include <ostream>

struct Vector2D
{
    double x;
    double y;

    Vector2D(const double x = 0.0, const double y = 0.0): x(x), y(y) {}
    Vector2D operator+(const Vector2D v) const { return Vector2D(x + v.x, y + v.y); }
    Vector2D operator-(const Vector2D v) const { return Vector2D(x - v.x, y - v.y); }

    friend std::ostream &operator<<(std::ostream &out, const Vector2D &v)
    {
        out << v.x << " " << v.y;
        return out;
    }
};
