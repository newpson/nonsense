#include <cmath>
#include <iostream>
#include <string>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "math.h"
#include "variable.h"

double eval(
    const std::string &expr,
    const Evaluator::library_t &lib,
    double x)
{
    dynamic_cast<Variable &>(*lib.at("x")).value = x;
    return Evaluator::eval(expr, lib);
}

struct Parabola
{
    double a;
    double b;
    double c;

    static double calculate_a(const Vector2D &v1, const Vector2D &v2, const Vector2D &v3)
    {
        return (v3.y - (v3.x*(v2.y - v1.y) + v2.x*v1.y - v1.x*v2.y) / (v2.x - v1.x))
               / (v3.x*(v3.x-v1.x-v2.x) + v1.x*v2.x);
    }

    static double calculate_b(const Vector2D &v1, const Vector2D &v2, const double a)
    {
        return (v2.y - v1.y) / (v2.x - v1.x) - a*(v1.x + v2.x);
    }

    static double calculate_c(const Vector2D &v1, const Vector2D &v2, const double a)
    {
        return (v2.x*v1.y - v1.x*v2.y) / (v2.x - v1.x) + a * v1.x * v2.x;
    }

    static Vector2D get_top(const double a, const double b, const double c)
    {
        const double x = -b/2.0/a;
        const double y = a*x*x + b*x + c;
        return Vector2D(x, y);
    }

    static Vector2D get_top(const Vector2D &v1, const Vector2D &v2, const Vector2D &v3)
    {
        const double a = calculate_a(v1, v2, v3);
        return get_top(a, calculate_b(v1, v2, a), calculate_c(v1, v2, a));
    }

    Vector2D get_top() const
    {
        return get_top(a, b, c);
    }

    Parabola(const double a = 0.0, const double b = 0.0, const double c = 0.0)
        : a(a), b(b), c(c)
    {}

    Parabola(const Vector2D &v1, const Vector2D &v2, const Vector2D &v3)
        : a(calculate_a(v1, v2, v3))
        , b(calculate_b(v1, v2, a))
        , c(calculate_c(v1, v2, a))
    {}

    friend std::ostream &operator<<(std::ostream &out, const Parabola &p)
    {
        out << p.a << "*x**2 + " << p.b << "*x + " << p.c;
        return out;
    }
};

Vector2D find_min(
    const std::string &expr,
    const Evaluator::library_t &lib,
    const double a, const double b, const double origin,
    double h)
{
    if (b < a) {
        std::cerr << "invalid bounds order" << std::endl;
        return Vector2D();
    }
    if (origin > b || origin < a) {
        std::cerr << "origin out of bounds" << std::endl;
        return Vector2D();
    }

    const Vector2D A(a, eval(expr, lib, a));
    const Vector2D B(b, eval(expr, lib, b));
    const Vector2D O(origin, eval(expr, lib, origin));

    Vector2D cur(origin + h, eval(expr, lib, origin + h));
    if (!is_inside(cur.x, a, b)) {
        std::cerr << "early out of bounds" << std::endl;
        const double px = Parabola::get_top(A, O, B).x;
        const Vector2D ptop(px, eval(expr, lib, px));
        return is_inside(ptop.x, a, b) ? min(A, O, B, ptop)
                                    : min(A, O, B);
    }

    if (cur.y > O.y) {
        std::cerr << "wrong direction" << std::endl;
        h = -h;
        cur = Vector2D(origin + h, eval(expr, lib, origin + h));
        if (cur.y > O.y) {
            std::cerr << "early convex" << std::endl;
            const Vector2D positive_hO(origin + h, eval(expr, lib, origin + h));
            const Vector2D negative_hO(origin - h, eval(expr, lib, origin - h));
            const double px = Parabola::get_top(positive_hO, O, negative_hO).x;
            const Vector2D ptop(px, eval(expr, lib, px));
            return is_inside(ptop.x, a, b) ? min(positive_hO, O, negative_hO, ptop)
                                        : min(positive_hO, O, negative_hO);
        }
    }

    Vector2D prev = O;
    Vector2D next;
    int i = 1;
    while (true) {
        next.x = origin + h*std::pow(2, i);
        next.y = eval(expr, lib, next.x);

        if (!is_inside(next.x, a, b)) {
            std::cerr << "out of bounds" << std::endl;
            const double px = Parabola::get_top(prev, cur, next).x;
            const Vector2D ptop(px, eval(expr, lib, px));
            return is_inside(ptop.x, a, b) ? min(A, B, cur, ptop)
                                           : min(A, B, cur);
        }

        if (is_convex(prev, cur, next)) {
            std::cerr << "convex" << std::endl;
            Parabola parabola(prev, cur, next);
            std::cerr << parabola << std::endl;
            const double px = parabola.get_top().x;
            return Vector2D(px, eval(expr, lib, px));
        }

        prev = cur;
        cur = next;
        ++i;
    }

    return Vector2D();
}

int main(const int argc, const char **argv)
{
    if (argc < 1 + 3)
    {
        std::cout << "<expression> <a> <b> [eps]" << std::endl;
        return 0;
    }

    std::shared_ptr<Variable> x(new Variable());
    Evaluator::library_t library = {
        {"x", x},
    };
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    const std::string expression(argv[1]);
    const double a = std::strtof(argv[2], nullptr);
    const double b = std::strtof(argv[3], nullptr);

    // TODO make iterative calculation of minimum until eps
    const double eps = (argc > 1 + 3) ? std::strtof(argv[4], nullptr) : 1e-15;

    double h = std::abs(b - a) / 1.5;

    Vector2D minimum = find_min(
        expression, library,
        a, b, (a + b)/2, h);

    std::cerr << minimum << std::endl;

    return 0;
}
