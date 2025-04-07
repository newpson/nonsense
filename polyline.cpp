#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "variable.h"

struct Vector2D
{
    double x;
    double y;

    Vector2D(const double x = 0.0, const double y = 0.0): x(x), y(y) {}
    Vector2D operator+(const Vector2D v) { return Vector2D(x + v.x, y + v.y); }
    Vector2D operator-(const Vector2D v) { return Vector2D(x - v.x, y - v.y); }

    friend std::ostream &operator<<(std::ostream &out, const Vector2D &v)
    {
        out << v.x << " " << v.y;
        return out;
    }
};

double eval(
    const std::string &expr,
    const Evaluator::library_t &lib,
    const double t)
{
    dynamic_cast<Variable &>(*lib.at("t")).value = t;
    return Evaluator::eval(expr, lib);
}

std::vector<Vector2D> subdivide(const std::string &expr_x,
                                const std::string &expr_y,
                                const Evaluator::library_t &lib,
                                const bool close,
                                const double ta = 0.0, const double tb = 1.0, const double h = 0.01)
{
    std::vector<Vector2D> points;
    for (double t = ta; t <= tb; t += h) {
        points.push_back(Vector2D(eval(expr_x, lib, t), eval(expr_y, lib, t)));
    }
    if (close)
        points.push_back(points.front());
    return points;
}

int winding_number(const Vector2D &target, const std::vector<Vector2D> &points)
{
    int w = 0;
    for (int i = 0; i < points.size() - 1; ++i) {
        if (points[i].x >= target.x) {
            if (points[i].y <= target.y && target.y < points[i+1].y)
                ++w;
            else if (points[i].y > target.y && target.y >= points[i+1].y)
                --w;
        }
    }
    return w;
}

void plot_line(Vector2D &begin, Vector2D &end)
{
    const Vector2D h(0.5, 0.5);
    const Vector2D linevector = end - begin;
    const Vector2D direction(std::copysign(1.0, linevector.x),
                             std::copysign(1.0, linevector.y));
    const double invertor = copysign(1.0, linevector.x * linevector.y); // инверсия осей в чётных квадрантах

    const Vector2D cell_begin(direction.x > 0 ? std::floor(begin.x/h.x)*h.x : std::ceil(begin.x/h.x)*h.x,
                              direction.y > 0 ? std::floor(begin.y/h.y)*h.y : std::ceil(begin.y/h.y)*h.y);
    const Vector2D cell_end(direction.x > 0 ? std::floor(end.x/h.x)*h.x : std::ceil(end.x/h.x)*h.x,
                            direction.y > 0 ? std::floor(end.y/h.y)*h.y : std::ceil(end.y/h.y)*h.y);

    int i = 0;
    Vector2D cell = cell_begin;
    while ((std::abs(cell.x - cell_end.x) > h.x/1000.0
            || std::abs(cell.y - cell_end.y) > h.y/1000.0) && i < 1000) {
        std::cout << "set object rect from "
                  << cell.x << "," << cell.y
                  << " to "
                  << cell.x + direction.x*h.x << "," << cell.y + direction.y*h.y
                  << std::endl;
        cell.x += h.x*direction.x;
        cell.y += h.y*direction.y;
        const double D = invertor*((cell.y - begin.y)*linevector.x - (cell.x - begin.x)*linevector.y);
        if (D > 0)
            cell.y -= h.y*direction.y;
        else if (D < 0)
            cell.x -= h.x*direction.x;
        ++i;
    }
    std::cout << "set object rect from "
              << cell.x << "," << cell.y
              << " to "
              << cell.x + direction.x*h.x << "," << cell.y + direction.y*h.y
              << std::endl;

    // if (i >= 1000) {
    //     std::cerr << "*EXTREME* Overflow at (" << x1 << ", " << y1 << ")->(" << x2 << ", " << y2 << ")" << std::endl;
    //     std::cerr << "Last should be (" << endx << ", " << endy << ")" << std::endl;
    // }
}

int main(int argc, char **argv)
{
    std::shared_ptr<Variable> t(new Variable());
    Evaluator::library_t library = {
        {"t", t},
    };
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    auto path = subdivide("t*sin(t)", "t*cos(t)", library, false, 0.0, 5*3.14);
    // auto frame = subdivide("12*cos(t)", "4*sin(t)", library, true, 0.0, 2*3.15);
    std::vector<Vector2D> frame = {{0.0, 0.0}, {20.0, 0.0}, {20.0, 20.0}, {0.0, 20.0}, {0.0}, {0.0}};

    for (int i = 0; i < path.size() - 1; ++i) {
        const bool is_inside = (winding_number(path[i], frame) != 0 &&
                                winding_number(path[i+1], frame) != 0);
        if (is_inside)
            plot_line(path[i], path[i+1]);
    }

    std::cout << "$path << e" << std::endl;
    for (int i = 0; i < path.size(); ++i)
        std::cout << path[i] << std::endl;
    std::cout << "e" << std::endl;
    std::cout << "$frame << e" << std::endl;
    for (int i = 0; i < frame.size(); ++i)
        std::cout << frame[i] << std::endl;
    std::cout << "e" << std::endl;
    std::cout << "plot $path with lines ls 1, $frame with lines ls 2" << std::endl;
    int trash;
    std::cin >> trash;

    return 0;
}

