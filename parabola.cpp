#include <cmath>
#include <iostream>
#include <string>

#include <format>
#include <vector>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "variable.h"

double eval(
    const std::string &expr,
    const Evaluator::library_t &lib,
    double x)
{
    dynamic_cast<Variable &>(*lib.at("x")).value = x;
    return Evaluator::eval(expr, lib);
}

void plot_point(const double x, const double y)
{
    std::cerr << "(" << x << ", " << y << ")" << std::endl;
    std::cout << x << " " << y << std::endl;
}

void point_push(std::vector<double> &points, double x, double y)
{
    points.push_back(x);
    points.push_back(y);
}

double get_top(const double x0, const double y0,
               const double x1, const double y1,
               const double x2, const double y2)
{
    const double a = (y2 - (x2*(y1 - y0) + x1*y0 - x0*y1)/(x1 - x0))/(x2*(x2-x0-x1) + x0*x1);
    std::cerr << "a = " << a << std::endl;
    const double b = (y1 - y0)/(x1 - x0) - a*(x0 + x1);
    std::cerr << "b = " << b << std::endl;
    return -b/a/2.0;
}

std::vector<double> find_min(
    const std::string &expr,
    const Evaluator::library_t &lib,
    double a, double b, double x_start, double h)
{
    std::vector<double> points;
    double dir = 1.0;

    double x0 = x_start;
    double y0 = eval(expr, lib, x0);
    point_push(points, x0, y0);

    int i = 1;
    double x1 = x_start + dir*i*h;
    double y1 = eval(expr, lib, x1);
    point_push(points, x1, y1);

    if (y1 > y0) {
        dir = -1.0;
        x1 = x_start + dir*i*h;
        y1 = eval(expr, lib, x1);
        point_push(points, x1, y1);
    }

    if (y1 > y0) {
        // std::cerr << "Already minimum" << std::endl;
        return points;
    }

    double x2 = 0.0;
    double y2 = 0.0;

    i = 2;
    while (true) {
        x2 = x_start + dir*h*std::pow(2, i-1);

        const bool isOutOfBounds = (x2 < a || x2 > b);
        if (isOutOfBounds) {
            std::cerr << "Out of bounds" << std::endl;
            double z = get_top(x0, y0, x1, y1, x2, y2);
            double v = eval(expr, lib, z);
            x2 = std::abs(x2 - a) < std::abs(x2 - b) ? a : b;
            y2 = eval(expr, lib, x2);
            if (v < y2) {
                std::cerr << "(but not really)" << std::endl;
                point_push(points, z, v);
            }
            point_push(points, x2, y2);
            break;
        }

        y2 = eval(expr, lib, x2);
        point_push(points, x2, y2);

        double delta_n = y0 - y1;
        double delta_p = y2 - y1;
        const bool isConvex = (delta_n >= 0) && (delta_p >= 0) && (delta_n + delta_p > 0);
        if (isConvex) {
            std::cerr << "Convex" << std::endl;
            // x2 = x1 + dir * h/2 * (delta_n - delta_p) / (delta_n + delta_p);
            x2 = get_top(x0, y0, x1, y1, x2, y2);
            y2 = eval(expr, lib, x2);
            point_push(points, x2, y2);
            break;
        }

        x0 = x1; x1 = x2;
        y0 = y1; y1 = y2;
        ++i;
    }

    return points;
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

    // {
    //     double x0 = -3.0;
    //     double y0 = 0.0;
    //     double x1 = -1.5;
    //     double y1 = -0.75;
    //     double x2 = 0.0;
    //     double y2 = 3.0;

    //     double delta_n = y0 - y1;
    //     double delta_p = y2 - y1;
    //     const bool isConvex = (delta_n >= 0) && (delta_p >= 0) && (delta_n + delta_p > 0);
    //     if (isConvex) {
    //         std::cerr << "Convex" << std::endl;
    //         std::cerr << get_top(x0, y0, x1, y1, x2, y2) << std::endl;
    //     }
    //     else
    //     {
    //         std::cerr << "Not convex" << std::endl;
    //     }
    // }
    // return 0;

    /* == Plotting pipeline ==
     *
     * plot f(x) ls 2, '-' ls 1, p(x) ls 3
     *
     * where f(x) - function
     */

    const double eps = (argc > 1 + 3) ? std::strtof(argv[4], nullptr) : 1e-15;

    std::vector<double> points;
    double h = std::abs(b - a) / 1000;

    points = find_min(
        expression, library,
        a, b, (a + b)/2, h);

    const double &ymin = points[points.size() - 1];
    const double &xmin = points[points.size() - 2];

    std::cout << "$points << e" << std::endl;
    for (int i = 0; i < points.size() - 2; i += 2)
        std::cout << points[i] << " " << points[i+1] << std::endl;
    std::cout << "e" << std::endl;

    std::cout << "$min << e" << std::endl;
    std::cout << xmin << " " << ymin << std::endl;
    std::cout << "e" << std::endl;

    const double &y3 = points[points.size() - 3];
    const double &x3 = points[points.size() - 4];
    const double &y2 = points[points.size() - 5];
    const double &x2 = points[points.size() - 6];
    const double &y1 = points[points.size() - 7];
    const double &x1 = points[points.size() - 8];

    std::cout << "plot " << gnuplot_convert(expression) << " ls 2";
    std::cout << std::format(", {} * (x - {})*(x - {}) / ({} - {})/({} - {})",
                             y1, x2, x3, x1, x2, x1, x3);
    std::cout << std::format("+ {} * (x - {})*(x - {}) / ({} - {})/({} - {})",
                             y2, x1, x3, x2, x1, x2, x3);
    std::cout << std::format("+ {} * (x - {})*(x - {}) / ({} - {})/({} - {})",
                             y3, x1, x2, x3, x1, x3, x2);
    std::cout << "ls 3 title \"parabola\"";
    std::cout << ", $points with linespoints ls 3 title \"points\", $min with points ls 1 title \"min\"" << std::endl;

    int trash;
    std::cin >> trash;

    return 0;
}
