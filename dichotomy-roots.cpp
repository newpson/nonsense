#include <cmath>
#include <iostream>
#include <set>
#include <string>
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

std::vector<std::pair<double, double>> search_intervals(
    const std::string &expr, Evaluator::library_t &lib,
    const double a, const double b,
    const double h = 1e-01, const double zero = 1e-06)
{
    std::vector<std::pair<double, double>> intervals;

    double x1 = a;
    double x2;
    double y1;
    double y2;
    do {
        y1 = eval(expr, lib, x1);
        x2 = x1 + h;
        y2 = eval(expr, lib, x2);
        if (y1 * y2 <= 0.0 || std::abs((y1 + y2) / 2.0) <= zero) {
            if (!intervals.empty() && intervals.back().second == x1)
                intervals.back().second = x2;
            else
                intervals.push_back({x1, x2});
        }
        x1 += h;
    } while (x2 < b);

    return intervals;
}

std::set<double> search_roots(
    const std::string &expr, Evaluator::library_t &lib,
    const std::vector<std::pair<double, double>> intervals,
    const double eps = 1e-15, const double zero = 1e-06)
{
    std::set<double> roots;

    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0, x3 = 0.0, y3 = 0.0;
    for (auto interval : intervals) {
        x1 = interval.first;
        x2 = interval.second;

        while (std::abs(x2 - x1) > 2 * eps) {
            y1 = eval(expr, lib, x1);
            y2 = eval(expr, lib, x2);

            x3 = (x1 + x2) / 2.0;
            y3 = eval(expr, lib, x3);

            if (y1 * y3 < y2 * y3)
                x2 = x3;
            else
                x1 = x3;
        }

        x3 = (x1 + x2) / 2.0;
        if (y1 * y2 < 0 || eval(expr, lib, x3) <= zero)
            roots.insert(x3);
    }

    return roots;
}

std::set<double> find_roots(
    const std::string &expr, Evaluator::library_t &lib,
    const std::shared_ptr<Variable> &x,
    const double a, const double b, const double eps = 1e-15,
    const double zero = 1e-06,
    const double split_h = 1e-02,
    const double split_eps = 1e-02)
{

    auto root_intervals = search_intervals(
        expr, lib, a, b, split_h, split_eps);

    return search_roots(
        expr, lib, root_intervals, eps, zero);
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

    std::set<double> roots;
    if (argc > 1 + 3)
    {
        const double eps = std::strtof(argv[4], nullptr);
        roots = find_roots(
            expression, library, x,
            a, b, eps);
    }
    else
        roots = find_roots(
            expression, library, x,
            a, b);

    for (double root: roots) {
        std::cerr << root << std::endl;
        std::cout << "set object circle at first " << root << ",0 radius char 0.25 fs solid fc rgb '#E41A1C' front" << std::endl;
    }

    std::cout << "plot " << gnuplot_convert(expression) << " ls 2" << std::endl;

    int trash;
    std::cin >> trash;

    return 0;
}
