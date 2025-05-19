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

Vector2D find_min(
    const std::string &expr,
    const Evaluator::library_t &lib,
    const double x0, const double eps = 1e-12)
{
    const double h = 0.1; // TODO dynamic step depending on epsilon (?)
    double x = x0;
    double err = 0.0;
    int i = 0;
    std::cerr << "info: running with x0 = " << x0 << ", eps = " << eps << std::endl;
    do {
        const double y_left = eval(expr, lib, x - h);
        const double y_mid = eval(expr, lib, x);
        const double y_right = eval(expr, lib, x + h);
        const double x_next = x - h/2.0 * (y_right - y_left) / (y_left - 2.0*y_mid + y_right);
        err = std::abs(x_next - x);
        x = x_next;
        ++i;
    } while (err > eps && i < 30);
    std::cerr << "info: i = " << i << std::endl;

    return Vector2D(x, eval(expr, lib, x));
}

int main(const int argc, const char **argv)
{
    if (argc < 1 + 2) {
        std::cout << "<expression> <x0> [eps]" << std::endl;
        return 0;
    }

    std::shared_ptr<Variable> x(new Variable());
    Evaluator::library_t library = {{"x", x}};
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    const std::string expression(argv[1]);
    const double x0 = std::strtod(argv[2], nullptr);
    Vector2D minimum;
    if (argc > 1 + 2) {
        const double eps = std::strtod(argv[3], nullptr);
        minimum = find_min(expression, library, x0, eps);
    } else {
        minimum = find_min(expression, library, x0);
    }

    std::cerr << minimum << std::endl;

    return 0;
}
