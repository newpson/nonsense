#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "variable.h"

static constexpr double sqrt5 = 2.23606797749978969640;

double eval(
    const std::string &expr,
    const Evaluator::library_t &lib,
    double x)
{
    dynamic_cast<Variable &>(*lib.at("x")).value = x;
    return Evaluator::eval(expr, lib);
}

double find_min(
    const std::string &expr,
    const Evaluator::library_t &lib,
    double a, double b, double eps = 1e-15)
{
    double c = (a + b) / 2;

    while (b - a > eps) {
        double f_c = eval(expr, lib, c);

        double ca = (a + c) / 2;
        double f_ca = eval(expr, lib, ca);

        double cb = (b + c) / 2;
        double f_cb = eval(expr, lib, cb);

        if (f_ca < f_c) {
            b = c;
            c = ca;
        } else if (f_cb < f_c) {
            a = c;
            c = cb;
        } else {
            a = ca;
            b = cb;
        }

        std::cerr << c << std::endl;
        std::cout << "set object circle at first "
                  << c
                  << ","
                  << eval(expr, lib, c)
                  << " radius char 0.5 fs solid fc rgb '#E41A1C' front"
                  << std::endl;
    }

    return c;
}

int main(const int argc, const char **argv)
{
    if (argc < 1 + 3)
    {
        std::cerr << "<expression> <a> <b> [eps]" << std::endl;
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

    double min_x;
    if (argc > 1 + 3)
    {
        const double eps = std::strtof(argv[4], nullptr);
        min_x = find_min(
            expression, library,
            a, b, eps);
    }
    else {
        min_x = find_min(
            expression, library,
            a, b);
    }

    std::cerr << min_x << std::endl;
    std::cout << "set object circle at first "
              << min_x
              << ","
              << eval(expression, library, min_x)
              << " radius char 0.5 fs solid fc rgb '#E41A1C' front"
              << std::endl;

    std::cout << "plot " << gnuplot_convert(expression) << " ls 2" << std::endl;

    int trash;
    std::cin >> trash;

    return 0;
}
