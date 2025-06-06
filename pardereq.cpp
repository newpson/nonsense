#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "math.h"
#include "variable.h"

#include "Eigen/Dense"

double eval(const std::string &expression,
            const Evaluator::library_t &library,
            const double x, const double t)
{
    dynamic_cast<Variable &>(*library.at("x")).value = x;
    dynamic_cast<Variable &>(*library.at("t")).value = t;
    return Evaluator::eval(expression, library);
}

void vector_print(const std::vector<double> &v)
{
    for (const double value: v)
        std::cout << value << " ";
    std::cout << std::endl;
}

int main(int argc, char **argv)
{
    if (argc < 1 + 8) {
        std::cerr << "<xa> <xb> <ta> <tb> <f(x,t)> <phi(x)> <xi(t)> <c>" << std::endl;
        return 0;
    }
    Evaluator::library_t library = {
        {"x", std::shared_ptr<Evaluable>(new Variable())},
        {"t", std::shared_ptr<Evaluable>(new Variable())},
    };
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    const double xa = std::strtod(argv[1], nullptr);
    const double xb = std::strtod(argv[2], nullptr);
    const double ta = std::strtod(argv[3], nullptr);
    const double tb = std::strtod(argv[4], nullptr);
    const std::string expr_f(argv[5]);
    const std::string expr_phi(argv[6]);
    const std::string expr_xi(argv[7]);
    const double c = std::strtod(argv[8], nullptr);

    const double num_steps = 10;
    const double num_layers = 10;
    const double h = (xb - xa) / (num_steps-1);
    const double tau = (tb - ta) / (num_layers-1);

    if (c < 0.0 && h/tau > -c)
        std::cerr << "warning: unstable" << std::endl;

    std::vector<double> v1(num_steps);
    std::vector<double> v2(num_steps);

    std::vector<double> *u_cur = &v1;
    std::vector<double> *u_prev = &v2;

    double x = xa;
    for (double &value: *u_prev) {
        value = eval(expr_phi, library, x, ta);
        std::cout << x << " " << ta << " " << value << std::endl;
        x += h;
    }

    double t = ta;
    for (int j = 0; j < num_layers; ++j) {
        (*u_cur)[0] = eval(expr_xi, library, xa, t);
        x = xa;
        std::cout << x << " " << t << " " << (*u_cur)[0] << std::endl;
        for (int n = 1; n < num_steps; ++n) {
            x += h;
            (*u_cur)[n] = (eval(expr_f, library, x, t) + c/h * (*u_cur)[n-1] + 1/tau * (*u_prev)[n]) / (c/h + 1/tau);
            std::cout << x << " " << t << " " << (*u_cur)[n] << std::endl;
        }
        t += tau;
        std::swap(u_prev, u_cur);
    }

    return 0;
}
