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

void solve(const Evaluator::library_t &lib,
           const double xa, const double xb, const int num_steps,
           const double ta, const double tb, const int num_layers,
           const std::string &f, const std::string &phi, const std::string &xi,
           const double c,
           const std::string *u = nullptr)
{
    const double h = (xb - xa) / (num_steps-1);
    const double tau = (tb - ta) / (num_layers-1);

    if (c < 0.0 && h/tau > -c)
        std::cerr << "warning: unstable" << std::endl;

    std::vector<double> v1(num_steps);
    std::vector<double> v2(num_steps);

    std::vector<double> *y_cur = &v1;
    std::vector<double> *y_prev = &v2;

    if (u != nullptr) {
        double err_layer_max = 0.0;
        double err_grid_max = err_layer_max;

        double x = xa;
        for (double &value: *y_prev) {
            value = eval(phi, lib, x, ta);
            const double exact_u = eval(*u, lib, x, ta);
            double err = std::abs(exact_u - value);
            if (err > err_layer_max)
                err_layer_max = err;
            std::cout << "j = " << 0 << "; "
                      << "x = " << x << "; "
                      << "y = " << value << "; "
                      << "u = " << exact_u << "; "
                      << "err = "  << err
                      << std::endl;
            // std::cout << x << " " << ta << " " << value << std::endl;
            x += h;
        }

        double t = ta + tau;
        for (int j = 1; j < num_layers; ++j) {
            x = xa;
            (*y_cur)[0] = eval(xi, lib, x, t);
            const double exact_u = eval(*u, lib, x, t);
            double err = std::abs(exact_u - (*y_cur)[0]);
            if (err > err_layer_max)
                err_layer_max = err;
            std::cout << "j = " << j << "; "
                      << "x = " << x << "; "
                      << "y = " << (*y_cur)[0] << "; "
                      << "u = " << exact_u << "; "
                      << "err = " << err
                      << std::endl;
            // std::cout << x << " " << t << " " << (*y_cur)[0] << std::endl;
            for (int n = 1; n < num_steps; ++n) {
                x += h;
                (*y_cur)[n] = (eval(f, lib, x, t) + c/h * (*y_cur)[n-1] + (*y_prev)[n] / tau)
                              / (c/h + 1/tau);
                const double exact_u = eval(*u, lib, x, t);
                err = std::abs(exact_u - (*y_cur)[n]);
                if (err > err_layer_max)
                    err_layer_max = err;
                std::cout << "j = " << j << "; "
                          << "x = " << x << "; "
                          << "y = " << (*y_cur)[n] << "; "
                          << "u = " << exact_u << "; "
                          << "err = " << err
                          << std::endl;
                // std::cout << x << " " << t << " " << (*y_cur)[n] << std::endl;
            }
            t += tau;
            std::swap(y_prev, y_cur);

            if (err_layer_max > err_grid_max)
                err_grid_max = err_layer_max;
            std::cout << "info: err_layer_max = " << err_layer_max << std::endl;
        }
        std::cout << "info: err_grid_max = " << err_grid_max << std::endl;
    } else {
        double x = xa;
        for (double &value: *y_prev) {
            value = eval(phi, lib, x, ta);
            std::cout << x << " " << ta << " " << value << std::endl;
            x += h;
        }

        double t = ta + tau;
        for (int j = 1; j < num_layers; ++j) {
            x = xa;
            (*y_cur)[0] = eval(xi, lib, x, t);
            std::cout << x << " " << t << " " << (*y_cur)[0] << std::endl;
            for (int n = 1; n < num_steps; ++n) {
                x += h;
                (*y_cur)[n] = (eval(f, lib, x, t) + c/h * (*y_cur)[n-1] + (*y_prev)[n] / tau)
                              / (c/h + 1/tau);
                std::cout << x << " " << t << " " << (*y_cur)[n] << std::endl;
            }
            t += tau;
            std::swap(y_prev, y_cur);
        }
    }
}

int main(int argc, char **argv)
{
    if (argc < 1 + 10) {
        std::cerr << "<xa> <xb> <num_steps> <ta> <tb> <num_layers> <f(x,t)> <phi(x)> <xi(t)> <c> [u(x, t)]" << std::endl;
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
    const double num_steps = std::strtol(argv[3], nullptr, 10);

    const double ta = std::strtod(argv[4], nullptr);
    const double tb = std::strtod(argv[5], nullptr);
    const double num_layers = std::strtol(argv[6], nullptr, 10);

    const std::string expr_f(argv[7]);
    const std::string expr_phi(argv[8]);
    const std::string expr_xi(argv[9]);

    const double c = std::strtod(argv[10], nullptr);

    if (argc == 1 + 10 + 1) {
        std::cerr << "info: exact solution was specified" << std::endl;
        const std::string expr_u(argv[11]);
        solve(library,
              xa, xb, num_steps,
              ta, tb, num_layers,
              expr_f, expr_phi, expr_xi,
              c, &expr_u);
    } else {
        std::cerr << "info: blind solving" << std::endl;
        solve(library,
              xa, xb, num_steps,
              ta, tb, num_layers,
              expr_f, expr_phi, expr_xi,
              c);
    }

    return 0;
}
