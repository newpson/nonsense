/**
 * Численное решение двухточечной краевой задачи для ЛОДУ второго порядка
 * на отрезке [a, b] с шагом (b - a)/N.
 */

#include "evaluator.h"
#include "variable.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

double eval(const std::string &expression,
            const Evaluator::library_t &library,
            const double x)
{
    dynamic_cast<Variable &>(*library.at("x")).value = x;
    return Evaluator::eval(expression, library);
}

int main(int argc, char **argv)
{
    if (argc < 1 + 7) {
        std::cerr << "<q(x)> <f(x)> <a> <b> <u_a> <u_b> <N>" << std::endl;
        return 1;
    }

    std::shared_ptr<Variable> x(new Variable());
    Evaluator::library_t library = {{"x", x}};
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    const std::string q(argv[1]);
    const std::string f(argv[2]);
    const double a = std::strtod(argv[3], nullptr);
    const double b = std::strtod(argv[4], nullptr);
    const double u_a = std::strtod(argv[5], nullptr);
    const double u_b = std::strtod(argv[6], nullptr);
    const int N = std::strtol(argv[7], nullptr, 10);

    const double h = (b - a) / (N - 1);

    std::vector<double> A(N);
    std::vector<double> B(N);
    std::vector<double> C(N);
    std::vector<double> F(N);
    for (int i = 0; i < N; ++i) {
        A[i] = (i > 0 && i < N-1) ? -1.0 : 0.0;
        B[i] = (i > 0 && i < N-1) ? -1.0 : 0.0;
        C[i] = (i > 0 && i < N-1) ? h*h*eval(q, library, a + i*h) + 2.0 : 1.0;
        F[i] = (i == 0) ? u_a : (i == N-1) ? u_b : h*h*eval(f, library, a + i*h);
    }

    std::vector<double> mu(N);
    std::vector<double> nu(N);
    // there are mu[-1] = nu[i-1] = 0.0;
    mu[0] = F[0] / C[0];
    nu[0] = -B[0] / C[0];

    for (int i = 1; i < mu.size(); ++i) {
        mu[i] = (F[i] - mu[i-1]*A[i]) / (C[i] + nu[i-1]*A[i]);
        nu[i] = -B[i] / (C[i] + nu[i-1]*A[i]);
    }
    for (int i = mu.size() - 2; i >= 0; --i)
        mu[i] = mu[i] + nu[i]*mu[i+1];

    for (int i = 0; i < N; ++i)
        std::cout << a + i*h << " " << mu[i] << std::endl;
    return 0;
}
