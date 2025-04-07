#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "variable.h"

#include <raylib.h>

double eval(
    const std::string &expr,
    const Evaluator::library_t &lib,
    const double t)
{
    dynamic_cast<Variable &>(*lib.at("t")).value = t;
    return Evaluator::eval(expr, lib);
}

std::vector<double> subdivide(const std::string &expr,
                              const Evaluator::library_t &lib,
                              const double a, const double b)
{
    const double h = 0.1;
    std::vector<double> points;
    for (double t = a; t <= b; t += h) {
        points.push_back(eval(expr, lib, t));
    }
    return points;
}

int winding_number(const std::vector<double> &values_x,
                    const std::vector<double> &values_y,
                    const double x0, const double y0)
{
    int w = 0;
    for (int i = 0; i < values_x.size() - 1; ++i) {
        if (values_x[i] >= x0) {
            if (values_y[i] <= y0 && y0 < values_y[i+1])
                ++w;
            else if (values_y[i] > y0 && y0 >= values_y[i+1])
                --w;
        }
    }
    return w;
}

void plot_line(double x1, double y1, double x2, double y2)
{
    const double hx = 0.5; // шаг по x
    const double hy = 0.5; // шаг по y
    const double vx = x2 - x1; // вектор приближаемой прямой (x)
    const double vy = y2 - y1; // вектор приближаемой прямой (y)
    const double dirx = copysign(1.0, vx); // направление вектора прямой (x)
    const double diry = copysign(1.0, vy); // направление вектора прямой (y)
    const double invertor = copysign(1.0, vx*vy); // инверсия осей в чётных квадрантах

    const double startx = dirx > 0 ? floor(x1/hx)*hx : ceil(x1/hx)*hx; // начальная ячейка (x)
    const double starty = diry > 0 ? floor(y1/hy)*hy : ceil(y1/hy)*hy; // начальная ячейка (y)
    const double endx = dirx > 0 ? floor(x2/hx)*hx : ceil(x2/hx)*hx; // конечная ячейка (x)
    const double endy = diry > 0 ? floor(y2/hy)*hy : ceil(y2/hy)*hy; // конечная ячейка (y)

    double x = startx; // вектор приближения (x)
    double y = starty; // вектор приближения (y)

    int i = 0; // счётчик итераций (для отладки)
    while ((fabs(x - endx) > hx/1000.0 || fabs(y - endy) > hy/1000.0) && i < 1000) {
        std::cout << "set object rect from " << x << "," << y << " to " << x + dirx*hx << "," << y + diry*hy << std::endl;
        x += hx*dirx;
        y += hy*diry;
        double D = invertor*((y - y1)*vx - (x - x1)*vy);
        if (D > 0)
            y -= hy*diry;
        else if (D < 0)
            x -= hx*dirx;
        ++i;
    }
    std::cout << "set object rect from " << x << "," << y << " to " << x + dirx*hx << "," << y + diry*hy << std::endl;

    if (i >= 1000) {
        std::cerr << "*EXTREME* Overflow at (" << x1 << ", " << y1 << ")->(" << x2 << ", " << y2 << ")" << std::endl;
        std::cerr << "Last should be (" << endx << ", " << endy << ")" << std::endl;
    }
}



int main(int argc, char **argv)
{
    std::shared_ptr<Variable> t(new Variable());
    Evaluator::library_t library = {
        {"t", t},
    };
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    const double ta = std::strtod(argv[3], nullptr);
    const double tb = std::strtod(argv[4], nullptr);
    auto values_x = subdivide("t*sin(t)", library, 0.0, 5*3.14);
    auto values_y = subdivide("t*cos(t)", library, 0.0, 5*3.14);
    auto frame_x = subdivide("12*cos(t)", library, 0.0, 2*3.15);
    auto frame_y = subdivide("4*sin(t)", library, 0.0, 2*3.15);

    for (int i = 0; i < values_x.size() - 1; ++i) {
        const bool is_inside = (winding_number(frame_x, frame_y, values_x[i], values_y[i]) != 0 &&
                                winding_number(frame_x, frame_y, values_x[i+1], values_y[i+1]) != 0);
        if (is_inside)
            plot_line(values_x[i], values_y[i], values_x[i+1], values_y[i+1]);
    }
    std::cout << "$points << e" << std::endl;
    for (int i = 0; i < values_x.size(); ++i) {
        std::cout << values_x[i] << " " << values_y[i] << std::endl;
    }
    std::cout << "e" << std::endl;
    std::cout << "$frame << e" << std::endl;
    for (int i = 0; i < frame_x.size(); ++i) {
        std::cout << frame_x[i] << " " << frame_y[i] << std::endl;
    }
    std::cout << "e" << std::endl;
    std::cout << "plot $points with lines ls 1, $frame with lines ls 2" << std::endl;

    int trash;
    std::cin >> trash;




    // ============
    //    RAYLIB
    // ============

    // const char * const &expr_x = "(t^3-3*t)*50+200";
    // const char * const &expr_y = "(t^2-3)*50+100";
    // auto values_x = subdivide(expr_x, library, -std::sqrt(3), std::sqrt(3));
    // auto values_y = subdivide(expr_y, library, -std::sqrt(3), std::sqrt(3));

    // InitWindow(400, 200, "Demo");
    // SetTargetFPS(30);

    // struct Vector2 from;
    // struct Vector2 to;
    // while (!WindowShouldClose())
    // {
    //     int winding_number = 0;
    //     double x0 = GetMouseX();
    //     double y0 = GetMouseY();
    //     BeginDrawing();
    //     ClearBackground(RAYWHITE);
    //     for (int i = 0; i < values_x.size() - 1; ++i) {
    //         if (values_x[i] >= x0) {
    //             if (values_y[i] <= y0 && y0 < values_y[i+1])
    //                 ++winding_number;
    //             else if (values_y[i] > y0 && y0 >= values_y[i+1])
    //                 --winding_number;
    //         }
    //         from.x = values_x[i];
    //         from.y = values_y[i];
    //         to.x = values_x[i+1];
    //         to.y = values_y[i+1];
    //         DrawLineEx(from, to, 2.0, RED);
    //     }
    //     DrawCircle(x0, y0, 2, BLUE);
    //     if (winding_number != 0)
    //         DrawText("inside!", 20, 20, 20, DARKGRAY);
    //     DrawText(TextFormat("%d", winding_number), 2, 2, 20, DARKGRAY);

    //     EndDrawing();
    // }

    // CloseWindow();

    return 0;
}

