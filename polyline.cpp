#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <vector>

#include "evaluator.h"
#include "gnuplot-helpers.h"
#include "math.h"
#include "variable.h"

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

std::vector<Vector2D> read_polygon()
{
    std::vector<Vector2D> points;
    std::string line;
    double x = 0.0;
    double y = 0.0;

    while (std::getline(std::cin, line)) {
        if (line == "e") {
            break;
        } else {
            std::stringstream stream(line);
            stream >> x >> y;
            points.push_back(Vector2D(x, y));
        }
    }
    return points;
}

int winding_number(const Vector2D &target, const std::vector<Vector2D> &frame)
{
    // modified David G. Alciatore conception, split halfs and wholes
    int w = 0;
    int hw = 0;
    for (int i = 0; i < frame.size() - 1; ++i) {
        if (frame[i] == target)
            return 1; // TODO (may be not TODO) modify final sum with this value somehow

        if (frame[i].x >= target.x) {
            if (frame[i].y < target.y) {
                if (frame[i+1].y > target.y)
                    ++w;
                else if (frame[i+1].y == target.y)
                    ++hw;
            } else if (frame[i].y > target.y) {
                if (frame[i+1].y < target.y)
                    --w;
                else if (frame[i+1].y == target.y)
                    --hw;
            } else if (frame[i].y == target.y) {
                if (frame[i+1].y > target.y)
                    ++hw;
                else if (frame[i+1].y < target.y)
                    --hw;
            }
        }
    }

    return w + hw/2 + hw%2;
}

void plot_line(const Vector2D &begin, const Vector2D &end, const Vector2D &h)
{
    enum Cardinal {
        NORTH_EAST, // 0b00
        NORTH_WEST, // 0b01
        SOUTH_EAST, // 0b10
        SOUTH_WEST, // 0b11
        NUM_CARDINALS,
    };

    // TODO integer Vector2D type (may be templates)

    // cell pivot shift
    const Vector2D shift[NUM_CARDINALS] = {
        [NORTH_EAST] = {1.0, 1.0},
        [NORTH_WEST] = {0.0, 1.0},
        [SOUTH_EAST] = {1.0, 0.0},
        [SOUTH_WEST] = {0.0, 0.0},
    };

    // move vector when det > 0
    const Vector2D move_p[NUM_CARDINALS] = {
        [NORTH_EAST] = { 0.0,  1.0},
        [NORTH_WEST] = {-1.0,  0.0},
        [SOUTH_EAST] = { 1.0,  0.0},
        [SOUTH_WEST] = { 0.0, -1.0},
    };

    // move vector when det < 0
    const Vector2D move_n[NUM_CARDINALS] = {
        [NORTH_EAST] = { 1.0,  0.0},
        [NORTH_WEST] = { 0.0,  1.0},
        [SOUTH_EAST] = { 0.0, -1.0},
        [SOUTH_WEST] = {-1.0,  0.0},
    };

    const Vector2D guide = end - begin;
    // TODO maybe more effective way to set direction using bit arithmetics?
    const Cardinal direction = guide.y > 0 ? guide.x > 0 ? NORTH_EAST : NORTH_WEST
                                           : guide.x > 0 ? SOUTH_EAST : SOUTH_WEST;

    const Vector2D end_floored = Vector2D(std::floor(end.x/h.x),
                                          std::floor(end.y/h.y));
    Vector2D cur(std::floor(begin.x/h.x),
                 std::floor(begin.y/h.y));
    std::cerr << cur*h << std::endl;
    // std::cout << "set object rect from "
    //           << cur.x*h.x << "," << cur.y*h.y
    //           << " to "
    //           << (cur.x + 1.0)*h.x << "," << (cur.y + 1.0)*h.y
    //           << std::endl;
    int i = 0;
    Vector2D aux = h * (cur + shift[direction]) - begin;
    while (std::abs(aux.x) <= std::abs(guide.x) && std::abs(aux.y) <= std::abs(guide.y) && i < 100) { // FIXME make integer comparison
        aux = h * (cur + shift[direction]) - begin;
        const double det = aux.x * guide.y - guide.x * aux.y;
        cur = cur + (det > 0 ? move_p[direction]
                   : det < 0 ? move_n[direction]
                             : move_p[direction] + move_n[direction]);
        std::cerr << cur*h << std::endl;
        // std::cout << "set object rect from "
        //           << cur.x*h.x << "," << cur.y*h.y
        //           << " to "
        //           << (cur.x + 1.0)*h.x << "," << (cur.y + 1.0)*h.y
        //           << std::endl;
        ++i;
    }
}

int main(int argc, char **argv)
{
    // TODO
    plot_line({0.1, 1.0}, {0.2, 0.9}, {0.1, 0.1});
    return 0;

    std::shared_ptr<Variable> t(new Variable());
    Evaluator::library_t library = {
        {"t", t},
    };
    library.insert(Evaluator::default_library.begin(),
                   Evaluator::default_library.end());

    Vector2D step_space(0.1, 0.1);
    std::vector<Vector2D> path;
    std::vector<Vector2D> frame;

    std::string expr_path_x;
    std::string expr_path_y;
    double path_ta = 0.0;
    double path_tb = 1.0;
    double path_tau = 0.01;
    std::string expr_frame_x;
    std::string expr_frame_y;
    double frame_ta = 0.0;
    double frame_tb = 1.0;
    double frame_tau = 0.01;
    for (int i = 1; i < argc; ++i) {
        const std::string par(argv[i]);
        if (par == "--path") {
            expr_path_x = argv[i+1];
            expr_path_y = argv[i+2];
            path_ta = std::strtod(argv[i+3], nullptr);
            path_tb = std::strtod(argv[i+4], nullptr);
            path_tau = std::strtod(argv[i+5], nullptr);
        }
        else if (par == "--frame") {
            expr_frame_x = argv[i+1];
            expr_frame_y = argv[i+2];
            frame_ta = std::strtod(argv[i+3], nullptr);
            frame_tb = std::strtod(argv[i+4], nullptr);
            frame_tau = std::strtod(argv[i+5], nullptr);
        }
        else if (par == "--cell") {
            step_space.x = std::strtod(argv[i+1], nullptr);
            step_space.y = std::strtod(argv[i+2], nullptr);
        }
    }

    path = expr_path_x.empty() ? read_polygon() : subdivide(expr_path_x, expr_path_y, library, false, path_ta, path_tb, path_tau);
    frame = expr_frame_x.empty() ? read_polygon() : subdivide(expr_frame_x, expr_frame_y, library, false, frame_ta, frame_tb, frame_tau);

    std::cout << "unset object" << std::endl;
    for (int i = 0; i < path.size() - 1; ++i) {
        const bool is_inside = (winding_number(path[i], frame) != 0 &&
                                winding_number(path[i+1], frame) != 0);
        if (is_inside)
            plot_line(path[i], path[i+1], step_space);
        else {
            std::cerr << "Out of bounds at " << path[i] << " -> " << path[i+1] << std::endl;
            break;
        }
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

    return 0;
}
