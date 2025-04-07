#include <string>
#include <sstream>

std::string gnuplot_convert(const std::string &expr)
{
    std::stringstream gnuplot_expr;

    auto iter = expr.begin();
    while (iter < expr.end()) {
        auto pos = expr.find('^', iter - expr.begin());
        bool found = (pos != std::string::npos);
        auto const caret = expr.begin() + (found ? pos : expr.length());
        gnuplot_expr << std::string_view(iter, caret);
        if (found)
            gnuplot_expr << "**";
        iter = caret + 1;
    }

    return gnuplot_expr.str();
}
