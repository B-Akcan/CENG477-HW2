#ifndef __COLOR_H__
#define __COLOR_H__

class Color
{
public:
    double r, g, b;

    Color();
    Color(double r, double g, double b);
    Color(const Color &other);
    Color operator+(Color );
    Color operator-(Color );
    Color operator*(double);
    Color operator/(double);
    friend std::ostream &operator<<(std::ostream &os, const Color &c);
    Color Color::round();
};

#endif