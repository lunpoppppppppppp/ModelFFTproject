#include <iostream>
#include <vector>
#include <array>

class complex {
    double real;
    double imag;
public:
    complex(double real, double imag) : real(real), imag(imag) {}
    complex() : real(0), imag(0) {}

    double get_real() const {
        return real;
    }
    double get_imag() const {
        return imag;
    }

    double norm() const {
        return real * real + imag * imag;
    }
    double abs() const {
        return std::sqrt(norm());
    }  
    complex conj() const {
        return complex(real, -imag);
    }
    complex round(std::size_t digits) {
        double coef = 1;
        for (std::size_t i = 0; i < digits; ++i) {
            coef *= 10;
        }
        this->real = std::round(this->real * coef) / coef;
        this->imag = std::round(this->imag * coef) / coef;
        return *this;
    }

    complex operator+(const complex& other) const {
        return complex(real + other.real, imag + other.imag);
    }
    complex operator-(const complex& other) const {
        return complex(real - other.real, imag - other.imag);
    }
    
    complex operator*(const double& other) const {
        return complex(real * other, imag * other);
    }
    complex operator*(const complex& other) const {
        return complex(real * other.real - imag * other.imag, real * other.imag + imag * other.real);
    }
    complex pow(std::size_t n) const {
        complex res = complex(1, 0);
        complex x = *this;
        while (n) {
            if (n & 1) {
                res = res * x;
            }
            x = x * x;
            n >>= 1;
        }
        return res;
    }
    complex operator*=(const complex& other) {
        return *this = *this * other;
    }

    complex operator/(const double& other) const {
        return complex(real / other, imag / other);
    }
    complex operator/(const complex& other) const {
        return *this * other.conj() / other.norm();
    }

    bool operator==(const complex& other) const {
        return real == other.real && imag == other.imag;
    }
    bool operator==(const double& other) const {
        return real == other && imag == 0;
    }
    bool operator!=(const complex& other) const {
        return !(*this == other);
    }
    bool operator!=(const double& other) const {
        return !(*this == other);
    }
    
    friend std::ostream& operator<< (std::ostream& stream, const complex& c) {
        return stream << "real: " << c.real << " imag: " << c.imag;
    }
};