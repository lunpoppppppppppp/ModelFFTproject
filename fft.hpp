#include <iostream>
#include "complex.hpp"
#include <cmath>
#include <cstdint>


class fft_class {
public:
    fft_class() {}

    static constexpr double eps = 1e-5;

    static std::vector<complex> fft_inner(const std::vector<complex>& poly, complex root) {
        std::size_t n = poly.size();
        if (n == 1) {
            return poly;
        }
        std::vector<complex> even(n / 2);
        std::vector<complex> odd(n / 2);
        for (std::size_t i = 0; i < n / 2; ++i) {
            even[i] = poly[2 * i];
            odd[i] = poly[2 * i + 1];
        }
        even = fft_inner(even, root.pow(2));
        odd = fft_inner(odd, root.pow(2));
        std::vector<complex> res(n);
        for (std::size_t i = 0; i < n / 2; ++i) {
            auto root_i = root.pow(i);
            res[i] = even[i] + root_i * odd[i];
            res[i + n / 2] = even[i] - root_i * odd[i];
        }
        return res;
    }

    static std::vector<complex> fft(std::vector<complex> poly) {
        std::size_t n = 1;
        while (n < poly.size()) {
            n <<= 1;
        }
        poly.resize(n);
        complex root = complex(std::cos(2 * M_PI / n), std::sin(2 * M_PI / n));

        auto res = fft_inner(poly, root);
        return res;
    }

    static std::vector<complex> fft_i(const std::vector<complex>& poly) {
        std::size_t n = 1;
        while (n < poly.size()) {
            n <<= 1;
        }
        complex root = complex(std::cos(2 * M_PI / n), -std::sin(2 * M_PI / n));

        auto res = fft_inner(poly, root);
        for (std::size_t i = 0; i < res.size(); ++i) {
            res[i] = res[i] / n;
        }
        while (res.size() && res.back().abs() < eps) {
            res.pop_back();
        }
        return res;
    }
};

std::vector<int64_t> fft_multiplication(std::vector<int64_t> a,
                                        std::vector<int64_t> b) {
    std::size_t n = a.size() + b.size() - 1;
    a.resize(n);
    b.resize(n);
    std::vector<complex> complex_a(n);
    std::vector<complex> complex_b(n);
    for (std::size_t i = 0; i < n; ++i) {
        complex_a[i] = complex(a[i], 0);
        complex_b[i] = complex(b[i], 0);
    }

    auto a_fft = fft_class::fft(complex_a);
    auto b_fft = fft_class::fft(complex_b);
    for (std::size_t i = 0; i < a_fft.size(); ++i) {
        a_fft[i] *= b_fft[i];
    }
    auto complex_res = fft_class::fft_i(a_fft);

    std::vector<int64_t> res(n);
    for (std::size_t i = 0; i < n; ++i) {
        res[i] = static_cast<int64_t>(std::round(complex_res[i].get_real()));
    }
    return res;
}
