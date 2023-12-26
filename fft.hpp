#include <iostream>
#include "complex.hpp"


class fft_class {
public:
    fft_class() {}

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

    static std::vector<complex> fft(const std::vector<int64_t>& poly) {
        std::size_t n = 1;
        while (n < poly.size()) {
            n <<= 1;
        }
        std::vector<complex> prepared_poly(n);
        for (std::size_t i = 0; i < poly.size(); ++i) {
            prepared_poly[i] = complex(poly[i], 0);
        }
        complex root = complex(std::cos(2 * M_PI / n), std::sin(2 * M_PI / n));

        auto res = fft_inner(prepared_poly, root);
        return res;
    }

    static std::vector<int64_t> fft_i(const std::vector<complex>& poly) {
        std::size_t n = 1;
        while (n < poly.size()) {
            n <<= 1;
        }
        complex root = complex(std::cos(2 * M_PI / n), -std::sin(2 * M_PI / n));

        auto complex_res = fft_inner(poly, root);
        std::vector<int64_t> res(n);
        for (std::size_t i = 0; i < n; ++i) {
            res[i] = static_cast<int64_t>(std::round(complex_res[i].get_real() / n));
        }
        while (res.back() == 0) {
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
    auto a_fft = fft_class::fft(a);
    auto b_fft = fft_class::fft(b);
    assert(a_fft.size() == b_fft.size());
    for (std::size_t i = 0; i < a_fft.size(); ++i) {
        a_fft[i] *= b_fft[i];
    }
    auto res = fft_class::fft_i(a_fft);
    return res;
}