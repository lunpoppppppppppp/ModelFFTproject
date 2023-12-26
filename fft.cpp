#include <vector>
#include <iostream>
#include <random>
#include <chrono>
#include <cassert>

#include "fft.hpp"

complex eval(std::vector<int64_t> poly, complex x) {
    complex res = complex(0, 0);
    complex x_pow = complex(1, 0);
    for (std::size_t i = 0; i < poly.size(); ++i) {
        res = res + x_pow * poly[i];
        x_pow = x_pow * x;
    }
    return res;
}

std::vector<complex> greedy_ft(std::vector<int64_t> poly, bool inverse = false) {
    std::size_t n = 1;
    while (n < poly.size()) {
        n <<= 1;
    }
    std::vector<complex> res(n);
    complex root = complex(std::cos(2 * M_PI / n), std::sin(2 * M_PI / n));
    if (inverse) root = root.conj();
    for (std::size_t i = 0; i < n; ++i) {
        res[i] = eval(poly, root.pow(i));
        if (inverse) res[i] = res[i] / n;
        res[i].round(7);
    }
    return res;
}

std::vector<int64_t> greedy_mul(std::vector<int64_t> a, std::vector<int64_t> b) {
    std::vector<int64_t> res(a.size() + b.size() - 1);
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            res[i + j] += a[i] * b[j];
        }
    }
    return res;
}

void check_complex() {
    complex val = {1e-5, 1e-10};
    std::cout << "Check round: " << val << ", " << val.round(7) << '\n';
    val = {-1e-1, -1e-2};
    std::cout << "Check round: " << val << ", " << val.round(1) << '\n';
}

void check_fft() {
    std::vector<int64_t> poly;

    int64_t val_for_dist = 100000000;
    std::size_t min_size = 1;
    std::size_t max_size = 100;
    std::size_t stress_size = 100;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(-val_for_dist, val_for_dist);

    for (std::size_t poly_size = min_size; poly_size < max_size; ++poly_size) {
        for (std::size_t times = 0; times < stress_size; ++times) {
            bool is_ok = true;
            for (std::size_t i = 0; i < poly_size; ++i) {
                poly.push_back(dist(gen));
            }
            auto poly_f = greedy_ft(poly);

            auto res_f = fft_class::fft(poly);
            if (res_f.size() != poly_f.size()) is_ok = false;
            for (std::size_t i = 0; i < res_f.size(); ++i) {
                if ((res_f[i] - poly_f[i]).norm() > 1e-5) is_ok = false;
            }
            if (!is_ok) {
                std::cout << "FFT size " << poly_size << " is not OK" << std::endl;
                std::cout << "Greedy: " << std::endl;
                for (std::size_t i = 0; i < poly_f.size(); ++i) {
                    std::cout << poly_f[i] << " ";
                }
                std::cout << std::endl;
                std::cout << "FFT: " << std::endl;
                for (std::size_t i = 0; i < res_f.size(); ++i) {
                    std::cout << res_f[i] << " ";
                }
                std::cout << std::endl;
                return;
            }

            auto res = fft_class::fft_i(poly_f);
            if (res.size() != poly.size()) is_ok = false;
            for (std::size_t i = 0; i < res.size(); ++i) {
                if (res[i] != poly[i]) is_ok = false;
            }
            if (!is_ok) {
                std::cout << "FFT_i size " << poly_size << " is not OK" << std::endl;
                std::cout << "Greedy: " << std::endl;
                for (std::size_t i = 0; i < poly.size(); ++i) {
                    std::cout << poly[i] << " ";
                }
                std::cout << std::endl;
                std::cout << "FFT_i: " << std::endl;
                for (std::size_t i = 0; i < res.size(); ++i) {
                    std::cout << res[i] << " ";
                }
                std::cout << std::endl;
                return;
            }

            poly.clear();
        }
        // if (poly_size % 20 == 0) std::cout << "FFT size " << poly_size << " is OK" << std::endl;
    }
    std::cout << "FFT transformation is OK" << std::endl;
}

void check_mul() {
    std::vector<int64_t> a;
    std::vector<int64_t> b;

    int64_t val_for_dist = 100000000;
    std::size_t min_size = 1;
    std::size_t max_size = 100;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(-val_for_dist, val_for_dist);

    for (std::size_t a_size = min_size; a_size < max_size; ++a_size) {
        for (std::size_t i = 0; i < a_size; ++i) {
            a.push_back(dist(gen));
        }
        while (a.back() == 0) a.pop_back();
        for (std::size_t b_size = min_size; b_size < max_size; ++b_size) {
            for (std::size_t i = 0; i < b_size; ++i) {
                b.push_back(dist(gen));
            }
            while (b.back() == 0) b.pop_back();

            auto greedy_res = greedy_mul(a, b);
            auto fft_res = fft_multiplication(a, b);
            
            assert(fft_res.size() == greedy_res.size());
            for (std::size_t i = 0; i < greedy_res.size(); ++i) {
                assert(fft_res[i] == greedy_res[i]);
            }
            b.clear();
        }
        a.clear();
    }
}

void check_mul_big() {
    
    std::vector<int64_t> poly_a;
    std::vector<int64_t> poly_b;
    std::vector<std::size_t> sizes = {1, 10, 100, 1000, 10000, 50000, 100000};
    std::size_t stress_size = 100;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(-100, 100);

    for (auto poly_size : sizes) {
        std::size_t greedy_cumulative_time = 0;
        std::size_t fft_cumulative_time = 0;
        for (std::size_t times = 0; times < stress_size; ++times) {
            for (std::size_t i = 0; i < poly_size; ++i) {
                poly_a.push_back(dist(gen));
            }
            while (poly_a.back() == 0) poly_a.pop_back();
            for (std::size_t i = 0; i < poly_size; ++i) {
                poly_b.push_back(dist(gen));
            }
            while (poly_b.back() == 0) poly_b.pop_back();
            auto start = std::chrono::high_resolution_clock::now();
            auto greedy_res = greedy_mul(poly_a, poly_b);
            auto end = std::chrono::high_resolution_clock::now();
            greedy_cumulative_time += static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
            start = std::chrono::high_resolution_clock::now();
            auto fft_res = fft_multiplication(poly_a, poly_b);
            end = std::chrono::high_resolution_clock::now();
            fft_cumulative_time += static_cast<std::size_t>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
             
            assert(fft_res.size() == greedy_res.size());
            for (std::size_t i = 0; i < greedy_res.size(); ++i) {
                assert(fft_res[i] == greedy_res[i]);
            }
            poly_b.clear();
            poly_a.clear();
        }
        std::cout << "Polynomial size: " << poly_size << '\n';
        std::cout << "Greedy mean time = " << greedy_cumulative_time / stress_size << " mcs" << std::endl;
        std::cout << "FFT mean time = " << fft_cumulative_time / stress_size << " mcs" << std::endl;
        std::cout << "Greedy / FFT: " << static_cast<double>(greedy_cumulative_time) / fft_cumulative_time << "\n\n";
    }
}

int main() {
    check_complex();
    check_fft();
    check_mul();
    check_mul_big();
    return 0;
}