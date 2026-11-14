#pragma once

#include <random>
#include <chrono>
#include <vector>
#include <unordered_set>
#include <algorithm>


inline std::mt19937 generate()
{
    std::random_device rd {};

    // Create seed_seq with clock and 7 random numbers from std::random_device
    std::seed_seq ss{
        static_cast<std::seed_seq::result_type>(std::chrono::steady_clock::now().time_since_epoch().count()),
            rd(), rd(), rd(), rd(), rd(), rd(), rd() };

    return std::mt19937{ ss };
}

inline std::mt19937& rng()
{
    thread_local std::mt19937 gen {generate()};
    return gen;
}

inline double random_uniform(double lower=0.0, double upper=1.0)
{
    std::uniform_real_distribution<double> dist (lower, upper);
    return dist(rng());
}

template <typename T>
T rand_int(T lower, T upper)
{
    std::uniform_int_distribution<T> dist(lower, upper);
    return dist(rng());
}

template<typename T>
std::vector<T> rand_int(T lower, T upper, size_t k)
{
    std::unordered_set<T> result;
    std::uniform_int_distribution<T> dist(lower, upper);
    while (result.size() < k) { result.insert(dist(rng())); }
    return std::vector<T>(result.begin(), result.end());
}

template <typename T>
T get_sample(const std::vector<T>& vec) {
    std::uniform_int_distribution<size_t> dist(0, vec.size() - 1);
    return vec[dist(rng())];
}


template <typename T>
T pop_random(std::vector<T>& vec)
{
    std::uniform_int_distribution<size_t> dist(0, vec.size() - 1);
    size_t index = dist(rng());

    T value = vec[index];
    std::swap(vec[index], vec.back());
    vec.pop_back();

    return value;
}

template <typename T>
std::size_t weighted_sample(const std::vector<T>& weights)
{
    std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
    return dist(rng());
}

template <typename T, typename U>
T weighted_sample(const std::vector<T>& items, const std::vector<U>& weights)
{
    std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
    return items[dist(rng())];
}

template <typename T>
std::vector<T> sample_without_replacement(T values, size_t k) {
    std::vector<T> pool;
    pool.reserve(values);

    for (T i = 0; i < values; ++i)
        pool.push_back(i);

    std::shuffle(pool.begin(), pool.end(), rng());
    pool.resize(k);
    return pool;
}

template <typename T>
std::vector<T> sample_without_replacement(const std::vector<T>& vec, size_t n) {
    if (n > vec.size()) {
        throw std::invalid_argument("Sample size cannot be larger than vector size");
    }

    std::vector<T> result;
    result.reserve(n);

    // std::sample will pick n unique elements without replacement
    std::sample(vec.begin(), vec.end(), std::back_inserter(result), n, rng());

    return result;
}

template <typename T>
std::vector<T> weighted_sample_without_replacement(const std::vector<T>& items, std::vector<double> weights, std::size_t k){
    if (items.size() != weights.size())
        throw std::invalid_argument("items/weights size mismatch");
    if (k > items.size())
        throw std::invalid_argument("k > items.size()");
    
    std::vector<T> result;
    result.reserve(k);

    for (std::size_t i {}; i < k; ++i) {
        auto ind {weighted_sample(weights)};
        weights[ind] = 0;
        result.push_back(items[ind]);
    }

    return result;
}