#ifndef _RANDOM_SAMPLING_HPP_
#define _RANDOM_SAMPLING_HPP_

#include <boost/dynamic_bitset.hpp>
#include <random>
#include <cassert>
#include <concepts>
#include <algorithm>

template <class T>
requires std::strict_weak_order<std::less<T>, T, T>
std::vector<std::size_t> SortPermutation(const std::vector<T>& v) {
  std::vector<std::size_t> p (v.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
    [&](std::size_t i, std::size_t j) {
      return v[i] < v[j];
    });
  return p;
};

template <class T, class Compare>
requires std::strict_weak_order<Compare, T, T>
std::vector<std::size_t> SortPermutation(
  const std::vector<T>& v,
  const Compare& compare) {
  std::vector<std::size_t> p (v.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
    [&](std::size_t i, std::size_t j) {
      return compare(v[i], v[j]);
    });
  return p;
};

template <class T>
requires std::copy_constructible<T>
std::vector<T> ApplyPermutation(
  const std::vector<T>& v,
  const std::vector<std::size_t>& p) {
  std::vector<T> sorted_v;
  sorted_v.reserve(v.size());
  std::transform(p.begin(), p.end(), std::back_inserter(sorted_v),
    [&](std::size_t i) {
      return v[i];
    });
  return sorted_v;
};

template <class T>
requires std::swappable<T>
void ApplyPermutationInPlace(
  std::vector<T>& v,
  const std::vector<std::size_t>& p) {
  std::vector<bool> done (v.size());
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (done[i]) {
      continue;
    };
    done[i] = true;
    std::size_t prev_j = i;
    std::size_t j = p[i];
    while (i != j) {
      std::swap(v[prev_j], v[j]);
      done[j] = true;
      prev_j = j;
      j = p[j];
    };
  };
};

// Uniformly sample an integer from the [0, n) range.
template <class G>
requires std::uniform_random_bit_generator<G>
std::size_t Sample(std::size_t n, G& prng) {
  assert(n > 0);
  if (n < 2) {
    return 0;
  };
  std::uniform_int_distribution<std::size_t> distribution (0, n - 1);
  return distribution(prng);
};

// Uniformly sample a set bit from a bitset.
template <class G>
requires std::uniform_random_bit_generator<G>
std::size_t Sample(const boost::dynamic_bitset<>& bits, G& prng) {
  std::size_t bit_position = Sample(bits.count(), prng);
  std::size_t bit_idx = bits.find_first();
  for (std::size_t i = 0; i < bit_position; ++i) {
    bit_idx = bits.find_next(bit_idx);
  };
  return bit_idx;
};

// Uniformly sample without replacement m integers from the [0, n) range.
template <class G>
requires std::uniform_random_bit_generator<G>
std::vector<std::size_t> Sample(std::size_t n, std::size_t m, G& prng) {
  assert(m <= n);
  std::vector<std::size_t> indices (n);
  std::vector<std::size_t> sample (m);
  std::iota(indices.begin(), indices.end(), 0);
  std::sample(indices.begin(), indices.end(), sample.begin(), m, prng);
  return sample;
};

// Perform weighted sampling of a single element.
template <class G>
requires std::uniform_random_bit_generator<G>
std::size_t WeightedSample(const std::vector<double>& w, G& prng) {
  return std::discrete_distribution<std::size_t>(w.cbegin(), w.cend())(prng);
};

// Given a vector with n strictly positive weights, and a sample size of m
// (where m <= n), take m weighted samples without replacement.
// Efraimidis, P. S. & Spirakis, P. G. Weighted random sampling with a reservoir.
// Information Processing Letters 97, 181-185 (2006)
template <class G>
requires std::uniform_random_bit_generator<G>
std::vector<std::size_t> WeightedSample(
  const std::vector<double>& w, std::size_t m, G& prng) {
	// Calculate the key values of each weight.
	std::size_t n = w.size();
	assert(m <= n);
  std::vector<double> keys (n);
  static std::uniform_real_distribution<double> zero_to_one (0.0, 1.0);
  std::transform(w.begin(), w.end(), keys.begin(),
    [&](double x) {
      return std::pow(zero_to_one(prng), 1.0 / x);
    });
  // Determine the keys sorted descending order.
  static std::greater<double> comparator;
  std::vector<std::size_t> p = SortPermutation(keys, comparator);
  // Return the m indices with the highest key values.
  std::vector<std::size_t>::const_iterator it = p.begin();
  std::vector<std::size_t> sample (it, it + m);
	return sample;
};

template <class T, class G>
std::vector<T> WeightedShuffle(
  const std::vector<T>& v, const std::vector<double>& w, G& prng) {
  assert(v.size() == w.size());
  std::vector<std::size_t> p = WeightedSample(w, w.size(), prng);
  return ApplyPermutation(v, p);
};

template <class T, class G>
void WeightedShuffleInPlace(
  std::vector<T>& v, const std::vector<double>& w, G& prng) {
  assert(v.size() == w.size());
  std::vector<std::size_t> p = WeightedSample(w, w.size(), prng);
  ApplyPermutationInPlace(v, p);
};

#endif // !_RANDOM_SAMPLING_HPP_
