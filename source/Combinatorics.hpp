#pragma once
#ifndef _COMBINATORICS_HPP_
#define _COMBINATORICS_HPP_

#include "RandomSampling.hpp"
#include <ostream>

typedef std::vector<bool> Mask;

std::ostream& operator<<(std::ostream& os, const Mask& mask) {
  for (bool bit : mask) {
    os << bit;
  };
  return os;
};

std::size_t nCr(std::size_t n, std::size_t r) {
  if (r > n) {
    return 0;
  };
  if (r * 2 > n) {
    r = n - r;
  };
  if (r == 0) {
    return 1;
  };
  std::size_t ncr = n;
  for (std::size_t i = 2; i <= r; ++i) {
    ncr *= n - i + 1;
    ncr /= i;
  };
  return ncr;
};

std::size_t nCrk(std::size_t n, std::size_t r, std::size_t k) {
  std::size_t ncr = nCr(n, r);
  std::size_t ncrk = ncr;
  for (std::size_t i = r + 1; i <= k; ++i) {
    ncr = (ncr * (n - i + 1)) / i;
    ncrk += ncr;
  };
  return ncrk;
};

std::vector<Mask> Combinations(
  std::size_t n, std::size_t r) {
  std::vector<Mask> masks;
  masks.reserve(nCr(n, r));
  Mask mask (n, false);
  for (std::size_t i = 0; i < r; ++i) {
    mask[i] = true;
  };
  masks.push_back(mask);
  while (std::prev_permutation(mask.begin(), mask.end())) {
    masks.push_back(mask);
  };
  return masks;
};

std::vector<Mask> Combinations(
  std::size_t n, std::size_t r, std::size_t k) {
  std::vector<Mask> masks;
  masks.reserve(nCrk(n, r, k));
  for (std::size_t i = r; i <= k; ++i) {
    std::vector<Mask> m = Combinations(n, i);
    masks.insert(masks.cend(), m.cbegin(), m.cend());
  };
  return masks;
};


struct Combinator {
  Mask mask;
  std::size_t min_set_bits, max_set_bits, set_bits;
  std::size_t revs = 0;

  Combinator(std::size_t min_set_bits, std::size_t max_set_bits) :
    min_set_bits(min_set_bits), max_set_bits(max_set_bits) {};

  Combinator(
    std::size_t n_bits,
    std::size_t min_bits,
    std::size_t max_bits) :
    mask(n_bits, false),
    min_set_bits(min_bits > n_bits ? n_bits : min_bits),
    max_set_bits(max_bits > n_bits ? n_bits : max_bits),
    set_bits(min_set_bits) {
    SetFirstBits(set_bits);
  };

  void SetFirstBits(std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {
      mask[i] = true;
    };
  };

  void UnsetBits() {
    std::fill(mask.begin(), mask.end(), false);
  };

  void operator++() {
    if (!std::prev_permutation(mask.begin(), mask.end())) {
      if (++set_bits > max_set_bits) {
        set_bits = min_set_bits;
        ++revs;
      };
      SetFirstBits(set_bits);
    };
  };
};


// Like a regular Combinator, but starting with a random mask as opposed to the
// lexicographically first mask.
struct RandomCombinator : Combinator {
  Mask init_mask;
  std::size_t init_set_bits;

  RandomCombinator(
    std::size_t n_bits,
    std::size_t min_bits,
    std::size_t max_bits,
    std::mt19937& prng) :
    Combinator(
      min_bits > n_bits ? n_bits : min_bits,
      max_bits > n_bits ? n_bits : max_bits),
    init_mask(n_bits, false) {
    std::uniform_int_distribution<std::size_t>
      distribution (min_set_bits, max_set_bits);
    init_set_bits = distribution(prng);
    set_bits = init_set_bits;
    for (std::size_t bit_idx : Sample(n_bits, init_set_bits, prng)) {
      init_mask[bit_idx] = true;
    };
    mask = init_mask;
  };

  void operator++() {
    if (!std::prev_permutation(mask.begin(), mask.end())) {
      if (++set_bits > max_set_bits) {
        set_bits = min_set_bits;
      };
      UnsetBits();
      SetFirstBits(set_bits);
    };
    if (set_bits == init_set_bits && mask == init_mask) {
      ++revs;
    };
  };
};


class Odometer {
  struct Wheel {
    std::size_t pos = 0, max = 0, revs = 0;
    Wheel(std::size_t max) : max(max) {};
    void operator++() {
      if (++pos > max) {
        pos = 0;
        ++revs;
      };
    };
  };

  std::vector<Wheel> wheels;

public:
  Odometer() = default;
  Odometer(std::size_t wheel_max, std::size_t n_wheels = 1) {
    AddWheels(wheel_max, n_wheels);
  };

  Odometer(const std::vector<std::size_t>& wheel_maxs) {
    wheels.reserve(wheel_maxs.size());
    for (std::size_t wheel_max : wheel_maxs) {
      wheels.emplace_back(wheel_max);
    };
  };

  void operator++() {
    for (Wheel& wheel : wheels) {
      std::size_t prev_revs = wheel.revs;
      ++wheel;
      if (wheel.revs == prev_revs) {
        return;
      };
    };
  };

  void AddWheels(std::size_t wheel_max, std::size_t n_wheels = 1) {
    wheels.reserve(wheels.size() + n_wheels);
    for (std::size_t i = 0; i < n_wheels; ++i) {
      wheels.emplace_back(wheel_max);
    };
  };

  std::size_t operator[](std::size_t wheel_idx) const {
    return wheels[wheel_idx].pos;
  };

  std::size_t revs() const {
    return wheels.back().revs;
  };

  std::size_t size() const {
    return wheels.size();
  };
};

#endif // !_COMBINATORICS_HPP_
