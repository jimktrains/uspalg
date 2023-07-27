/*
 *  uspal - micro spatial algorithms
 *  Copyright (C) 2023  James Keener <jim@jimkeener.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
//#include <bitset>
#include <cstdint>
//#include <iostream>

// I'd template this, but I'm not sure how I'd go about getting the double
// sized type for multiplication.
//
// ----------------------------------
//
// Qs10.21 is a bit of a weird size, but I need to store spherical coördinates
// so I want the integral part to be (-720, 720) (to allow for being able
// to add 2 and then normalize them and the fractional part to be able to
// handle at least 0.1-arcseconds.
//   l(2^22)/l(10)
//   6.62265990460758629472
//   l(36000)/l(10)
//   4.55630250076728726502
//   l(36000)/l(2)
//   15.13570928610439940666
//
// Single precision floating point only has a 23-bit mantissa.
//   l(2^23)/l(10)
//   6.92368990027156748993
//   (l(2)+l(720)+l(36000))/l(10)
//   7.71466499286253692048
//   (l(2)+l(720)+l(36000))/l(2)
//   25.62756238243407411754
// So, that won't work. Double-repcision would, but it's 8-bytes long,
// which seems a bit wasteful in a contrainted environment like an avr
// when there's going to be a lot of these being stored.
//
// ---------------------------
//
// I also wasn't able to understand how to use <stdfix.h> :(
// However, I'm not sure _Accum and _Fract are large enough
// for what I want anyway?
struct Qs10d21 {
  uint32_t value;

  Qs10d21() : value{0} {}

  // Converts between a the signed-magnitude and twos-complement
  // representations of the number.
  Qs10d21 twosComplementIfNegative() const {
    // This bit shifting maddness is all in an attempt
    // to avoid if-statements for negative numbers.
    //
    // It's OK to cast to int32 because we're just using this to grab
    // the sign bit and not a value.
    uint32_t signedness = (((int32_t)value) >> 31);
    uint32_t ifneg = (signedness & (~value + 1 + (((uint32_t)1) << 31)));
    uint32_t ifpos = ((~signedness) & value);
    return Qs10d21::rawInit(ifneg | ifpos);
  }

  static Qs10d21 rawInit(uint32_t x) {
    Qs10d21 t;
    t.value = x;
    return t;
  }

  explicit Qs10d21(double x) {
    int s = 0;
    if (x < 0) {
      x *= -1;
      s = ((uint32_t)1) << 31;
    }

    uint32_t i = ((uint32_t)x);
    uint32_t f = ((uint32_t)((x - i) * (((uint32_t)1) << 21)));

    value = s + (i << 21) + f;
  }

  Qs10d21 operator-() const {
    return Qs10d21::rawInit(value ^ (((uint32_t)1) << 31));
  }

  Qs10d21 operator+(const Qs10d21 &other) const {
    return Qs10d21::rawInit((twosComplementIfNegative().value +
                             other.twosComplementIfNegative().value))
        .twosComplementIfNegative();
  }

  Qs10d21 operator-(const Qs10d21 &other) const { return *this + -other; }

  Qs10d21 operator*(const Qs10d21 &other) const {
    uint32_t n =
        (value & (((uint32_t)1) << 31)) ^ (other.value & (((uint32_t)1) << 31));
    uint32_t t = value & ~(((uint32_t)1) << 31);
    uint32_t o = other.value & ~(((uint32_t)1) << 31);
    uint32_t m = (((int64_t)t) * o) >> 21;
    m |= n;
    return Qs10d21::rawInit(m);
  }

  Qs10d21 operator/(const Qs10d21 &other) const {
    uint32_t n =
        (value & (((uint32_t)1) << 31)) ^ (other.value & (((uint32_t)1) << 31));
    uint32_t t = value & ~(((uint32_t)1) << 31);
    uint32_t o = other.value & ~(((uint32_t)1) << 31);
    if (o == 0) {
      return Qs10d21::rawInit(~(((uint32_t)1) << 31) | n);
    }
    uint32_t m = ((((uint64_t)t) << 21) / o);
    return Qs10d21::rawInit(n | m);
  }

  bool operator==(const Qs10d21 &other) const { return value == other.value; }

  bool operator<(const Qs10d21 &other) const {
    return (int32_t)twosComplementIfNegative().value <
           (int32_t)other.twosComplementIfNegative().value;
  }

  bool operator<=(const Qs10d21 &other) const {
    return (int32_t)twosComplementIfNegative().value <=
           (int32_t)other.twosComplementIfNegative().value;
  }

  bool operator>(const Qs10d21 &other) const {
    return (int32_t)twosComplementIfNegative().value >
           (int32_t)other.twosComplementIfNegative().value;
  }

  bool operator>=(const Qs10d21 &other) const {
    return (int32_t)twosComplementIfNegative().value >=
           (int32_t)other.twosComplementIfNegative().value;
  }

  explicit operator double() const {
    uint32_t x = value;
    int32_t s = 1;
    if (x & (((uint32_t)1) << 31)) {
      s = -1;
    }
    // Clear the sign bit since it's been taken care of.
    x &= ~(((uint32_t)1) << 31);
    uint32_t i = x >> 21;
    double f = (x / (static_cast<double>(((uint32_t)1) << 21))) - i;

    return s * (i + f);
  }
};

static const Qs10d21 MAX_Qs10d21 = -Qs10d21::rawInit(~((uint32_t)0));
static const Qs10d21 MIN_Qs10d21 = Qs10d21::rawInit(~((uint32_t)0));

/*
int main() {
  double xx = 10.5;
  double yy = -3.25;
  Qs10d21 x(xx);
  Qs10d21 y(yy);
  std::cout << xx << std::endl;
  std::cout << ((double)x) << std::endl;
  std::cout << std::bitset<32>(x.value) << std::endl;
  std::cout << std::endl;
  std::cout << yy << std::endl;
  std::cout << ((double)y) << std::endl;
  std::cout << std::bitset<32>(y.value) << std::endl;

  std::cout << std::endl << "x-y" << std::endl;
  std::cout << (xx-yy) << std::endl;
  std::cout << (double)(x-y) << std::endl;
  std::cout << std::bitset<32>((x-y).value) << std::endl;

  std::cout << std::endl << "y-x" << std::endl;
  std::cout << (yy-xx) << std::endl;
  std::cout << (double)(y-x) << std::endl;
  std::cout << std::bitset<32>((y-x).value) << std::endl;

  std::cout << std::endl << "x+y" << std::endl;
  std::cout << (xx+yy) << std::endl;
  std::cout << (double)(x+y) << std::endl;
  std::cout << std::bitset<32>((x+y).value) << std::endl;

  std::cout << std::endl << "x*y" << std::endl;
  std::cout << (xx*yy) << std::endl;
  std::cout << (double)(x*y) << std::endl;
  std::cout << std::bitset<32>((x*y).value) << std::endl;

  std::cout << std::endl << "y*x" << std::endl;
  std::cout << (yy*xx) << std::endl;
  std::cout << (double)(y*x) << std::endl;
  std::cout << std::bitset<32>((y*x).value) << std::endl;

  std::cout << std::endl << "x/y" << std::endl;
  std::cout << (xx/yy) << std::endl;
  std::cout << (double)(x/y) << std::endl;
  std::cout << std::bitset<32>((x/y).value) << std::endl;

  std::cout << std::endl << "y/x" << std::endl;
  std::cout << (yy/xx) << std::endl;
  std::cout << (double)(y/x) << std::endl;
  std::cout << std::bitset<32>((y/x).value) << std::endl;


  return 0;
};
*/
