#include <bitset>
#include <cstdint>
#include <iostream>
// I'd template this, but I'm not sure how I'd go about getting the double
// sized type for multiplication.
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
// I also wasn't able to understand how to use <stdfix.h> :(
class Qs10d21 {

  Qs10d21(uint32_t v) : value{v} {};

  // Converts between a the signed-magnitude and twos-complement
  // representations of the number.
  Qs10d21 twosComplementIfNegative() const {
    // This bit shifting maddness is all in an attempt
    // to avoid if-statements for negative numbers.
    //
    // It's OK to cast to int32 because we're just using this to grab
    // the sign bit and not a value.
    return Qs10d21{((((int32_t)value) >> 31) & (~value + 1 + (1 << 31))) |
                   (((~((int32_t)value) >> 31) & value))};
  };

public:
  uint32_t value;

  static Qs10d21 rawInit(uint32_t x) { return Qs10d21(x); }

  Qs10d21(double x) {
    int s = 0;
    if (x < 0) {
      x *= -1;
      s = 1 << 31;
    }

    uint32_t i = ((uint32_t)x);
    uint32_t f = ((uint32_t)((x - i) * (1 << 21)));

    value = s + (i << 21) + f;
  };

  Qs10d21 operator-() const { return Qs10d21{value ^ (1 << 31)}; };

  Qs10d21 operator+(const Qs10d21 &other) const {
    return Qs10d21{(twosComplementIfNegative().value +
                    other.twosComplementIfNegative().value)}
        .twosComplementIfNegative();
  }

  Qs10d21 operator-(const Qs10d21 &other) const { return *this + -other; }

  Qs10d21 operator*(const Qs10d21 &other) const {
    uint32_t n = (value & (1 << 31)) ^ (other.value & (1 << 31));
    uint32_t t = value & ~(1 << 31);
    uint32_t o = other.value & ~(1 << 31);
    uint32_t m = (((int64_t)t) * o) >> 21;
    m |= n;
    return Qs10d21(m);
  }

  Qs10d21 operator/(const Qs10d21 &other) const {
    uint32_t n = (value & (1 << 31)) ^ (other.value & (1 << 31));
    uint32_t t = value & ~(1 << 31);
    uint32_t o = other.value & ~(1 << 31);
    uint32_t m = ((((uint64_t)t) << 21) / o);
    return Qs10d21(n | m);
  }

  bool operator==(const Qs10d21 &other) const { return value == other.value; };

  bool operator<(const Qs10d21 &other) const {
    return (int32_t)twosComplementIfNegative().value <
           (int32_t)other.twosComplementIfNegative().value;
  };

  bool operator<=(const Qs10d21 &other) const {
    return (int32_t)twosComplementIfNegative().value <=
           (int32_t)other.twosComplementIfNegative().value;
  };

  bool operator>(const Qs10d21 &other) const {
    return (int32_t)twosComplementIfNegative().value >
           (int32_t)other.twosComplementIfNegative().value;
  };

  bool operator>=(const Qs10d21 &other) const {
    return (int32_t)twosComplementIfNegative().value >=
           (int32_t)other.twosComplementIfNegative().value;
  };

  explicit operator double() const {
    uint32_t x = value;
    int32_t s = 1;
    if (x & (1 << 31)) {
      s = -1;
    }
    // Clear the sign bit since it's been taken care of.
    x &= ~(1 << 31);
    uint32_t i = x >> 21;
    double f = (x / ((double)(1 << 21))) - i;

    return s * (i + f);
  };
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
