#include <iostream>
#include <string>

struct Time {
  long days;
  int hour;
  int minutes;
  int seconds;

  Time(long d, int h, int m, int s)
    : days{d}
    , hour{h}
    , minutes{m}
    , seconds{s} {

    // Leap Seconds aren't handled -_- and I don't care :o -- JSK
    if (-60 >= seconds || seconds >= 60) {
      minutes += seconds / 60;
      seconds %= 60;
    }

    if (-60 >= minutes || minutes >= 60) {
      hour += minutes / 60;
      minutes %= 60;
    }

    if (-24 >= hour || hour >= 24) {
      days += hour / 24;
      hour %= 24;
    }
  }

  Time operator-(const Time& other) const {
    return Time(days - other.days,
                hour - other.hour,
                minutes - other.minutes,
                seconds - other.seconds);
  };

  Time operator+(const Time& other) const {
    return Time(days + other.days,
                hour + other.hour,
                minutes + other.minutes,
                seconds + other.seconds);
  };

  double decimalDay() const {
    return days + (((((seconds/60.0) + minutes)/60.0) + hour) / 24.0);
  };
};

// > In many books we read 'Julian Date' instead of 'Julian Day'.
// > For ys, a Julian date is a date in the Julian Calendar, just
// > as a Gregorian Date refers to the Gregorian Calendar. The JD
// > has nothing to do with the Julian Calendar.
//
// Meeus, Jean. Astronomical Algorithms. First English Edition,
//    Willmann-Bell Inc, 1991, pp59
struct JulianDay {
  double jd;

  // Meeus, Jean. Astronomical Algorithms. First English Edition,
  //    Willmann-Bell Inc, 1991, pp65
  int dayOfWeek() const {
    auto jd = (long)(this->jd + 1.5);
    return jd % 7;
  };

  std::string  namedDayOfWeekShort() const {
    switch(dayOfWeek()) {
      case 0: return "Sun";
      case 1: return "Mon";
      case 2: return "Tue";
      case 3: return "Wed";
      case 4: return "Thu";
      case 5: return "Fri";
      case 6: return "Sat";
      default: return "Unk";
    };
  }

  Time toTime() const {
    auto j = jd;
    long d = (long)j;
    j -= d;
    j *= 24;
    int h = (int)j;
    j -= h;
    j *= 60;
    int m = (int)j;
    j -= m;
    j *= 60;
    int s = (int)j;
    return Time{d, h, m, s};
  };
};


struct GregorianDate {
  int year;
  unsigned int month;
  // Float is ~6-7 digits, so 2 for the day and 4-5 for the decimal (time),
  // which doesn't seem precise enough? 1h/3600s*1d/24h = .00001157407407407
  // log_{10}(1h/3600s*1d/24h) = -4.936, so I don't know that it'd work
  // super well.
  // - JSK
  double day;

  // Meeus, Jean. Astronomical Algorithms. First English Edition,
  //    Willmann-Bell Inc, 1991, pp60-61
  //
  // I'm not sure how to get this to return a JulianDay. If I forward
  // declare it with `struct JulianDay;`, g++ says I can't use the type
  // as a return type, nor can I construct it at the bottom of the
  // function. - JSK
  double toJulianDay() const {
    long y = year;
    long m = month;
    double d = day;
    if (m <= 2) {
      y -= 1;
      m += 12;
    }

    long a = y/100;
    long b = 2 - a + (a/4);

    double jd = (long)(365.25 * (y + 4716)) + (long)(30.6001 * (m+1)) + (d + b - 1524.5);

    return jd;
  };

  std::string monthNameShort() {
    switch (month) {
      case  1: return "Jan";
      case  2: return "Feb";
      case  3: return "Mar";
      case  4: return "Apr";
      case  5: return "May";
      case  6: return "Jun";
      case  7: return "Jul";
      case  8: return "Aug";
      case  9: return "Sep";
      case 10: return "Oct";
      case 11: return "Nov";
      case 12: return "Dec";
      default: return "Unk";
    };
  };

  // Meeus, Jean. Astronomical Algorithms. First English Edition,
  //    Willmann-Bell Inc, 1991, pp62
  bool isLeapYear() const {
    if (year % 4 == 0) {
      if (year % 100 == 0) {
        if (year % 400 == 0) {
          return true;
        }
        return false;
      }
      return true;
    }
    return false;
  };

  // Meeus, Jean. Astronomical Algorithms. First English Edition,
  //    Willmann-Bell Inc, 1991, pp65
  int dayOfYear() const {
    auto k = 2;
    if (isLeapYear()) {
      k = 1;
    }
    auto n = (275 * month)/9 - (k * ((month + 9)/12)) + day - 30;
    return n;
  };

  // Meeus, Jean. Astronomical Algorithms. First English Edition,
  //    Willmann-Bell Inc, 1991, pp63
  //
  // N.B.: "The following method is valid for positive as well as for
  // negative years, but not for negative Julian Day numbers."
  GregorianDate static toGregorianDate(const JulianDay& j) {
    auto jd = j.jd;

    jd += 0.5;
    long z = (long)jd;
    auto f = jd - z;
    auto a = z;

    if (z >= 2299161) {
      auto α = (long)((z - 1867216.25)/36524.25);
      a = z + 1 + α - (long)(α / 4.0);
    }

    auto b = a + 1524;
    auto c = (long)((b - 122.1)/365.25);
    auto d = (long)(365.25 * c);
    auto e = (long)((b - d)/30.6001);

    auto day = (double)(b - d - (long)(30.6001 * e) + f);
    auto month = (unsigned int)(e - 1);
    if (e >= 14) {
      month = e - 13;
    }
    auto year = (int)c - 4716;
    if (month <= 2) {
      year = c - 4715;
    }

    return GregorianDate{year, month, day};
  };
};
