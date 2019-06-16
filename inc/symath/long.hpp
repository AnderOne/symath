/**
 * Copyright (c) 2019 Andrey Baranov <armath123@gmail.com>
 *
 * This file is part of SyMath++.
 *
 * SyMath++ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * SyMath++ is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SyMath++;
 * if not, see <http://www.gnu.org/licenses/>
**/

#ifndef __INCLUDE_LONG_H
#define __INCLUDE_LONG_H

#include <gmpxx.h>
#include <string>

namespace symath {

#define __DEF_MATH_OPERATOR(TYPE, op) \
TYPE operator op(const TYPE &rhs) const { return TYPE(val op rhs.val); }\
template <typename T>\
TYPE operator op(const T &rhs) const { return TYPE(val op rhs); }

#define __DEF_SET_OPERATOR(TYPE, op) \
template <typename T>\
TYPE operator op##=(const T &rhs) { return *this = *this op rhs; }

#define __DEF_CMP_OPERATOR(TYPE, op) \
bool operator op(const TYPE &rhs) const { return val op rhs.val; }\
template <typename T>\
bool operator op(const T &rhs) const { return val op rhs; }

#define __DEF_OPERATOR(TYPE) \
TYPE operator-() const { return - val; }\
__DEF_MATH_OPERATOR(TYPE, -)\
__DEF_MATH_OPERATOR(TYPE, +)\
__DEF_MATH_OPERATOR(TYPE, *)\
__DEF_MATH_OPERATOR(TYPE, /)\
__DEF_SET_OPERATOR(TYPE,  -)\
__DEF_SET_OPERATOR(TYPE,  +)\
__DEF_SET_OPERATOR(TYPE,  *)\
__DEF_SET_OPERATOR(TYPE,  /)\
__DEF_CMP_OPERATOR(TYPE, ==)\
__DEF_CMP_OPERATOR(TYPE, !=)\
__DEF_CMP_OPERATOR(TYPE, <=)\
__DEF_CMP_OPERATOR(TYPE, >=)\
__DEF_CMP_OPERATOR(TYPE, <)\
__DEF_CMP_OPERATOR(TYPE, >)

struct t_long {

	friend std::ostream &operator << (std::ostream &out, const t_long &src);
	friend std::istream &operator >> (std::istream &inp, t_long &dst);
	friend t_long pow(const t_long &src, unsigned deg);
	friend t_long abs(const t_long &);

	inline operator std::string() const { return val.get_str(); }
	inline operator long() const { return val.get_si(); }
	inline mpz_class mpz() const { return val; }

	inline t_long(const t_long &src): val(src.val) {}
	template <typename ... T>
	inline t_long(const T &... arg): val(arg ...) {}
	inline t_long(mpz_class mpz): val(mpz) {}
	inline t_long() {}

	__DEF_MATH_OPERATOR(t_long, %)
	__DEF_SET_OPERATOR(t_long, %)
	__DEF_OPERATOR(t_long)

private:
	mpz_class val;
};

struct t_frac {

	friend std::ostream &operator << (std::ostream &out, const t_frac &src);
	friend std::istream &operator >> (std::istream &inp, t_frac &dst);
	friend t_frac pow(const t_frac &src, long deg);
	friend t_frac abs(const t_frac &);

	inline operator std::string() const { return val.get_str(); }
	inline operator double() const { return val.get_d(); }
	inline mpq_class mpq() const { return val; }

	inline t_long upper() const { return val.get_num(); }
	inline t_long lower() const { return val.get_den(); }
	inline bool iszero() const {
		return val.get_num() == 0;
	}
	inline bool isint() const {
		return val.get_den() == 1;
	}

	inline t_frac(const t_long &num, const t_long &den):
	              val(num.mpz(), den.mpz()) {}
	inline t_frac(const t_frac &src): val(src.val) {}
	template <typename ... T>
	inline t_frac(const T &... arg): val(arg ...) {
		val.canonicalize();
	}
	inline t_frac(mpq_class mpq): val(mpq) {
		val.canonicalize();
	}
	inline t_frac() {}

	__DEF_OPERATOR(t_frac)

private:
	mpq_class val;
};

#undef __DEF_MATH_OPERATOR
#undef __DEF_CMP_OPERATOR
#undef __DEF_SET_OPERATOR
#undef __DEF_OPERATOR

inline std::ostream &operator << (std::ostream &out, const t_long &src) {
	return out << src.val;
}

inline std::ostream &operator << (std::ostream &out, const t_frac &src) {
	return out << src.val;
}

inline std::istream &operator >> (std::istream &inp, t_long &dst) {
	return inp >> dst.val;
}

inline std::istream &operator >> (std::istream &inp, t_frac &dst) {
	return inp >> dst.val;
}

inline t_long pow(const t_long &src, unsigned n) {
	mpz_class b = src.val, a = 1;
	while (n) {
		if (n & 1) a *= b;
		b *= b;
		n /= 2;
	}
	return a;
}

inline t_frac pow(const t_frac &src, long n) {
	t_frac a = t_frac(
		pow(src.upper(), abs(n)), pow(src.lower(), abs(n))
	);
	return (n < 0)? (t_frac(1) / a): (a);
}

inline t_long abs(const t_long &src) {
	return abs(src.val);
}

inline t_frac abs(const t_frac &src) {
	return abs(src.val);
}

//Format of conversion to string:
enum t_form {
FRM_RAT = 0, FRM_RED = 1, FRM_DOT = 2, FRM_ALL = 3
};

std::string format(
const t_frac &val, t_form frm = FRM_ALL
);

//...

}

#endif //__INCLUDE_LONG_H
