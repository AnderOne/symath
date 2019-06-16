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

struct t_long: public mpz_class {
	template <typename ... T> inline t_long(const T &... arg): mpz_class(arg ...) {}
	inline operator std::string() const { return get_str(); }
	inline operator long() const { return get_si(); }
};

struct t_frac {

	friend std::ostream &operator << (std::ostream &out, const t_frac &src);
	friend std::istream &operator >> (std::istream &inp, t_frac &dst);
	friend t_frac pow(const t_frac &src, long deg);
	friend t_frac abs(const t_frac &);

	inline operator std::string() const { return val.get_str(); }
	inline operator double() const { return val.get_d(); }

	inline t_long upper() const { return val.get_num(); }
	inline t_long lower() const { return val.get_den(); }
	inline bool iszero() const {
		return val.get_num() == 0;
	}
	inline bool isint() const {
		return val.get_den() == 1;
	}

	inline t_frac(const t_frac &src): val(src.val) {}
	template <typename ... T>
	inline t_frac(const T &... arg): val(arg ...) {
		val.canonicalize();
	}
	inline t_frac(mpq_class mpq): val(mpq) {
		val.canonicalize();
	}
	inline t_frac() {}

	#define __DEF_MATH_OPERATOR(op) \
	t_frac operator op(const t_frac &rhs) const {\
		return t_frac(val op rhs.val);\
	}\
	template <typename T>\
	t_frac operator op(const T &rhs) const {\
		return t_frac(val op rhs);\
	}

	#define __DEF_CMP_OPERATOR(op) \
	bool operator op(const t_frac &rhs) const {\
		return val op rhs.val;\
	}\
	template <typename T>\
	bool operator op(const T &rhs) const {\
		return val op rhs;\
	}

	#define __DEF_SET_OPERATOR(op) \
	template <typename T>\
	t_frac operator op##= (const T &rhs) {\
		return (*this) = (*this) op rhs;\
	}

	t_frac operator-() const {
		return - val;
	}
	__DEF_MATH_OPERATOR(-)
	__DEF_MATH_OPERATOR(+)
	__DEF_MATH_OPERATOR(*)
	__DEF_MATH_OPERATOR(/)

	__DEF_CMP_OPERATOR(==)
	__DEF_CMP_OPERATOR(!=)
	__DEF_CMP_OPERATOR(<=)
	__DEF_CMP_OPERATOR(>=)
	__DEF_CMP_OPERATOR(<)
	__DEF_CMP_OPERATOR(>)

	__DEF_SET_OPERATOR(-)
	__DEF_SET_OPERATOR(+)
	__DEF_SET_OPERATOR(*)
	__DEF_SET_OPERATOR(/)

	#undef __DEF_MATH_OPERATOR
	#undef __DEF_CMP_OPERATOR
	#undef __DEF_SET_OPERATOR

private:
	mpq_class val;
};

inline std::ostream &operator << (std::ostream &out, const t_frac &src) {
	return out << src.val;
}

inline std::istream &operator >> (std::istream &inp, t_frac &dst) {
	return inp >> dst.val;
}

inline t_frac pow(const t_frac &src, long deg) {
	unsigned long n = std::abs(deg);
	mpq_class b = src.val, a = 1;
	while (n) {
		if (n & 1) a *= b;
		b *= b;
		n /= 2;
	}
	if (deg < 0) a = 1 / a;
	return t_frac(a);
}

inline t_frac abs(const t_frac &src) {
	return abs(src.val);
}

#endif //__INCLUDE_LONG_H
