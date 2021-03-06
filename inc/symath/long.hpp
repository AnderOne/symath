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

#define __DEF_RIGHT_MATH_OPERATOR(TYPE, op) \
template <typename L>\
TYPE operator op(const L &lhs, const TYPE &rhs) { return TYPE(lhs) op rhs; }

#define __DEF_RIGHT_CMP_OPERATOR(TYPE, op) \
template <typename L>\
bool operator op(const L &lhs, const TYPE &rhs) { return TYPE(lhs) op rhs; }

#define __DEF_RIGHT_OPERATOR(TYPE) \
__DEF_RIGHT_MATH_OPERATOR(TYPE, -)\
__DEF_RIGHT_MATH_OPERATOR(TYPE, +)\
__DEF_RIGHT_MATH_OPERATOR(TYPE, *)\
__DEF_RIGHT_MATH_OPERATOR(TYPE, /)\
__DEF_RIGHT_CMP_OPERATOR(TYPE, ==)\
__DEF_RIGHT_CMP_OPERATOR(TYPE, !=)\
__DEF_RIGHT_CMP_OPERATOR(TYPE, <=)\
__DEF_RIGHT_CMP_OPERATOR(TYPE, >=)\
__DEF_RIGHT_CMP_OPERATOR(TYPE, <)\
__DEF_RIGHT_CMP_OPERATOR(TYPE, >)

#define __DEF_LEFT_MATH_OPERATOR(TYPE, op) \
TYPE operator op(const TYPE &rhs) const { return TYPE(val op rhs.val); }\
template <typename T>\
TYPE operator op(const T &rhs) const { return TYPE(val op rhs); }

#define __DEF_LEFT_SET_OPERATOR(TYPE, op) \
template <typename T>\
TYPE operator op##=(const T &rhs) { return *this = *this op rhs; }

#define __DEF_LEFT_CMP_OPERATOR(TYPE, op) \
bool operator op(const TYPE &rhs) const { return val op rhs.val; }\
template <typename T>\
bool operator op(const T &rhs) const { return val op rhs; }

#define __DEF_LEFT_OPERATOR(TYPE) \
TYPE operator-() const { return - val; }\
__DEF_LEFT_MATH_OPERATOR(TYPE, -)\
__DEF_LEFT_MATH_OPERATOR(TYPE, +)\
__DEF_LEFT_MATH_OPERATOR(TYPE, *)\
__DEF_LEFT_MATH_OPERATOR(TYPE, /)\
__DEF_LEFT_SET_OPERATOR(TYPE,  -)\
__DEF_LEFT_SET_OPERATOR(TYPE,  +)\
__DEF_LEFT_SET_OPERATOR(TYPE,  *)\
__DEF_LEFT_SET_OPERATOR(TYPE,  /)\
__DEF_LEFT_CMP_OPERATOR(TYPE, ==)\
__DEF_LEFT_CMP_OPERATOR(TYPE, !=)\
__DEF_LEFT_CMP_OPERATOR(TYPE, <=)\
__DEF_LEFT_CMP_OPERATOR(TYPE, >=)\
__DEF_LEFT_CMP_OPERATOR(TYPE, <)\
__DEF_LEFT_CMP_OPERATOR(TYPE, >)

struct t_long {

	friend std::ostream &operator << (std::ostream &out, const t_long &src);
	friend std::istream &operator >> (std::istream &inp, t_long &dst);

	inline operator std::string() const { return val.get_str(); }
	inline operator mpz_class() const { return val; }

	inline std::string str() const { return val.get_str(); }
	inline long get() const { return val.get_si(); }

	t_long abs() const { return (val < 0)? (- val): (val); }
	t_long pow(unsigned n) const;

	inline t_long(const t_long &src): val(src.val) {}
	template <typename ... T>
	inline t_long(const T &... arg): val(arg ...) {}
	inline t_long(mpz_class mpz): val(mpz) {}
	inline t_long() {}

	__DEF_LEFT_MATH_OPERATOR(t_long, %)
	__DEF_LEFT_SET_OPERATOR(t_long, %)
	__DEF_LEFT_OPERATOR(t_long)

private:
	mpz_class val;
};

__DEF_RIGHT_MATH_OPERATOR(t_long, %)
__DEF_RIGHT_OPERATOR(t_long)

struct t_frac {

	friend std::ostream &operator << (std::ostream &out, const t_frac &src);
	friend std::istream &operator >> (std::istream &inp, t_frac &dst);

	//Format of conversion to string:
	enum t_form { FRM_RAT = 0, FRM_RED = 1, FRM_DOT = 2, FRM_ALL = 3 };

	inline operator std::string() const { return val.get_str(); }
	inline operator mpq_class() const { return val; }
	inline double get() const { return val.get_d(); }
	std::string str(t_form frm = FRM_ALL);

	inline t_long upper() const { return val.get_num(); }
	inline t_long lower() const { return val.get_den(); }
	inline bool iszero() const {
		return val.get_num() == 0;
	}
	inline bool isint() const {
		return val.get_den() == 1;
	}

	t_frac abs() const { return (val < 0)? (- val): (val); }
	t_frac pow(long n) const;

	inline t_frac(const t_long &num, const t_long &den):
	                                  val(num, den) {}
	inline t_frac(const t_frac &src): val(src.val) {}
	template <typename ... T>
	inline t_frac(const T &... arg): val(arg ...) {
		val.canonicalize();
	}
	inline t_frac(mpq_class mpq): val(mpq) {
		val.canonicalize();
	}
	inline t_frac() {}

	__DEF_LEFT_OPERATOR(t_frac)

private:
	mpq_class val;
};

__DEF_RIGHT_OPERATOR(t_frac)

#undef __DEF_RIGHT_MATH_OPERATOR
#undef __DEF_RIGHT_CMP_OPERATOR
#undef __DEF_RIGHT_OPERATOR

#undef __DEF_LEFT_MATH_OPERATOR
#undef __DEF_LEFT_CMP_OPERATOR
#undef __DEF_LEFT_SET_OPERATOR
#undef __DEF_LEFT_OPERATOR

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
	mpq_class val;
	inp >> val;
	dst = t_frac(val);
	return inp;
}

//...

}

#endif //__INCLUDE_LONG_H
