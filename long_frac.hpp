#ifndef __INCLUDE_LONG_FRAC_H
#define __INCLUDE_LONG_FRAC_H

#include <gmpxx.h>
#include <string>

typedef mpz_class t_long;

struct t_long_frac {

	friend std::ostream &operator << (std::ostream &out, const t_long_frac &src);
	friend std::istream &operator >> (std::istream &inp, t_long_frac &dst);

	inline operator std::string() const { return val.get_str(); }
	inline operator double() const { return val.get_d(); }

	inline const t_long &upper() const {
		return val.get_num();
	}
	inline const t_long &lower() const {
		return val.get_den();
	}
	inline bool iszero() const {
		return val.get_num() == 0;
	}
	inline bool isint() const {
		return val.get_den() == 1;
	}

	inline t_long_frac(const t_long_frac &src): val(src.val) {}
	template <typename ... T>
	inline t_long_frac(const T &... arg): val(arg ...) {
		val.canonicalize();
	}
	inline t_long_frac(mpq_class mpq): val(mpq) {
		val.canonicalize();
	}
	inline t_long_frac() {}

	#define DEF_MATH_OPERATOR(op) \
	t_long_frac operator op(const t_long_frac &rhs) const {\
		return t_long_frac(val op rhs.val);\
	}\
	template <typename T>\
	t_long_frac operator op(const T &rhs) const {\
		return t_long_frac(val op rhs);\
	}

	#define DEF_CMP_OPERATOR(op) \
	bool operator op(const t_long_frac &rhs) const {\
		return val op rhs.val;\
	}\
	template <typename T>\
	bool operator op(const T &rhs) const {\
		return val op rhs;\
	}

	#define DEF_SET_OPERATOR(op) \
	template <typename T>\
	t_long_frac operator op##= (const T &rhs) {\
		return (*this) = (*this) op rhs;\
	}

	t_long_frac operator-() const {
		return - val;
	}
	DEF_MATH_OPERATOR(-)
	DEF_MATH_OPERATOR(+)
	DEF_MATH_OPERATOR(*)
	DEF_MATH_OPERATOR(/)

	DEF_CMP_OPERATOR(==)
	DEF_CMP_OPERATOR(!=)
	DEF_CMP_OPERATOR(<=)
	DEF_CMP_OPERATOR(>=)
	DEF_CMP_OPERATOR(<)
	DEF_CMP_OPERATOR(>)

	DEF_SET_OPERATOR(-)
	DEF_SET_OPERATOR(+)
	DEF_SET_OPERATOR(*)
	DEF_SET_OPERATOR(/)

	#undef DEF_MATH_OPERATOR
	#undef DEF_CMP_OPERATOR

private:
	friend t_long_frac abs(const t_long_frac &);
	mpq_class val;
};

inline std::ostream &operator << (std::ostream &out, const t_long_frac &src) {
	return out << src.val;
}

inline std::istream &operator >> (std::istream &inp, t_long_frac &dst) {
	return inp >> dst.val;
}

inline t_long_frac abs(const t_long_frac &src) {
	return abs(src.val);
}

#endif //__INCLUDE_LONG_FRAC_H
