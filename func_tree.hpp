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

#ifndef __INCLUDE_FUNC_TREE
#define __INCLUDE_FUNC_TREE

#include "long_frac.hpp"

#include <type_traits>
#include <exception>
#include <cassert>
#include <ostream>
#include <istream>
#include <memory>
#include <regex>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <stack>
#include <map>
#include <set>

struct t_func_tree {

	inline const double &operator[] (char var) const { return DATA[(var - 'a')]; }
	inline double &operator[] (char var) { return DATA[(var - 'a')]; }

	inline operator std::string() const { return root->str(); }
	inline double get() const { return root->get(); }

	t_func_tree &operator=(const t_func_tree &rhs);

	t_func_tree(const t_func_tree &src);
	t_func_tree(std::string str);
	t_func_tree();

	bool create(std::string str);
	void derive(char var);
	void reduce();

protected:

	struct h_item;
	struct t_item;

	//Фабричные методы для создания разных типов узлов:
	h_item gener(std::string str, const h_item &lhs, const h_item &rhs) const;
	h_item gener(std::string str, const h_item &rhs) const;
	h_item gener(std::string str) const;
	h_item gener(t_long_frac val) const;
	h_item gener(long val) const;

	struct h_item: public std::shared_ptr<t_item> {

		template <typename ... T>
		h_item(T ... args): std::shared_ptr<t_item> (args ...) {}

		#define __DEF_ITEM_BIN(p) \
		inline h_item operator p(const h_item &rhs) const {\
			return get()->own.gener(#p, *this, rhs);\
		}\
		inline h_item operator p(t_long_frac val) const {\
			return get()->own.gener(\
			#p, *this, get()->own.gener(val)\
			);\
		}
		inline h_item operator-() const {
			return get()->own.gener("-", *this);
		}

		__DEF_ITEM_BIN(-)
		__DEF_ITEM_BIN(+)
		__DEF_ITEM_BIN(*)
		__DEF_ITEM_BIN(/)
		__DEF_ITEM_BIN(^)

		#undef __DEF_ITEM_BIN

	};

	//Математические функции:
	#define __DEF_ITEM_ONE(p) \
	inline h_item p(const h_item &arg) const { return gener(#p, arg); }

	__DEF_ITEM_ONE(sqrt)
	__DEF_ITEM_ONE(exp)
	__DEF_ITEM_ONE(log)
	__DEF_ITEM_ONE(cos)
	__DEF_ITEM_ONE(sin)

	#undef __DEF_ITEM_ONE

	struct t_item {
		explicit t_item(const t_func_tree &_own, t_long_frac _num = 1): own(_own), num(_num) {}
		virtual ~t_item();
		virtual std::string str() const = 0;	//Переводит выражение в строку;
		virtual h_item dif(char) const = 0;	//Возвращает производную;
		virtual h_item red() const = 0;	//Сокращает выражение;
		virtual double get() const = 0;	//Вычисляет значение;
		virtual h_item cpy(const t_func_tree &) const = 0;
		virtual h_item cpy() const = 0;
		//...
		h_item mul(t_long_frac _fac) { h_item ptr = ref(); ptr->num *= _fac; return ptr; }
		h_item ref() { return own.LINK[this].lock(); }
		//...
		const t_func_tree &own;
		t_long_frac num;
	};

	struct t_item_num:
	public t_item {
		explicit t_item_num(const t_func_tree &_own, t_long_frac _num):
		         t_item{_own, _num} {}
		std::string str() const override;
		h_item dif(char) const override;
		h_item red() const override;
		double get() const override;
		h_item cpy(const t_func_tree &) const override;
		h_item cpy() const override;
	};

	struct t_item_var:
	public t_item {
		explicit t_item_var(const t_func_tree &_own, char _ind):
		         t_item(_own), ind(_ind) {}
		std::string str() const override;
		h_item dif(char) const override;
		h_item red() const override;
		double get() const override;
		h_item cpy(const t_func_tree &) const override;
		h_item cpy() const override;
		const char ind;
	};

	struct t_item_bin:
	public t_item {
		explicit t_item_bin(const t_func_tree &_own, const h_item &_lhs, const h_item &_rhs):
		         t_item(_own),arg{_lhs, _rhs} {}
		virtual h_item split(t_long_frac &_num) const = 0;
		h_item arg[2];
	};

	#define __DEF_ITEM_BIN(p) \
	struct t_item_##p:\
	public t_item_bin {\
		explicit t_item_##p(const t_func_tree &_own, const h_item &_lhs, const h_item &_rhs):\
		         t_item_bin(_own, _lhs, _rhs) {}\
		h_item split(t_long_frac &_num) const override;\
		std::string str() const override;\
		h_item dif(char) const override;\
		h_item red() const override;\
		double get() const override;\
		h_item cpy(const t_func_tree &)\
		const override;\
		h_item cpy() const override;\
	};

	__DEF_ITEM_BIN(add)
	__DEF_ITEM_BIN(mul)
	__DEF_ITEM_BIN(div)
	__DEF_ITEM_BIN(pow)

	#undef __DEF_ITEM_BIN

	struct t_item_one:
	public t_item {
		explicit t_item_one(const t_func_tree &_own, const h_item &_arg):
		         t_item(_own), arg{_arg} {}
		h_item arg;
	};

	#define __DEF_ITEM_ONE(p) \
	struct t_item_##p:\
	public t_item_one {\
		explicit t_item_##p(const t_func_tree &_own, const h_item &_arg):\
		         t_item_one(_own, _arg) {}\
		std::string str() const override;\
		h_item dif(char) const override;\
		h_item red() const override;\
		double get() const override;\
		h_item cpy(const t_func_tree &)\
		const override;\
		h_item cpy() const override;\
	};

	__DEF_ITEM_ONE(exp)
	__DEF_ITEM_ONE(log)
	__DEF_ITEM_ONE(cos)
	__DEF_ITEM_ONE(sin)

	#undef __DEF_ITEM_ONE

	friend struct t_item;

private:
	h_item store(t_item *&& _ptr) const;

	mutable std::map<
	const t_item *,
	std::weak_ptr<t_item> > LINK;

	double DATA[26];

	h_item root;
};

inline std::ostream &operator<<(std::ostream &out, const t_func_tree &tree) {
	return out << (std::string) tree;
}

inline std::istream &operator>>(std::istream &inp, t_func_tree &tree) {
	std::string str;
	std::getline(inp, str);
	tree.create(str);
	return inp;
}

#endif //__INCLUDE_FUNC_TREE
