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

#ifndef __INCLUDE_TREE
#define __INCLUDE_TREE

#include "long.hpp"

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

struct t_tree {

	inline const double &operator[] (char var) const { return DATA[(var - 'a')]; }
	inline double &operator[] (char var) { return DATA[(var - 'a')]; }

	inline operator std::string() const { return root->str(); }
	inline double get() const { return root->get(); }

	t_tree &operator=(const t_tree &rhs);

	t_tree(const t_tree &src);
	t_tree(std::string str);
	t_tree();

	bool create(std::string str);
	void derive(char var);
	void reduce();

protected:

	struct h_node;
	struct t_node;

	//Фабричные методы для создания разных типов узлов:
	h_node gener(std::string str, const h_node &lhs, const h_node &rhs) const;
	h_node gener(std::string str, const h_node &rhs) const;
	h_node gener(std::string str) const;
	h_node gener(t_frac val) const;
	h_node gener(long val) const;

	struct h_node: public std::shared_ptr<t_node> {

		template <typename ... T>
		h_node(T ... args): std::shared_ptr<t_node> (args ...) {}

		#define __DEF_ITEM_BIN(p) \
		inline h_node operator p(const h_node &rhs) const {\
			return get()->own.gener(#p, *this, rhs);\
		}\
		inline h_node operator p(t_frac val) const {\
			return get()->own.gener(\
			#p, *this, get()->own.gener(val)\
			);\
		}
		inline h_node operator-() const {
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
	inline h_node p(const h_node &arg) const { return gener(#p, arg); }

	__DEF_ITEM_ONE(sqrt)
	__DEF_ITEM_ONE(exp)
	__DEF_ITEM_ONE(log)
	__DEF_ITEM_ONE(cos)
	__DEF_ITEM_ONE(sin)

	#undef __DEF_ITEM_ONE

	struct t_node {
		explicit t_node(const t_tree &_own, t_frac _num = 1): own(_own), num(_num) {}
		virtual ~t_node();
		virtual std::string str() const = 0;	//Переводит выражение в строку;
		virtual h_node dif(char) const = 0;	//Возвращает производную;
		virtual h_node red() const = 0;	//Сокращает выражение;
		virtual double get() const = 0;	//Вычисляет значение;
		virtual h_node cpy(const t_tree &) const = 0;
		virtual h_node cpy() const = 0;
		//...
		h_node mul(t_frac _fac) { h_node ptr = ref(); ptr->num *= _fac; return ptr; }
		h_node ref() { return own.LINK[this].lock(); }
		//...
		const t_tree &own;
		t_frac num;
	};

	struct t_node_num:
	public t_node {
		explicit t_node_num(const t_tree &_own, t_frac _num):
		         t_node{_own, _num} {}
		std::string str() const override;
		h_node dif(char) const override;
		h_node red() const override;
		double get() const override;
		h_node cpy(const t_tree &) const override;
		h_node cpy() const override;
	};

	struct t_node_var:
	public t_node {
		explicit t_node_var(const t_tree &_own, char _ind):
		         t_node(_own), ind(_ind) {}
		std::string str() const override;
		h_node dif(char) const override;
		h_node red() const override;
		double get() const override;
		h_node cpy(const t_tree &) const override;
		h_node cpy() const override;
		const char ind;
	};

	struct t_node_bin:
	public t_node {
		explicit t_node_bin(const t_tree &_own, const h_node &_lhs, const h_node &_rhs):
		         t_node(_own),arg{_lhs, _rhs} {}
		virtual h_node split(t_frac &_num) const = 0;
		h_node arg[2];
	};

	#define __DEF_ITEM_BIN(p) \
	struct t_node_##p:\
	public t_node_bin {\
		explicit t_node_##p(const t_tree &_own, const h_node &_lhs, const h_node &_rhs):\
		         t_node_bin(_own, _lhs, _rhs) {}\
		h_node split(t_frac &_num) const override;\
		std::string str() const override;\
		h_node dif(char) const override;\
		h_node red() const override;\
		double get() const override;\
		h_node cpy(const t_tree &)\
		const override;\
		h_node cpy() const override;\
	};

	__DEF_ITEM_BIN(add)
	__DEF_ITEM_BIN(mul)
	__DEF_ITEM_BIN(div)
	__DEF_ITEM_BIN(pow)

	#undef __DEF_ITEM_BIN

	struct t_node_one:
	public t_node {
		explicit t_node_one(const t_tree &_own, const h_node &_arg):
		         t_node(_own), arg{_arg} {}
		h_node arg;
	};

	#define __DEF_ITEM_ONE(p) \
	struct t_node_##p:\
	public t_node_one {\
		explicit t_node_##p(const t_tree &_own, const h_node &_arg):\
		         t_node_one(_own, _arg) {}\
		std::string str() const override;\
		h_node dif(char) const override;\
		h_node red() const override;\
		double get() const override;\
		h_node cpy(const t_tree &)\
		const override;\
		h_node cpy() const override;\
	};

	__DEF_ITEM_ONE(exp)
	__DEF_ITEM_ONE(log)
	__DEF_ITEM_ONE(cos)
	__DEF_ITEM_ONE(sin)

	#undef __DEF_ITEM_ONE

	friend struct t_node;

private:
	h_node store(t_node *&& _ptr) const;

	mutable std::map<
	const t_node *,
	std::weak_ptr<t_node> > LINK;

	double DATA[26];

	h_node root;
};

inline std::ostream &operator<<(std::ostream &out, const t_tree &tree) {
	return out << (std::string) tree;
}

inline std::istream &operator>>(std::istream &inp, t_tree &tree) {
	std::string str;
	std::getline(inp, str);
	tree.create(str);
	return inp;
}

#endif //__INCLUDE_TREE
