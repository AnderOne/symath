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

#include "func_tree.hpp"

#include <exception>
#include <cassert>
#include <cmath>
#include <regex>
#include <algorithm>
#include <vector>
#include <string>
#include <stack>
#include <set>

//Методы клонирования узлов:

t_func_tree::h_item t_func_tree::t_item_var::cpy(const t_func_tree &_own) const { return _own.gener(std::string(1, ind))->mul(num); }
t_func_tree::h_item t_func_tree::t_item_var::cpy() const { return own.gener(std::string(1, ind))->mul(num); }

t_func_tree::h_item t_func_tree::t_item_num::cpy(const t_func_tree &_own) const { return _own.gener(num); }
t_func_tree::h_item t_func_tree::t_item_num::cpy() const { return own.gener(num); }

#define __DEF_ITEM_BIN_CPY(ind, op) \
t_func_tree::h_item t_func_tree::t_item_##ind::cpy(const t_func_tree &_own) const { return _own.gener(#op, arg[0], arg[1])->mul(num); }\
t_func_tree::h_item t_func_tree::t_item_##ind::cpy() const { return own.gener(#op, arg[0], arg[1])->mul(num); }

__DEF_ITEM_BIN_CPY(add, +)
__DEF_ITEM_BIN_CPY(mul, *)
__DEF_ITEM_BIN_CPY(div, /)
__DEF_ITEM_BIN_CPY(pow, ^)

#undef __DEF_ITEM_BIN_CPY

#define __DEF_ITEM_ONE_CPY(ind) \
t_func_tree::h_item t_func_tree::t_item_##ind::cpy(const t_func_tree &_own) const { return _own.gener(#ind, arg)->mul(num); }\
t_func_tree::h_item t_func_tree::t_item_##ind::cpy() const { return own.gener(#ind, arg)->mul(num); }

__DEF_ITEM_ONE_CPY(exp)
__DEF_ITEM_ONE_CPY(log)
__DEF_ITEM_ONE_CPY(cos)
__DEF_ITEM_ONE_CPY(sin)

#undef __DEF_ITEM_ONE_CPY

//Перевод в строку:

std::string t_func_tree::t_item_var::str() const {
	return ((num < 0)? ("- "): ("")) + ((abs(num) != 1)? (std::string(abs(num)) + " * "): ("")) + std::string(1, ind);
}

std::string t_func_tree::t_item_num::str() const {
	return ((num < 0)? ("- "): ("")) + std::string(abs(num));
}

std::string t_func_tree::t_item_add::str() const {
	std::string tmp;
	tmp = ((num < 0)? ("- "): ("")) + ((abs(num) != 1)? (std::string(abs(num)) + " * "): (""));
	tmp += "(";
	tmp += arg[0]->str();
	tmp += (arg[1]->num >= 0)?
	       (" + "):
	       (" ");
	tmp += arg[1]->str();
	tmp += ")";
	return tmp;
}

#define __DEF_ITEM_BIN_STR(ind, op) \
std::string t_func_tree::t_item_##ind::str() const {\
	std::string tmp;\
	tmp = ((num < 0)? ("- "): ("")) + ((abs(num) != 1)? (std::string(abs(num)) + " * "): (""));\
	tmp += "((";\
	tmp += arg[0]->str();\
	tmp += ") "#op" (";\
	tmp += arg[1]->str();\
	tmp += "))";\
	return tmp;\
}

__DEF_ITEM_BIN_STR(mul, *)
__DEF_ITEM_BIN_STR(div, /)
__DEF_ITEM_BIN_STR(pow, ^)

#undef __DEF_ITEM_BIN_STR

#define __DEF_ITEM_ONE_STR(ind) \
std::string t_func_tree::t_item_##ind::str() const {\
	std::string tmp;\
	tmp = ((num < 0)? ("- "): ("")) + ((abs(num) != 1)? (std::string(abs(num)) + " * "): (""));\
	tmp += #ind"(";\
	tmp += arg->str();\
	tmp += ")";\
	return tmp;\
}

__DEF_ITEM_ONE_STR(exp)
__DEF_ITEM_ONE_STR(log)
__DEF_ITEM_ONE_STR(cos)
__DEF_ITEM_ONE_STR(sin)

#undef __DEF_ITEM_ONE_STR

//Вычисление выражения:

double t_func_tree::t_item_var::get() const { return num * own[ind]; }
double t_func_tree::t_item_num::get() const { return num; }

double t_func_tree::t_item_pow::get() const {
	return num * std::pow(arg[0]->get(), arg[1]->get());
}

#define __DEF_ITEM_BIN_GET(ind, op) \
double t_func_tree::t_item_##ind::get() const {\
	return num * (arg[0]->get() op arg[1]->get());\
}

__DEF_ITEM_BIN_GET(add, +)
__DEF_ITEM_BIN_GET(mul, *)
__DEF_ITEM_BIN_GET(div, /)

#undef __DEF_ITEM_BIN_GET

#define __DEF_ITEM_ONE_GET(ind) \
double t_func_tree::t_item_##ind::get() const {\
	return num * std::ind(arg->get());\
}

__DEF_ITEM_ONE_GET(exp)
__DEF_ITEM_ONE_GET(log)
__DEF_ITEM_ONE_GET(cos)
__DEF_ITEM_ONE_GET(sin)

#undef __DEF_ITEM_ONE_GET

//Методы дифференцирования узлов:

t_func_tree::h_item t_func_tree::t_item_pow::dif(char var) const {
	h_item hand = arg[1]->cpy() * arg[0]->dif(var) * (arg[0]->cpy() ^ (arg[1]->cpy() - own.gener(1))) +
	              arg[1]->dif(var) * own.log(arg[0]->cpy()) * cpy();
	return hand->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_div::dif(char var) const {
	h_item hand = arg[0]->dif(var) / arg[1]->cpy() - arg[0]->cpy() * arg[1]->dif(var) / arg[1]->cpy();
	return hand->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_mul::dif(char var) const {
	h_item hand = arg[0]->dif(var) * arg[1]->cpy() + arg[0]->cpy() * arg[1]->dif(var);
	return hand->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_add::dif(char var) const {
	h_item hand = arg[0]->dif(var) + arg[1]->dif(var);
	return hand->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_cos::dif(char var) const {
	h_item hand = - (arg->dif(var) * own.sin(arg->cpy()));
	return hand->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_sin::dif(char var) const {
	h_item hand = arg->dif(var) * own.cos(arg->cpy());
	return hand->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_log::dif(char var) const {
	h_item hand = arg->dif(var) / arg->cpy();
	return hand->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_exp::dif(char var) const {
	h_item hand = arg->dif(var) * cpy();
	return hand->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_var::dif(char var) const {
	h_item hand = own.gener((ind == var)? num: t_long_frac(0));
	return hand;
}

t_func_tree::h_item t_func_tree::t_item_num::dif(char var) const {
	h_item hand = own.gener(0);
	return hand;
}

//Методы сокращения выражений:

t_func_tree::h_item t_func_tree::t_item_add::split(t_long_frac &_num) const {

	t_item_num *ptr = dynamic_cast<t_item_num *> (arg[0].get());
	if (ptr) {
		_num = ptr->num;
		return arg[1];
	}
	ptr = dynamic_cast<t_item_num *> (arg[1].get());
	if (ptr) {
		_num = ptr->num;
		return arg[0];
	}
	return nullptr;
}

t_func_tree::h_item t_func_tree::t_item_mul::split(t_long_frac &_num) const {
	return nullptr;
}

t_func_tree::h_item t_func_tree::t_item_div::split(t_long_frac &_num) const {
	return nullptr;
}

t_func_tree::h_item t_func_tree::t_item_pow::split(t_long_frac &_num) const {
	return nullptr;
}

t_func_tree::h_item t_func_tree::t_item_var::red() const { return cpy(); }

t_func_tree::h_item t_func_tree::t_item_num::red() const { return cpy(); }

t_func_tree::h_item t_func_tree::t_item_add::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	t_item_num *ln = dynamic_cast<t_item_num *> (lhs.get());
	t_item_num *rn = dynamic_cast<t_item_num *> (rhs.get());
	if (ln && rn) {
		return own.gener(ln->num + rn->num)->mul(num);
	}
	if (ln && (ln->num == 0)) {
		return rhs->mul(num);
	}
	if (rn && (rn->num == 0)) {
		return lhs->mul(num);
	}

	t_item_bin *lb = dynamic_cast<t_item_bin *> (lhs.get());
	t_item_bin *rb = dynamic_cast<t_item_bin *> (rhs.get());
	t_long_frac lv, rv;
	h_item l2, r2;
	if (lb) l2 = lb->split(lv);
	if (rb) r2 = rb->split(rv);
	if (l2 && r2) {
		return (own.gener(lv * lb->num + rv * rb->num) + (l2->mul(lb->num) + r2->mul(rb->num)))->mul(num);
	}
	if (l2 && rn) {
		return (own.gener(lv * lb->num + rn->num) + l2->mul(lb->num))->mul(num);
	}
	if (ln && r2) {
		return (own.gener(ln->num + rv * rb->num) + r2->mul(rb->num))->mul(num);
	}

	return (lhs + rhs)->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_mul::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	t_item_num *ln = dynamic_cast<t_item_num *> (lhs.get());
	t_item_num *rn = dynamic_cast<t_item_num *> (rhs.get());
	if (ln && rn) {
		return own.gener(ln->num * rn->num)->mul(num);
	}
	if (ln && (std::abs(ln->num) == 1)) {
		return rhs->mul(ln->num * num);
	}
	if (rn && (std::abs(rn->num) == 1)) {
		return lhs->mul(rn->num * num);
	}
	if (ln && (ln->num == 0)) {
		return own.gener(0);
	}
	if (rn && (rn->num == 0)) {
		return own.gener(0);
	}

	t_item_bin *lb = dynamic_cast<t_item_bin *> (lhs.get());
	t_item_bin *rb = dynamic_cast<t_item_bin *> (rhs.get());
	t_long_frac lv, rv;
	h_item l2, r2;
	if (lb) l2 = lb->split(lv);
	if (rb) r2 = rb->split(rv);
	if (l2 && rn) {
		return (own.gener(lv * lb->num * rn->num) + l2->mul(lb->num * rn->num))->mul(num);
	}
	if (ln && r2) {
		return (own.gener(ln->num * rv * rb->num) + r2->mul(rb->num * ln->num))->mul(num);
	}
	if (ln) {
		return (rhs->mul(ln->num))->mul(num);
	}
	if (rn) {
		return (lhs->mul(rn->num))->mul(num);
	}
	t_long_frac n =
	     num * lhs->num * rhs->num;
	lhs->num = 1;
	rhs->num = 1;
	return (lhs * rhs)->mul(n);
}

t_func_tree::h_item t_func_tree::t_item_div::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	t_item_num *ln = dynamic_cast<t_item_num *> (lhs.get());
	t_item_num *rn = dynamic_cast<t_item_num *> (rhs.get());
	if (ln && (ln->num == 0)) {
		return own.gener(0);
	}
	if (ln && rn) {
		return own.gener(ln->num / rn->num)->mul(num);
	}

	t_item_bin *lb = dynamic_cast<t_item_bin *> (lhs.get());
	t_long_frac lv;
	h_item l2;
	if (lb) l2 = lb->split(lv);
	if (l2 && rn) {
		return (own.gener((lv * lb->num) / rn->num) + l2->mul(lb->num / rn->num))->mul(num);
	}
	if (rn) {
		return lhs->mul(num / rn->num);
	}

	t_long_frac n =
	     num * lhs->num / rhs->num;
	lhs->num = 1;
	rhs->num = 1;
	return (lhs / rhs)->mul(n);
}

t_func_tree::h_item t_func_tree::t_item_pow::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	t_item_num *ln = dynamic_cast<t_item_num *> (lhs.get());
	t_item_num *rn = dynamic_cast<t_item_num *> (rhs.get());
	if (ln && rn) {
		//NOTE: We use constraint on a maximum degree!
		if (rn->num.isint() && (abs(rn->num) < 100)) {
			return own.gener(
				pow(ln->num, rn->num.upper())
			)->mul(num);
		}
	}
	if (rn && (rn->num == 0)) { return own.gener(1); }
	if (ln && (ln->num == 1)) { return own.gener(1); }
	if (rn && (rn->num == 1)) return lhs->mul(num);

	return (lhs ^ rhs)->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_exp::red() const {

	h_item lhs = arg->red();
	if (lhs->num == 0) { return own.gener(num); }
	return
	own.exp(lhs)->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_log::red() const {

	h_item lhs = arg->red();
	if (dynamic_cast<t_item_num *> (lhs.get()) &&
	    lhs->num == 1) {
		return own.gener(0);
	}
	return
	own.log(lhs)->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_cos::red() const {

	h_item lhs = arg->red();
	if (lhs->num == 0) {
		return own.gener(num);
	}
	return
	own.cos(lhs)->mul(num);
}

t_func_tree::h_item t_func_tree::t_item_sin::red() const {

	h_item lhs = arg->red();
	if (lhs->num == 0) {
		return own.gener(0);
	}
	return
	own.sin(lhs)->mul(num);
}

//Генераторы узлов:

t_func_tree::h_item t_func_tree::store(t_item *&& _ptr) const { h_item hand(_ptr); LINK[_ptr] = hand; return hand; }

t_func_tree::h_item t_func_tree::gener(std::string str, const h_item &lhs, const h_item &rhs) const {

	if (str == "-") return store(new t_item_add(*this, lhs, rhs->mul(-1)));
	if (str == "+") return store(new t_item_add(*this, lhs, rhs));
	if (str == "*") return store(new t_item_mul(*this, lhs, rhs));
	if (str == "/") return store(new t_item_div(*this, lhs, rhs));
	if (str == "^") return store(new t_item_pow(*this, lhs, rhs));

	return nullptr;
}

t_func_tree::h_item t_func_tree::gener(std::string str, const h_item &rhs) const {

	if (str == "sqrt") return store(new t_item_pow(*this, rhs, gener(t_long_frac(1, 2))));
	if (str == "exp") return store(new t_item_exp(*this, rhs));
	if (str == "log") return store(new t_item_log(*this, rhs));
	if (str == "cos") return store(new t_item_cos(*this, rhs));
	if (str == "sin") return store(new t_item_sin(*this, rhs));
	if (str == "-u") return rhs->mul(-1);
	if (str == "-") return rhs->mul(-1);

	return nullptr;
}

t_func_tree::h_item t_func_tree::gener(std::string str) const {

	if ((str.length() == 1) && isalpha(str[0]) && islower(str[0])) {
		return store(new t_item_var(*this, str[0]));
	}
	return nullptr;
}

t_func_tree::h_item t_func_tree::gener(t_long_frac val) const {

	return store(new t_item_num(*this, val));
}

t_func_tree::h_item t_func_tree::gener(long val) const {

	return gener(t_long_frac(val));
}

//...

bool t_func_tree::create(std::string str) {

	static const std::regex expr("^\\s*(((\\d+)(\\.(\\d+))?)|([a-z]+)|(\\()|([\\+\\-\\/\\^]|\\*|\\)))\\s*");
	//...
	static const int IND_CONST_MAN = 5;
	static const int IND_CONST_INT = 3;
	static const int IND_CONST = 2;
	static const int IND_INDEX = 6;
	static const int IND_BRACE = 7;
	static const int IND_BSIGN = 8;
	//...
	static const int IS_START = 0;
	static const int IS_VALUE = 1;
	static const int IS_PFUNC = 2;
	static const int IS_BSIGN = 3;
	static const int IS_ERROR = 4;
	//...
	static std::map<std::string, int> LEVEL;
	static std::set<std::string> FUNC1;
	static std::set<std::string> SIGN1;
	static std::set<std::string> SIGN2;
	static bool start = true;
	if (start) {
		FUNC1.insert("sqrt");
		FUNC1.insert("log");
		FUNC1.insert("exp");
		FUNC1.insert("cos");
		FUNC1.insert("sin");
		//...
		SIGN1.insert("-u");
		//...
		SIGN2.insert("-");
		SIGN2.insert("+");
		SIGN2.insert("*");
		SIGN2.insert("/");
		SIGN2.insert("^");
		//...
		LEVEL["sqrt"] = LEVEL["log"] = LEVEL["exp"] = LEVEL["cos"] = LEVEL["sin"] = 6;
		LEVEL["-u"] = 6;
		LEVEL["^"] = 5;
		LEVEL["/"] = LEVEL["*"] = 4;
		LEVEL["+"] = LEVEL["-"] = 3;
		LEVEL[")"] = 2;
		LEVEL[","] = 1;
		LEVEL["("] = -1;
		LEVEL[""] = 0;
		//...
		start = false;
	}

	std::stack<h_item> OPER; std::stack<std::string> SIGN; std::string tmp; std::smatch res;
	int state = IS_START;
	root = nullptr;
	while (true) {
		if (str.size()) {
			if (!std::regex_search(str, res, expr)) {
				throw std::logic_error("There is incorrect token in the string!");
			}
			tmp = res[1].str();
		}
		else {
			tmp = "";
		}
		//Знаки арифметических операций или закрывающая скобка:
		if (res[IND_BSIGN].length() || (tmp == "")) {

			if (state != IS_VALUE && (state != IS_START || tmp != "-")) {
				throw std::logic_error("Incorrect location of the operations!");
			}
			int cur = LEVEL[tmp];
			while (!SIGN.empty() && (LEVEL[SIGN.top()] >= cur)) {
				//Для бинарных операторов:
				if (SIGN2.count(SIGN.top())) {
					h_item rhs = OPER.top();
					OPER.pop();
					h_item lhs = OPER.top();
					OPER.pop();
					OPER.push(gener(SIGN.top(), lhs, rhs));
				}
				//Для унарных операторов:
				if (SIGN1.count(SIGN.top())) {
					h_item rhs = OPER.top();
					OPER.pop();
					OPER.push(gener(SIGN.top(), rhs));
				}
				SIGN.pop();
			}
			if (tmp == "-" && state == IS_START) {
				SIGN.push("-u");
			}
			else
			if (tmp == ")") {
				state = IS_VALUE;
				if (!SIGN.size()) {
					throw std::logic_error("There is extra closing bracket!");
				}
				SIGN.pop();
			}
			else
			if (tmp == "") {
				if (SIGN.size()) {
					throw std::logic_error("There is extra opening bracket!");
				}
				break;
			}
			else {
				state = IS_BSIGN;
				SIGN.push(tmp);
			}
		}
		//Имена переменных и функции:
		else if (res[IND_INDEX].length()) {

			if (state == IS_VALUE) {
				throw std::logic_error("Incorrect location of the operands!");
			}
			if (!FUNC1.count(tmp) && !SIGN2.count(tmp)) {
				OPER.push(gener(tmp));
				if (!OPER.top()) {
				throw std::logic_error("Incorrect name of the operand!");
				}
				state = IS_VALUE;
			}
			else {
				state = IS_PFUNC;
				SIGN.push(tmp);
			}
		}
		//Открывающая скобка:
		else if (res[IND_BRACE].length()) {
			if (state != IS_START &&
			    state != IS_PFUNC &&
			    state != IS_BSIGN) {
				throw std::logic_error("Incorrect location of bracket!");
			}
			state = IS_START;
			SIGN.push(tmp);
		}
		//Числовые константы:
		else if (res[IND_CONST].length()) {
			t_long_frac val;
			std::string a = res[IND_CONST_INT].str();
			std::string b = res[IND_CONST_MAN].str();
			if (b.size()) {
				std::string s = "1";
				for (int i = 0; i < b.size(); ++ i) {
					s += "0";
				}
				val = t_long_frac(a + b + "/" + s);
			}
			else {
				val = t_long_frac(a);
			}
			OPER.push(gener(val));
			state = IS_VALUE;
		}
		str = res.suffix().str();
	}
	root = OPER.top();

	return true;
}

void t_func_tree::derive(char var) {
	root = root->dif(var);
}

void t_func_tree::reduce() {
	root = root->red();
}

//...

t_func_tree &t_func_tree::operator=(const t_func_tree &rhs) {

	std::copy(rhs.DATA, rhs.DATA + 26, DATA);
	root = rhs.root->cpy(*this);
	return *this;
}

t_func_tree::t_func_tree(const t_func_tree &src) {
	*this = src;
}

t_func_tree::t_func_tree(std::string str) {

	std::fill(DATA, DATA + 26, 0);
	create(str);
}

t_func_tree::t_func_tree() {

	std::fill(DATA, DATA + 26, 0);
}

//...

t_func_tree::t_item::~t_item() {
	own.LINK.erase(this);
}

//...
