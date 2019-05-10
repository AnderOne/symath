#include "func_tree.hpp"

#include <type_traits>
#include <exception>
#include <cassert>
#include <cerrno>
#include <climits>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <regex>
#include <algorithm>
#include <vector>
#include <string>
#include <stack>
#include <map>
#include <set>

//Методы клонирования узлов:

t_func_tree::h_item t_func_tree::t_item_var::cpy(const t_func_tree &_own) const { return _own.gener(std::string(1, ind)); }
t_func_tree::h_item t_func_tree::t_item_var::cpy() const { return own.gener(std::string(1, ind)); }

t_func_tree::h_item t_func_tree::t_item_num::cpy(const t_func_tree &_own) const { return _own.gener(val); }
t_func_tree::h_item t_func_tree::t_item_num::cpy() const { return own.gener(val); }

#define __DEF_ITEM_BIN_CPY(ind, op) \
t_func_tree::h_item t_func_tree::t_item_##ind::cpy(const t_func_tree &_own) const { return _own.gener(#op, arg[0], arg[1]); }\
t_func_tree::h_item t_func_tree::t_item_##ind::cpy() const { return own.gener(#op, arg[0], arg[1]); }

__DEF_ITEM_BIN_CPY(sub, -)
__DEF_ITEM_BIN_CPY(add, +)
__DEF_ITEM_BIN_CPY(mul, *)
__DEF_ITEM_BIN_CPY(div, /)
__DEF_ITEM_BIN_CPY(pow, ^)

#undef __DEF_ITEM_BIN_CPY

#define __DEF_ITEM_ONE_CPY(ind) \
t_func_tree::h_item t_func_tree::t_item_##ind::cpy(const t_func_tree &_own) const { return _own.gener(#ind, arg); }\
t_func_tree::h_item t_func_tree::t_item_##ind::cpy() const { return own.gener(#ind, arg); }

__DEF_ITEM_ONE_CPY(exp)
__DEF_ITEM_ONE_CPY(log)
__DEF_ITEM_ONE_CPY(cos)
__DEF_ITEM_ONE_CPY(sin)

#undef __DEF_ITEM_ONE_CPY

//Перевод в строку:

std::string t_func_tree::t_item_var::str() const { return std::string(1, ind); }
std::string t_func_tree::t_item_num::str() const { return val; }

#define __DEF_ITEM_BIN_STR(ind, op) \
std::string t_func_tree::t_item_##ind::str() const {\
	return "(" + arg[0]->str() + ") "#op" (" + arg[1]->str() + ")";\
}

__DEF_ITEM_BIN_STR(sub, -)
__DEF_ITEM_BIN_STR(add, +)
__DEF_ITEM_BIN_STR(mul, *)
__DEF_ITEM_BIN_STR(div, /)
__DEF_ITEM_BIN_STR(pow, ^)

#undef __DEF_ITEM_BIN_STR

#define __DEF_ITEM_ONE_STR(ind) \
std::string t_func_tree::t_item_##ind::str() const {\
	return #ind"(" + arg->str() + ")";\
}

__DEF_ITEM_ONE_STR(exp)
__DEF_ITEM_ONE_STR(log)
__DEF_ITEM_ONE_STR(cos)
__DEF_ITEM_ONE_STR(sin)

#undef __DEF_ITEM_ONE_STR

//Вычисление выражения:

double t_func_tree::t_item_var::get() const { return own[ind]; }
double t_func_tree::t_item_num::get() const { return val; }

double t_func_tree::t_item_pow::get() const {
	return std::pow(arg[0]->get(), arg[1]->get());
}
#define __DEF_ITEM_BIN_GET(ind, op) \
double t_func_tree::t_item_##ind::get() const {\
	return arg[0]->get() op arg[1]->get();\
}

__DEF_ITEM_BIN_GET(sub, -)
__DEF_ITEM_BIN_GET(add, +)
__DEF_ITEM_BIN_GET(mul, *)
__DEF_ITEM_BIN_GET(div, /)

#undef __DEF_ITEM_BIN_GET

#define __DEF_ITEM_ONE_GET(ind) \
double t_func_tree::t_item_##ind::get() const {\
	return std::ind(arg->get());\
}

__DEF_ITEM_ONE_GET(exp)
__DEF_ITEM_ONE_GET(log)
__DEF_ITEM_ONE_GET(cos)
__DEF_ITEM_ONE_GET(sin)

#undef __DEF_ITEM_ONE_GET

//Методы дифференцирования узлов:

t_func_tree::h_item t_func_tree::t_item_pow::dif(char var) const {
	return arg[1]->cpy() * arg[0]->dif(var) * (arg[0]->cpy() ^ (arg[1]->cpy() - own.gener(1))) +
	       arg[1]->dif(var) * own.log(arg[0]->cpy()) * cpy();
}

t_func_tree::h_item t_func_tree::t_item_div::dif(char var) const {
	return arg[0]->dif(var) / arg[1]->cpy() - arg[0]->cpy() * arg[1]->dif(var) / arg[1]->cpy();
}

t_func_tree::h_item t_func_tree::t_item_mul::dif(char var) const {
	return arg[0]->dif(var) * arg[1]->cpy() + arg[0]->cpy() * arg[1]->dif(var);
}

t_func_tree::h_item t_func_tree::t_item_add::dif(char var) const {
	return arg[0]->dif(var) + arg[1]->dif(var);
}

t_func_tree::h_item t_func_tree::t_item_sub::dif(char var) const {
	return arg[0]->dif(var) - arg[1]->dif(var);
}

t_func_tree::h_item t_func_tree::t_item_cos::dif(char var) const {
	return - (arg->dif(var) * own.sin(arg->cpy()));
}

t_func_tree::h_item t_func_tree::t_item_sin::dif(char var) const {
	return arg->dif(var) * own.cos(arg->cpy());
}

t_func_tree::h_item t_func_tree::t_item_log::dif(char var) const {
	return arg->dif(var) / arg->cpy();
}

t_func_tree::h_item t_func_tree::t_item_exp::dif(char var) const {
	return arg->dif(var) * cpy();
}

t_func_tree::h_item t_func_tree::t_item_var::dif(char var) const {
	return own.gener((ind == var)? 1: 0);
}

t_func_tree::h_item t_func_tree::t_item_num::dif(char var) const {
	return own.gener(0);
}

//Методы сокращения выражений:

t_func_tree::h_item t_func_tree::t_item_sub::split(t_long_frac &_num) const {

	t_item_num *num = dynamic_cast<t_item_num *> (arg[0].get());
	if (num) {
		_num =   num->val; return own.gener("-", arg[1]);
	}
	num = dynamic_cast<t_item_num *> (arg[1].get());
	if (num) {
		_num = - num->val;
		return arg[0];
	}
	return nullptr;
}

t_func_tree::h_item t_func_tree::t_item_add::split(t_long_frac &_num) const {

	t_item_num *num = dynamic_cast<t_item_num *> (arg[0].get());
	if (num) {
		_num = num->val;
		return arg[1];
	}
	num = dynamic_cast<t_item_num *> (arg[1].get());
	if (num) {
		_num = num->val;
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

t_func_tree::h_item t_func_tree::t_item_sub::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	t_item_num *ln = dynamic_cast<t_item_num *> (lhs.get());
	t_item_num *rn = dynamic_cast<t_item_num *> (rhs.get());
	if (ln && rn) {
		return own.gener(ln->val - rn->val);
	}
	if (rn && (rn->val == 0)) {
		return lhs;
	}

	t_item_bin *lb = dynamic_cast<t_item_bin *> (lhs.get());
	t_item_bin *rb = dynamic_cast<t_item_bin *> (rhs.get());
	t_long_frac lv, rv;
	h_item l2, r2;
	if (lb) l2 = lb->split(lv);
	if (rb) r2 = rb->split(rv);
	if (l2 && r2) {
		return own.gener(lv - rv) + (l2 - r2);
	}
	if (l2 && rn) {
		return own.gener(lv - rn->val) + l2;
	}
	if (ln && r2) {
		return own.gener(ln->val - rv) - r2;
	}

	return lhs - rhs;
}

t_func_tree::h_item t_func_tree::t_item_add::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	t_item_num *ln = dynamic_cast<t_item_num *> (lhs.get());
	t_item_num *rn = dynamic_cast<t_item_num *> (rhs.get());
	if (ln && rn) {
		return own.gener(ln->val + rn->val);
	}
	if (ln && (ln->val == 0)) {
		return rhs;
	}
	if (rn && (rn->val == 0)) {
		return lhs;
	}

	t_item_bin *lb = dynamic_cast<t_item_bin *> (lhs.get());
	t_item_bin *rb = dynamic_cast<t_item_bin *> (rhs.get());
	t_long_frac lv, rv;
	h_item l2, r2;
	if (lb) l2 = lb->split(lv);
	if (rb) r2 = rb->split(rv);
	if (l2 && r2) {
		return own.gener(lv + rv) + (l2 + r2);
	}
	if (l2 && rn) {
		return own.gener(lv + rn->val) + l2;
	}
	if (ln && r2) {
		return own.gener(ln->val + rv) + r2;
	}

	return lhs + rhs;
}

t_func_tree::h_item t_func_tree::t_item_mul::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	t_item_num *ln = dynamic_cast<t_item_num *> (lhs.get());
	t_item_num *rn = dynamic_cast<t_item_num *> (rhs.get());
	if (ln && rn) {
		return own.gener(ln->val * rn->val);
	}
	if (ln && (ln->val == 0)) {
		return own.gener(0);
	}
	if (rn && (rn->val == 0)) {
		return own.gener(0);
	}
	if (ln && (ln->val == 1)) {
		return rhs;
	}
	if (rn && (rn->val == 1)) {
		return lhs;
	}

	t_item_bin *lb = dynamic_cast<t_item_bin *> (lhs.get());
	t_item_bin *rb = dynamic_cast<t_item_bin *> (rhs.get());
	t_long_frac lv, rv;
	h_item l2, r2;
	if (lb) l2 = lb->split(lv);
	if (rb) r2 = rb->split(rv);
	if (l2 && rn) {
		return own.gener(lv * rn->val) + l2 * rn->val;
	}
	if (ln && r2) {
		return own.gener(ln->val * rv) + r2 * ln->val;
	}

	return lhs * rhs;
}

t_func_tree::h_item t_func_tree::t_item_div::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	t_item_num *ln = dynamic_cast<t_item_num *> (lhs.get());
	t_item_num *rn = dynamic_cast<t_item_num *> (rhs.get());
	if (ln && rn) {
		return own.gener((ln->val) /
		                 (rn->val));
	}
	if (ln && (ln->val == 0)) {
		return own.gener(0);
	}
	if (rn && (rn->val == 1)) {
		return lhs;
	}

	t_item_bin *lb = dynamic_cast<t_item_bin *> (lhs.get());
	t_long_frac lv;
	h_item l2;
	if (lb) l2 = lb->split(lv);
	if (l2 && rn) {
		return own.gener(lv / rn->val) + l2 / rn->val;
	}

	return lhs / rhs;
}

t_func_tree::h_item t_func_tree::t_item_pow::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	t_item_num *ln = dynamic_cast<t_item_num *> (lhs.get());
	t_item_num *rn = dynamic_cast<t_item_num *> (rhs.get());
	//TODO: ...
	/*if (ln && rn) {
		//...
	}*/
	if (rn && (rn->val == 0)) {
		return own.gener(1);
	}
	if (ln && (ln->val == 1)) {
		return own.gener(1);
	}
	if (rn && (rn->val == 1)) {
		return lhs;
	}

	return lhs ^ rhs;
}

t_func_tree::h_item t_func_tree::t_item_exp::red() const {

	h_item lhs = arg->red();
	//TODO: ...
	return own.exp(lhs);
}

t_func_tree::h_item t_func_tree::t_item_log::red() const {

	h_item lhs = arg->red();
	//TODO: ...
	return own.log(lhs);
}

t_func_tree::h_item t_func_tree::t_item_cos::red() const {

	h_item lhs = arg->red();
	//TODO: ...
	return own.cos(lhs);
}

t_func_tree::h_item t_func_tree::t_item_sin::red() const {

	h_item lhs = arg->red();
	//TODO: ...
	return own.sin(lhs);
}

//Генераторы узлов:

t_func_tree::h_item t_func_tree::store(t_item *&& _ptr) const { h_item hand(_ptr); LINK[_ptr] = hand; return hand; }

t_func_tree::h_item t_func_tree::gener(std::string str, const h_item &lhs, const h_item &rhs) const {

	if (str == "-") return store(new t_item_sub(*this, lhs, rhs));
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
	if (str == "-") return gener("-", gener(0), rhs);

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
	static const int IND_CONST = 2;
	static const int IND_CONST_INT = 3;
	static const int IND_CONST_MAN = 5;
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
	static std::set<std::string> FUNC2;
	static bool start = true;
	if (start) {
		FUNC1.insert("sqrt");
		FUNC1.insert("log");
		FUNC1.insert("exp");
		FUNC1.insert("cos");
		FUNC1.insert("sin");
		//...
		FUNC2.insert("^");
		FUNC2.insert("/");
		FUNC2.insert("*");
		FUNC2.insert("+");
		FUNC2.insert("-");
		//...
		LEVEL["sqrt"] = LEVEL["log"] = LEVEL["exp"] = LEVEL["cos"] = LEVEL["sin"] = 6;
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
				if (FUNC2.count(SIGN.top())) {
					h_item rhs = OPER.top();
					OPER.pop();
					h_item lhs = OPER.top();
					OPER.pop();
					OPER.push(gener(SIGN.top(), lhs, rhs));
				}
				//Для унарных операторов:
				if (FUNC1.count(SIGN.top())) {
					h_item rhs = OPER.top();
					OPER.pop();
					OPER.push(gener(SIGN.top(), rhs));
				}
				SIGN.pop();
			}
			if (tmp == "-" && state == IS_START) {
				OPER.push(gener(0));
			}
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
			if (!FUNC1.count(tmp) && !FUNC2.count(tmp)) {
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
