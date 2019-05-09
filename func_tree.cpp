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

t_func_tree::h_item t_func_tree::t_item_index::cpy(const t_func_tree &_own) const { return _own.gener(std::string(1, ind)); }

t_func_tree::h_item t_func_tree::t_item_const::cpy(const t_func_tree &_own) const { return _own.gener(val); }

t_func_tree::h_item t_func_tree::t_item_index::cpy() const { return own.gener(std::string(1, ind)); }

t_func_tree::h_item t_func_tree::t_item_const::cpy() const { return own.gener(val); }

#define __DEF_METH_CPY(func, dim)\
t_func_tree::h_item t_func_tree::t_item_##func::cpy(const t_func_tree &_own) const {\
	std::shared_ptr<t_item_##func> h(new t_item_##func(_own));\
	for (int i = 0; i < dim; ++ i) {\
		h->arg[i] = this->arg[i]->cpy(_own);\
	}\
	return h;\
}\
t_func_tree::h_item t_func_tree::t_item_##func::cpy() const {\
	std::shared_ptr<t_item_##func> h(new t_item_##func(*this));\
	for (int i = 0; i < dim; ++ i) {\
		h->arg[i] = this->arg[i]->cpy();\
	}\
	return h;\
}

__DEF_METH_CPY(sqrt, 1)
__DEF_METH_CPY(log, 1)
__DEF_METH_CPY(exp, 1)
__DEF_METH_CPY(cos, 1)
__DEF_METH_CPY(sin, 1)
__DEF_METH_CPY(neg, 1)
__DEF_METH_CPY(sub, 2)
__DEF_METH_CPY(add, 2)
__DEF_METH_CPY(mul, 2)
__DEF_METH_CPY(div, 2)
__DEF_METH_CPY(pow, 2)

//Перевод в строку:

std::string t_func_tree::t_item_index::str() const { return std::string(1, ind); }
std::string t_func_tree::t_item_const::str() const { return val; }

std::string t_func_tree::t_item_pow::str() const { return "(" + arg[0]->str() + ") ^ (" + arg[1]->str() + ")"; }
std::string t_func_tree::t_item_div::str() const { return "(" + arg[0]->str() + ") / (" + arg[1]->str() + ")"; }
std::string t_func_tree::t_item_mul::str() const { return "(" + arg[0]->str() + ") * (" + arg[1]->str() + ")"; }
std::string t_func_tree::t_item_add::str() const { return "(" + arg[0]->str() + ") + (" + arg[1]->str() + ")"; }
std::string t_func_tree::t_item_sub::str() const { return "(" + arg[0]->str() + ") - (" + arg[1]->str() + ")"; }

#define __DEF_METH_STR(func)\
std::string t_func_tree::t_item_##func::str() const { return #func"(" + arg[0]->str() + ")"; }

__DEF_METH_STR(sqrt)
__DEF_METH_STR(log)
__DEF_METH_STR(exp)
__DEF_METH_STR(cos)
__DEF_METH_STR(sin)

std::string t_func_tree::t_item_neg::str() const { return "-(" + arg[0]->str() + ")"; }

//Вычисление выражения:

double t_func_tree::t_item_index::get() const { return own[ind]; }
double t_func_tree::t_item_const::get() const { return val; }

double t_func_tree::t_item_pow::get() const { return pow(arg[0]->get(), arg[1]->get()); }
double t_func_tree::t_item_div::get() const { return arg[0]->get() / arg[1]->get(); }
double t_func_tree::t_item_mul::get() const { return arg[0]->get() * arg[1]->get(); }
double t_func_tree::t_item_add::get() const { return arg[0]->get() + arg[1]->get(); }
double t_func_tree::t_item_sub::get() const { return arg[0]->get() - arg[1]->get(); }

#define __DEF_METH_GET(func)\
double t_func_tree::t_item_##func::get() const { return std::func(arg[0]->get()); }

__DEF_METH_GET(sqrt)
__DEF_METH_GET(log)
__DEF_METH_GET(exp)
__DEF_METH_GET(cos)
__DEF_METH_GET(sin)

double t_func_tree::t_item_neg::get() const { return - arg[0]->get(); }

//Методы дифференцирования узлов:

t_func_tree::h_item t_func_tree::t_item_index::dif(char var) const {
	return own.gener((var == ind)? 1: 0);
}

t_func_tree::h_item t_func_tree::t_item_const::dif(char var) const {
	return own.gener(0);
}

t_func_tree::h_item t_func_tree::t_item_sqrt::dif(char var) const {
	return arg[0]->dif(var) / (own.gener(2) * cpy());
}

t_func_tree::h_item t_func_tree::t_item_log::dif(char var) const {
	return arg[0]->dif(var) / arg[0]->cpy();
}

t_func_tree::h_item t_func_tree::t_item_exp::dif(char var) const {
	return arg[0]->dif(var) * cpy();
}

t_func_tree::h_item t_func_tree::t_item_cos::dif(char var) const {
	return - (arg[0]->dif(var) * own.sin(arg[0]->cpy()));
}

t_func_tree::h_item t_func_tree::t_item_sin::dif(char var) const {
	return (arg[0]->dif(var) * own.cos(arg[0]->cpy()));
}

t_func_tree::h_item t_func_tree::t_item_neg::dif(char var) const {
	return - arg[0]->dif(var);
}

t_func_tree::h_item t_func_tree::t_item_sub::dif(char var) const {
	return arg[0]->dif(var) - arg[1]->dif(var);
}

t_func_tree::h_item t_func_tree::t_item_add::dif(char var) const {
	return arg[0]->dif(var) + arg[1]->dif(var);
}

t_func_tree::h_item t_func_tree::t_item_mul::dif(char var) const {
	return arg[0]->dif(var) * arg[1]->cpy() +
	       arg[0]->cpy() * arg[1]->dif(var);
}

t_func_tree::h_item t_func_tree::t_item_div::dif(char var) const {
	return arg[0]->dif(var) / arg[1]->cpy() -
	       arg[0]->cpy() * arg[1]->dif(var) / arg[1]->cpy();
}

t_func_tree::h_item t_func_tree::t_item_pow::dif(char var) const {
	return arg[1]->cpy() * arg[0]->dif(var) * (arg[0]->cpy() ^ (arg[1]->cpy() - own.gener(1))) +
	       arg[1]->dif(var) * own.log(arg[0]->cpy()) * cpy();
}

//Методы сокращения выражений:

t_func_tree::h_item t_func_tree::t_item_index::red() const {
	return cpy();
}

t_func_tree::h_item t_func_tree::t_item_const::red() const {
	return cpy();
}

#define __DEF_METH_RED(func)\
t_func_tree::h_item t_func_tree::t_item_##func::red() const {\
	h_item lhs = arg[0]->red();\
	return own.func(lhs);\
}

__DEF_METH_RED(sqrt)
__DEF_METH_RED(log)
__DEF_METH_RED(exp)
__DEF_METH_RED(cos)
__DEF_METH_RED(sin)

t_func_tree::h_item t_func_tree::t_item_neg::red() const {

	h_item lhs = arg[0]->red();

	const t_item_const *lc = dynamic_cast<const t_item_const *> (lhs.get());
	if (lc != nullptr) {
		return own.gener(- lc->val);
	}

	const t_item_neg *ln = dynamic_cast<const t_item_neg *> (lhs.get());
	if (ln != nullptr) {
		return ln->arg[0];
	}

	return - lhs;
}

t_func_tree::h_item t_func_tree::t_item_sub::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	const t_item_const *lc = dynamic_cast<const t_item_const *> (lhs.get());
	const t_item_const *rc = dynamic_cast<const t_item_const *> (rhs.get());
	if (lc && rc) {
		return own.gener(lc->val - rc->val);
	}
	if (lc && (lc->val == 0)) {
		return (- rhs)->red();
	}
	if (rc && (rc->val == 0)) {
		return lhs;
	}

	const t_item_neg *rn = dynamic_cast<const t_item_neg *> (rhs.get());
	if (rn) {
		return lhs + rn->arg[0];
	}

	return lhs - rhs;
}

t_func_tree::h_item t_func_tree::t_item_add::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	const t_item_const *lc = dynamic_cast<const t_item_const *> (lhs.get());
	const t_item_const *rc = dynamic_cast<const t_item_const *> (rhs.get());
	if (lc && rc) {
		return own.gener(lc->val + rc->val);
	}
	if (lc && (lc->val == 0)) {
		return rhs;
	}
	if (rc && (rc->val == 0)) {
		return lhs;
	}

	const t_item_neg *rn = dynamic_cast<const t_item_neg *> (rhs.get());
	if (rn) {
		return (lhs - rn->arg[0])->red();
	}

	return lhs + rhs;
}

t_func_tree::h_item t_func_tree::t_item_mul::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	const t_item_const *lc = dynamic_cast<const t_item_const *> (lhs.get());
	const t_item_const *rc = dynamic_cast<const t_item_const *> (rhs.get());
	if (lc && rc) {
		return own.gener(lc->val * rc->val);
	}
	if (lc && (lc->val == 0)) {
		return own.gener(0);
	}
	if (rc && (rc->val == 0)) {
		return own.gener(0);
	}
	if (lc && (lc->val == -1)) {
		return - rhs;
	}
	if (rc && (rc->val == -1)) {
		return - lhs;
	}
	if (lc && (lc->val == 1)) {
		return rhs;
	}
	if (rc && (rc->val == 1)) {
		return lhs;
	}

	const t_item_neg *ln = dynamic_cast<const t_item_neg *> (lhs.get());
	const t_item_neg *rn = dynamic_cast<const t_item_neg *> (rhs.get());
	if (ln && rn) {
		return (ln->arg[0] * rn->arg[0])->red();
	}
	if (ln) {
		return (- (ln->arg[0] * rhs))->red();
	}
	if (rn) {
		return (- (lhs * rn->arg[0]))->red();
	}

	return lhs * rhs;
}

t_func_tree::h_item t_func_tree::t_item_div::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	const t_item_const *lc = dynamic_cast<const t_item_const *> (lhs.get());
	const t_item_const *rc = dynamic_cast<const t_item_const *> (rhs.get());
	if (lc && rc) {
		return own.gener(lc->val / rc->val);
	}
	if (lc && (lc->val == 0)) {
		return own.gener(0);
	}
	if (rc && (rc->val == - 1)) {
		return - lhs;
	}
	if (rc && (rc->val == 1)) {
		return lhs;
	}

	const t_item_neg *ln = dynamic_cast<const t_item_neg *> (lhs.get());
	const t_item_neg *rn = dynamic_cast<const t_item_neg *> (rhs.get());
	if (ln && rn) {
		return (ln->arg[0] / rn->arg[0])->red();
	}
	if (ln) {
		return (- (ln->arg[0] / rhs))->red();
	}
	if (rn) {
		return (- (lhs / rn->arg[0]))->red();
	}

	return lhs / rhs;
}

t_func_tree::h_item t_func_tree::t_item_pow::red() const {

	h_item lhs = arg[0]->red();
	h_item rhs = arg[1]->red();

	const t_item_const *lc = dynamic_cast<const t_item_const *> (lhs.get());
	const t_item_const *rc = dynamic_cast<const t_item_const *> (rhs.get());
	//TODO: ...
	/*if (lc && rc) {
		return own.gener(this->get());
	}*/
	if (lc && (lc->val == 1)) {
		return own.gener(1);
	}
	if (rc && (rc->val == 1)) {
		return lhs;
	}

	return lhs ^ rhs;
}

//Генераторы узлов:

t_func_tree::h_item t_func_tree::store(t_item *&& _ptr) const { h_item hand(_ptr); LINK[_ptr] = hand; return hand; }

t_func_tree::h_item t_func_tree::gener(std::string str, const h_item &lhs, const h_item &rhs) const {

	if (str == "^") return store(new t_item_pow(*this, lhs, rhs));
	if (str == "/") return store(new t_item_div(*this, lhs, rhs));
	if (str == "*") return store(new t_item_mul(*this, lhs, rhs));
	if (str == "+") return store(new t_item_add(*this, lhs, rhs));
	if (str == "-") return store(new t_item_sub(*this, lhs, rhs));

	return nullptr;
}

t_func_tree::h_item t_func_tree::gener(std::string str, const h_item &rhs) const {

	if (str == "sqrt") return store(new t_item_sqrt(*this, rhs));
	if (str == "log") return store(new t_item_log(*this, rhs));
	if (str == "exp") return store(new t_item_exp(*this, rhs));
	if (str == "cos") return store(new t_item_cos(*this, rhs));
	if (str == "sin") return store(new t_item_sin(*this, rhs));
	if (str == "-") return store(new t_item_neg(*this, rhs));

	return nullptr;
}

t_func_tree::h_item t_func_tree::gener(std::string str) const {

	if ((str.length() == 1) && isalpha(str[0]) && islower(str[0])) {
		return store(new t_item_index(*this, str[0]));
	}
	return nullptr;
}

t_func_tree::h_item t_func_tree::gener(t_long_frac val) const {

	return store(new t_item_const(*this, val));
}

t_func_tree::h_item t_func_tree::gener(long val) const {

	return store(
	new t_item_const(*this, t_long_frac(val))
	);
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
