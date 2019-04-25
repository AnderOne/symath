#include "func_tree.hpp"

#include <type_traits>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <regex>
#include <algorithm>
#include <vector>
#include <string>
#include <stack>
#include <map>
#include <set>

//Методы клонирования узлов:

t_func_tree::h_item t_func_tree::t_item_const::cpy(const t_func_tree &_own) const { return h_item(new t_item_const(_own, val)); }

t_func_tree::h_item t_func_tree::t_item_index::cpy(const t_func_tree &_own) const { return h_item(new t_item_index(_own, ind)); }

t_func_tree::h_item t_func_tree::t_item_const::cpy() const { return h_item(new t_item_const(*this)); }

t_func_tree::h_item t_func_tree::t_item_index::cpy() const { return h_item(new t_item_index(*this)); }

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

std::string t_func_tree::t_item_const::str() const { return std::to_string(val); }
std::string t_func_tree::t_item_index::str() const { return std::string() + ind; }

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

double t_func_tree::t_item_const::get() const { return (double) val; }
double t_func_tree::t_item_index::get() const { return own[ind]; }

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

//Генераторы узлов:

t_func_tree::h_item t_func_tree::gener(std::string str, const h_item &lhs, const h_item &rhs) const {

	if (str == "^") return h_item(new t_item_pow(*this, lhs, rhs));
	if (str == "/") return h_item(new t_item_div(*this, lhs, rhs));
	if (str == "*") return h_item(new t_item_mul(*this, lhs, rhs));
	if (str == "+") return h_item(new t_item_add(*this, lhs, rhs));
	if (str == "-") return h_item(new t_item_sub(*this, lhs, rhs));

	return nullptr;
}

t_func_tree::h_item t_func_tree::gener(std::string str, const h_item &rhs) const {

	if (str == "sqrt") return h_item(new t_item_sqrt(*this, rhs));
	if (str == "log") return h_item(new t_item_log(*this, rhs));
	if (str == "exp") return h_item(new t_item_exp(*this, rhs));
	if (str == "cos") return h_item(new t_item_cos(*this, rhs));
	if (str == "sin") return h_item(new t_item_sin(*this, rhs));
	if (str == "-") return h_item(new t_item_neg(*this, rhs));

	return nullptr;
}

t_func_tree::h_item t_func_tree::gener(std::string str) const {

	return h_item(new t_item_index(*this, str[0]));
}

t_func_tree::h_item t_func_tree::gener(double val) const {

	return h_item(new t_item_const(*this, val));
}

//...

bool t_func_tree::create(std::string str) {

	static const std::regex expr("\\s*((\\d+(\\.\\d+)?)|([a-z]+)|(\\()|([\\+\\-\\/\\^]|\\*|\\)|,))");
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
		LEVEL["("] = 0;
		LEVEL[""] = -1;
		//...
		start = false;
	}
	std::stack<h_item> OPER; std::stack<std::string> SIGN; std::smatch res; bool first = true;
	while (true) {
		std::string tmp = ""; if (std::regex_search(str, res, expr)) tmp = res[1].str();
		//Бинарные операции, запятая или закрывающая скобка:
		if (res[6].length() || (tmp == "")) {
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
			if (tmp != "," && tmp != ")") SIGN.push(tmp);
			if (tmp == "-" && first) OPER.push(gener(0));
			if (tmp == ")") SIGN.pop();
			first = (tmp == ",");
			if (tmp == "") break;
		}
		//Имена переменных и функции:
		else if (res[4].length()) {
			if (!FUNC1.count(tmp) && !FUNC2.count(tmp)) {
				OPER.push(gener(tmp));
				first = false;
			}
			else {
				SIGN.push(tmp);
				first = true;
			}
		}
		//Открывающая скобка:
		else if (res[5].length()) {
			SIGN.push(tmp);
			first = true;
		}
		//Числовые константы:
		else if (res[2].length()) {
			OPER.push(
			gener(
			atof(res[2].str().c_str())
			)
			);
			first = false;
		}
		//...
		if (tmp != "") {
			str = res.suffix().str();
			continue;
		}
		break;
	}
	root = OPER.top();

	return true;
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
