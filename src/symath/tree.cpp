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

#include <symath/tree.hpp>

//Генераторы узлов:

t_tree::h_node t_tree::store(t_node *&& _ptr) const { h_node hand(_ptr); LINK[_ptr] = hand; return hand; }

t_tree::h_node t_tree::gener(std::string str, const h_node &lhs, const h_node &rhs) const {

	if (str == "-") return store(new t_node_add(*this, lhs, rhs->mul(-1)));
	if (str == "+") return store(new t_node_add(*this, lhs, rhs));
	if (str == "*") return store(new t_node_mul(*this, lhs, rhs));
	if (str == "/") return store(new t_node_div(*this, lhs, rhs));
	if (str == "^") return store(new t_node_pow(*this, lhs, rhs));

	return nullptr;
}

t_tree::h_node t_tree::gener(std::string str, const h_node &rhs) const {

	if (str == "sqrt") return store(new t_node_pow(*this, rhs, gener(t_frac(1, 2))));
	if (str == "exp") return store(new t_node_exp(*this, rhs));
	if (str == "log") return store(new t_node_log(*this, rhs));
	if (str == "cos") return store(new t_node_cos(*this, rhs));
	if (str == "sin") return store(new t_node_sin(*this, rhs));
	if (str == "-u") return rhs->mul(-1);
	if (str == "-") return rhs->mul(-1);

	return nullptr;
}

t_tree::h_node t_tree::gener(std::string str) const {

	if ((str.length() == 1) && isalpha(str[0]) && islower(str[0])) {
		return store(new t_node_var(*this, str[0]));
	}
	return nullptr;
}

t_tree::h_node t_tree::gener(t_frac val) const {

	return store(new t_node_num(*this, val));
}

t_tree::h_node t_tree::gener(long val) const {

	return gener(t_frac(val));
}

//...

bool t_tree::create(std::string str) {

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

	std::stack<h_node> OPER; std::stack<std::string> SIGN; std::string tmp; std::smatch res;
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
				//Для унарных операторов и функций:
				if (FUNC1.count(SIGN.top()) || SIGN1.count(SIGN.top())) {
					h_node rhs = OPER.top();
					OPER.pop();
					OPER.push(gener(SIGN.top(), rhs));
				}
				//Для бинарных операторов:
				if (SIGN2.count(SIGN.top())) {
					h_node rhs = OPER.top();
					OPER.pop();
					h_node lhs = OPER.top();
					OPER.pop();
					OPER.push(gener(SIGN.top(), lhs, rhs));
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
			t_frac val;
			std::string a = res[IND_CONST_INT].str();
			std::string b = res[IND_CONST_MAN].str();
			if (b.size()) {
				std::string s = "1";
				for (int i = 0; i < b.size(); ++ i) {
					s += "0";
				}
				val = t_frac(a + b + "/" + s);
			}
			else {
				val = t_frac(a);
			}
			OPER.push(gener(val));
			state = IS_VALUE;
		}
		str = res.suffix().str();
	}
	root = OPER.top();

	return true;
}

void t_tree::derive(char var) {
	root = root->dif(var);
}

void t_tree::reduce() {
	root = root->red();
}

//...

t_tree &t_tree::operator=(const t_tree &rhs) {

	std::copy(rhs.DATA, rhs.DATA + 26, DATA);
	root = rhs.root->cpy(*this);
	return *this;
}

t_tree::t_tree(const t_tree &src) {
	*this = src;
}

t_tree::t_tree(std::string str) {

	std::fill(DATA, DATA + 26, 0);
	create(str);
}

t_tree::t_tree() {

	std::fill(DATA, DATA + 26, 0);
}

//...

t_tree::t_node::~t_node() {
	own.LINK.erase(this);
}

//...
