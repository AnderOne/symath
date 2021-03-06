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

namespace symath {

//Перевод в строку:

std::string t_tree::t_node_var::str() const {
	return ((num < 0)? ("- "): ("")) + ((num.abs() != 1)? (num.abs().str() + " * "): ("")) + std::string(1, ind);
}

std::string t_tree::t_node_num::str() const {
	return ((num < 0)? ("- "): ("")) + num.abs().str();
}

std::string t_tree::t_node_add::str() const {
	std::string tmp;
	tmp = ((num < 0)? ("- "): ("")) + ((num.abs() != 1)? (num.abs().str() + " * "): (""));
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
std::string t_tree::t_node_##ind::str() const {\
	std::string tmp;\
	tmp = ((num < 0)? ("- "): ("")) + ((num.abs() != 1)? (num.abs().str() + " * "): (""));\
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

#define __DEF_ITEM_ONE_STR(ind) \
std::string t_tree::t_node_##ind::str() const {\
	std::string tmp;\
	tmp = ((num < 0)? ("- "): ("")) + ((num.abs() != 1)? (num.abs().str() + " * "): (""));\
	tmp += #ind"(";\
	tmp += arg->str();\
	tmp += ")";\
	return tmp;\
}

__DEF_ITEM_ONE_STR(exp)
__DEF_ITEM_ONE_STR(log)
__DEF_ITEM_ONE_STR(cos)
__DEF_ITEM_ONE_STR(sin)

//...

}
