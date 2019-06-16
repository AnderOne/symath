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

//Методы клонирования узлов:

t_tree::h_item t_tree::t_item_var::cpy(const t_tree &_own) const {
	return _own.gener(std::string(1, ind))->mul(num);
}
t_tree::h_item t_tree::t_item_var::cpy() const {
	return own.gener(std::string(1, ind))->mul(num);
}

t_tree::h_item t_tree::t_item_num::cpy(const t_tree &_own) const {
	return _own.gener(num);
}
t_tree::h_item t_tree::t_item_num::cpy() const {
	return own.gener(num);
}

#define __DEF_ITEM_BIN_CPY(ind, op) \
t_tree::h_item \
t_tree::t_item_##ind::cpy(const t_tree &_own) const {\
	return _own.gener(#op, arg[0], arg[1])->mul(num);\
}\
t_tree::h_item \
t_tree::t_item_##ind::cpy() const {\
	return own.gener(#op, arg[0], arg[1])->mul(num);\
}

__DEF_ITEM_BIN_CPY(add, +)
__DEF_ITEM_BIN_CPY(mul, *)
__DEF_ITEM_BIN_CPY(div, /)
__DEF_ITEM_BIN_CPY(pow, ^)

#define __DEF_ITEM_ONE_CPY(ind) \
t_tree::h_item \
t_tree::t_item_##ind::cpy(const t_tree &_own) const {\
	return _own.gener(#ind, arg)->mul(num);\
}\
t_tree::h_item \
t_tree::t_item_##ind::cpy() const {\
	return own.gener(#ind, arg)->mul(num);\
}

__DEF_ITEM_ONE_CPY(exp)
__DEF_ITEM_ONE_CPY(log)
__DEF_ITEM_ONE_CPY(cos)
__DEF_ITEM_ONE_CPY(sin)

//...
