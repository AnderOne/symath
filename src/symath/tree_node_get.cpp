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

//Вычисление выражения:

double t_tree::t_node_var::get() const { return num.get() * own[ind]; }
double t_tree::t_node_num::get() const { return num.get(); }

double t_tree::t_node_pow::get() const {
	return num.get() * std::pow(arg[0]->get(), arg[1]->get());
}

#define __DEF_ITEM_BIN_GET(ind, op) \
double t_tree::t_node_##ind::get() const {\
	return num.get() * (arg[0]->get() op arg[1]->get());\
}

__DEF_ITEM_BIN_GET(add, +)
__DEF_ITEM_BIN_GET(mul, *)
__DEF_ITEM_BIN_GET(div, /)

#define __DEF_ITEM_ONE_GET(ind) \
double t_tree::t_node_##ind::get() const {\
	return num.get() * std::ind(arg->get());\
}

__DEF_ITEM_ONE_GET(exp)
__DEF_ITEM_ONE_GET(log)
__DEF_ITEM_ONE_GET(cos)
__DEF_ITEM_ONE_GET(sin)

//...

}
