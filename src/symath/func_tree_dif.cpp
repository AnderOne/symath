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

#include <symath/func_tree.hpp>

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
	h_item hand = own.gener((ind == var)? num: t_frac(0));
	return hand;
}

t_func_tree::h_item t_func_tree::t_item_num::dif(char var) const {
	h_item hand = own.gener(0);
	return hand;
}

//...
