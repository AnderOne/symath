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

//Методы дифференцирования узлов:

t_tree::h_node t_tree::t_node_pow::dif(char var) const {
	h_node hand = arg[1]->cpy() * arg[0]->dif(var) * (arg[0]->cpy() ^ (arg[1]->cpy() - own.gener(1))) +
	              arg[1]->dif(var) * own.log(arg[0]->cpy()) * cpy();
	return hand->mul(num);
}

t_tree::h_node t_tree::t_node_div::dif(char var) const {
	h_node hand = arg[0]->dif(var) / arg[1]->cpy() -
	              arg[0]->cpy() * (arg[1]->dif(var) / (arg[1]->cpy() ^ own.gener(2)));
	return hand->mul(num);
}

t_tree::h_node t_tree::t_node_mul::dif(char var) const {
	h_node hand = arg[0]->dif(var) * arg[1]->cpy() + arg[0]->cpy() * arg[1]->dif(var);
	return hand->mul(num);
}

t_tree::h_node t_tree::t_node_add::dif(char var) const {
	h_node hand = arg[0]->dif(var) + arg[1]->dif(var);
	return hand->mul(num);
}

t_tree::h_node t_tree::t_node_cos::dif(char var) const {
	h_node hand = - (arg->dif(var) * own.sin(arg->cpy()));
	return hand->mul(num);
}

t_tree::h_node t_tree::t_node_sin::dif(char var) const {
	h_node hand = arg->dif(var) * own.cos(arg->cpy());
	return hand->mul(num);
}

t_tree::h_node t_tree::t_node_log::dif(char var) const {
	h_node hand = arg->dif(var) / arg->cpy();
	return hand->mul(num);
}

t_tree::h_node t_tree::t_node_exp::dif(char var) const {
	h_node hand = arg->dif(var) * cpy();
	return hand->mul(num);
}

t_tree::h_node t_tree::t_node_var::dif(char var) const {
	h_node hand = own.gener(ind == var? num: 0);
	return hand;
}

t_tree::h_node t_tree::t_node_num::dif(char var) const {
	h_node hand = own.gener(0);
	return hand;
}

//...

}
