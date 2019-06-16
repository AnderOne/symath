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

//Методы сокращения выражений:

t_func_tree::h_item t_func_tree::t_item_add::split(t_frac &_num) const {

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

t_func_tree::h_item t_func_tree::t_item_mul::split(t_frac &_num) const {
	return nullptr;
}

t_func_tree::h_item t_func_tree::t_item_div::split(t_frac &_num) const {
	return nullptr;
}

t_func_tree::h_item t_func_tree::t_item_pow::split(t_frac &_num) const {
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
	t_frac lv, rv;
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
	t_frac lv, rv;
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
	t_frac n =
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
	t_frac lv;
	h_item l2;
	if (lb) l2 = lb->split(lv);
	if (l2 && rn) {
		return (own.gener((lv * lb->num) / rn->num) + l2->mul(lb->num / rn->num))->mul(num);
	}
	if (rn) {
		return lhs->mul(num / rn->num);
	}

	t_frac n =
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
	if (rn && (rn->num == 0)) {
		return own.gener(1);
	}
	if (ln && (ln->num == 1)) {
		return own.gener(1);
	}
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

//...