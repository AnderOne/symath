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

#include <symath/long.hpp>

namespace symath {

std::string format(const t_frac &val, t_form frm) {

	std::string str = (val < 0)? ("- "): (""); if (val.isint()) return str + std::string(abs(val));
	t_long a = abs(val).upper();
	t_long b = abs(val).lower(); if (b == 0) return str + "INF";
	if (frm & FRM_DOT) {
		t_long q = b;
		int n, m, c;
		int n5 = 0; while (q % 5 == 0) { ++ n5; q /= 5; }
		int n2 = 0; while (q % 2 == 0) { ++ n2; q /= 2; }
		if (q == 1) {
			if (n5 < n2) {
				n = (m = n2) - n5; c = 5;
			}
			else {
				n = (m = n5) - n2; c = 2;
			}
			if (a > b) {
				str += std::string(t_long(a / b));
				a = (a % b);
			}
			else {
				str += '0';
			}
			std::string tmp =
			std::string(a * pow(t_long(c), n));
			str += '.';
			for (m -= tmp.size(); m; -- m) {
			str += '0';
			}
			str += tmp;
			return str;
		}
	}
	if (frm & FRM_RED) {
		if (a > b) {
			str += std::string(t_long(a / b));
			a = (a % b);
			str += "[" +
			       std::string(a) + "/" +
			       std::string(b) +
			       "]";
			return str;
		}
	}
	str += std::string(a) + "/" +
	       std::string(b);
	return str;
}

//...

}
