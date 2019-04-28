#include "func_tree.hpp"

#include <iostream>
#include <fstream>

int main() {

	t_func_tree TREE;
	std::string str; std::getline(std::cin, str);
	TREE.create(str);
	std::cout << TREE << std::endl;
	TREE.derive('x');
	std::cout << TREE << std::endl;
	TREE.reduce();
	std::cout << TREE << std::endl;

	return 0;
}
