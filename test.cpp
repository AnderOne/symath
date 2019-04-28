#include "func_tree.hpp"

#include <iostream>
#include <fstream>

int main() {

	t_func_tree TREE;
	std::string str; std::getline(std::cin, str);
	TREE.create(str);
	std::cout << (std::string) TREE << std::endl;
	TREE.derive('x');
	std::cout << (std::string) TREE << std::endl;
	TREE.reduce();
	std::cout << (std::string) TREE << std::endl;

	return 0;
}
