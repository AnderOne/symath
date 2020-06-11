#define BOOST_TEST_MODULE test_tree

#include <boost/test/unit_test.hpp>
#include <symath/tree.hpp>
#include <exception>
#include <string>
#include <vector>

BOOST_AUTO_TEST_SUITE(suite_of_tree_tests)

BOOST_AUTO_TEST_CASE(test_parser_error) {

	const std::vector<std::pair<std::string, std::string>> TEST = {
		{"((123 - 4)", "There is extra opening bracket!"},
		{"(123 - 4))", "There is extra closing bracket!"},
		{"abc", "Incorrect name of the operand!"},
		{"2 x", "Incorrect location of operand!"},
		{"x 2", "Incorrect location of operand!"},
		{"2 (", "Incorrect location of bracket!"},
		{"x (", "Incorrect location of bracket!"},
		{"+", "Incorrect location of operation!"},
		{"@", "There is incorrect token in the string!"}
	};

	for (int i = 0; i < TEST.size(); ++ i) {
		symath::t_tree TREE;
		auto check = [&TEST, i](const std::logic_error &err) {
			return (err.what() == TEST[i].second);
		};
		BOOST_REQUIRE_EXCEPTION(
			TREE.create(TEST[i].first),
			std::logic_error,
			check
		);
	}

}

//...

BOOST_AUTO_TEST_SUITE_END()
