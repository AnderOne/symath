#ifndef __INCLUDE_FUNC_TREE
#define __INCLUDE_FUNC_TREE

#include <type_traits>
#include <ostream>
#include <istream>
#include <string>
#include <memory>

struct t_func_tree {

	inline const double &operator[] (char var) const { return DATA[(var - 'a')]; }
	inline double &operator[] (char var) { return DATA[(var - 'a')]; }

	inline operator std::string() const { return root->str(); }
	inline double get() const { return root->get(); }

	t_func_tree &operator=(const t_func_tree &rhs);

	t_func_tree(const t_func_tree &src);
	t_func_tree(std::string str);
	t_func_tree();

	bool create(std::string str);
	void derive(char var);
	void reduce();

protected:

	struct h_item;
	struct t_item;

	//Фабричные методы для создания разных типов узлов:
	h_item gener(std::string str, const h_item &lhs, const h_item &rhs) const;
	h_item gener(std::string str, const h_item &rhs) const;
	h_item gener(std::string str) const;
	h_item gener(double val) const;

	struct h_item: public std::shared_ptr<t_item> {

		template <typename ... T>
		h_item(T ... args): std::shared_ptr<t_item> (args ...) {}

		inline h_item operator^(const h_item &rhs) const {
			return get()->own.gener("^", *this, rhs);
		}
		inline h_item operator/(const h_item &rhs) const {
			return get()->own.gener("/", *this, rhs);
		}
		inline h_item operator*(const h_item &rhs) const {
			return get()->own.gener("*", *this, rhs);
		}
		inline h_item operator+(const h_item &rhs) const {
			return get()->own.gener("+", *this, rhs);
		}
		inline h_item operator-(const h_item &rhs) const {
			return get()->own.gener("-", *this, rhs);
		}
		inline h_item operator-() const {
			return get()->own.gener("-", *this);
		}
	};

	//Математические функции:
	inline h_item sqrt(const h_item &arg) const {
		return gener("sqrt", arg);
	}
	inline h_item log(const h_item &arg) const {
		return gener("log", arg);
	}
	inline h_item exp(const h_item &arg) const {
		return gener("exp", arg);
	}
	inline h_item cos(const h_item &arg) const {
		return gener("cos", arg);
		
	}
	inline h_item sin(const h_item &arg) const {
		return gener("sin", arg);
	}

	struct t_item {
		explicit t_item(const t_func_tree &_own): own(_own) {}
		virtual std::string str() const = 0;	//Переводит выражение в строку;
		virtual h_item dif(char) const = 0;	//Возвращает производную;
		virtual h_item red() const = 0;	//Сокращает выражение;
		virtual double get() const = 0;	//Вычисляет значение;
		virtual h_item cpy(const t_func_tree &) const = 0;
		virtual h_item cpy() const = 0;
		virtual ~t_item() {}
		const t_func_tree &own;
	};

	struct t_item_const:
	public t_item {
		explicit t_item_const(const t_func_tree &_own, double _val):
		         t_item(_own), val(_val) {}
		std::string str() const override;
		h_item dif(char) const override;
		h_item red() const override;
		double get() const override;
		h_item cpy(const t_func_tree &) const override;
		h_item cpy() const override;
		const double val;
	};

	struct t_item_index:
	public t_item {
		explicit t_item_index(const t_func_tree &_own, char _ind):
		         t_item(_own), ind(_ind) {}
		std::string str() const override;
		h_item dif(char) const override;
		h_item red() const override;
		double get() const override;
		h_item cpy(const t_func_tree &) const override;
		h_item cpy() const override;
		const char ind;
	};

	template <int N>
	struct t_item_arg: public t_item {
		template<typename ... T,
		         typename =
		         typename std::enable_if<sizeof ... (T) == N>
		         ::type>
		explicit t_item_arg(const t_func_tree &_own,
		                    const T & ... _arg):
		         t_item(_own), arg{_arg ...}
		         {}
		explicit t_item_arg(const t_func_tree &_own):
		         t_item(_own)
		         {}
		h_item arg[N];
	};

	//Для каждого вида операции свой тип узла:
	#define __DECL_ITEM_TYPE(func, dim) \
	struct t_item_##func: public t_item_arg<dim> {\
		template <typename ... T>\
		explicit t_item_##func(const T &... _arg):\
		         t_item_arg<dim>(_arg ...) {}\
		std::string str() const override;\
		h_item dif(char) const override;\
		h_item red() const override;\
		double get() const override;\
		h_item cpy(const t_func_tree &)\
		const override;\
		h_item cpy() const override;\
	};

	__DECL_ITEM_TYPE(sqrt, 1)
	__DECL_ITEM_TYPE(log, 1)
	__DECL_ITEM_TYPE(exp, 1)
	__DECL_ITEM_TYPE(cos, 1)
	__DECL_ITEM_TYPE(sin, 1)
	__DECL_ITEM_TYPE(neg, 1)
	__DECL_ITEM_TYPE(sub, 2)
	__DECL_ITEM_TYPE(add, 2)
	__DECL_ITEM_TYPE(mul, 2)
	__DECL_ITEM_TYPE(div, 2)
	__DECL_ITEM_TYPE(pow, 2)

	#undef __DECL_ITEM_TYPE

private:
	double DATA[26];
	h_item root;
};

inline std::ostream &operator<<(std::ostream &out, const t_func_tree &tree) {
	return out << (std::string) tree;
}

inline std::istream &operator>>(std::istream &inp, t_func_tree &tree) {
	std::string str;
	std::getline(inp, str);
	tree.create(str);
	return inp;
}

#endif //__INCLUDE_FUNC_TREE
