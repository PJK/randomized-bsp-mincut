#ifndef PARALLEL_MINIMUM_CUT_INPUTITERATOR_HPP
#define PARALLEL_MINIMUM_CUT_INPUTITERATOR_HPP

#include <string>
#include <iostream>
#include <fstream>
#include "AdjacencyListGraph.hpp"

class InputIterator {
	std::ifstream file_;
	unsigned lines_, read_;
	unsigned vertices_;
	std::string name_;

	void open();
public:
	InputIterator(std::string name) : name_(name)
	{
		file_.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		open();
	}

	~InputIterator()
	{
		file_.close();
	}

	// Model of http://en.cppreference.com/w/cpp/concept/InputIterator concept
	class Iterator {
		unsigned position_ = 0;
	public:
		bool end_;
		InputIterator & parent_;
		AdjacencyListGraph::Edge edge_;

		Iterator(bool end, InputIterator & parent, AdjacencyListGraph::Edge edge) : end_(end), parent_(parent), edge_(edge) {}

		AdjacencyListGraph::Edge operator*()
		{
			return edge_;
		}

		AdjacencyListGraph::Edge * operator->()
		{
			return &edge_;
		}

		unsigned position() const {
			return position_;
		}

		void operator++() {
			if (parent_.read_ < parent_.lines_)
				edge_ = parent_.read();
			else
				end_ = true;
			position_++;
		}

		bool operator!=(Iterator & other) {
			return other.end_ != end_;
		}
	};

	Iterator begin();

	Iterator end();

	/**
	 * Reset the file stream in order to read the input again
	 */
	void reopen();

	unsigned vertexCount() { return vertices_; }
	unsigned edgeCount() { return lines_; }

protected:
	AdjacencyListGraph::Edge read();
};

namespace std {
	template<>
	struct iterator_traits<InputIterator::Iterator> {
		typedef std::ptrdiff_t difference_type; //almost always ptrdif_t
		typedef AdjacencyListGraph::Edge value_type; //almost always T
		typedef AdjacencyListGraph::Edge &reference; //almost always T& or const T&
		typedef AdjacencyListGraph::Edge *pointer; //almost always T* or const T*
		typedef std::input_iterator_tag iterator_category;  //usually std::forward_iterator_tag or similar
	};
}


#endif //PARALLEL_MINIMUM_CUT_INPUTITERATOR_HPP
