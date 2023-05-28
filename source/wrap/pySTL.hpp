#ifndef _PY_STL_HPP_
#define _PY_STL_HPP_

#include <boost/python.hpp>

namespace python = boost::python;

template<class T>
std::vector<T> to_vector(const python::object& iterable) {
  return std::vector<T>(
    python::stl_input_iterator<T>(iterable), python::stl_input_iterator<T>());
};

template <class T>
python::list to_list(const std::vector<T>& vector) {
  python::list list;
  for (const T& value : vector) {
    list.append(value);
  };
  return list;
};

template <class Map>
python::dict to_dict(const Map& map) {
  python::dict dictionary;
  for (const auto& [key, value] : map) {
    dictionary[key] = value;
  };
  return dictionary;
};

#endif // !_PY_STL_HPP_