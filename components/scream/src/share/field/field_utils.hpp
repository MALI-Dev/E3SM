#ifndef SCREAM_FIELD_UTILS_HPP
#define SCREAM_FIELD_UTILS_HPP

#include "share/field/field_utils_impl.hpp"

namespace scream {

// Check that two fields store the same entries.
// NOTE: if the field is padded, padding entries are NOT checked.
inline bool views_are_equal(const Field& f1, const Field& f2) {
  const auto& dt = f1.get_header().get_identifier().data_type();
  if (dt=="int") {
    return impl::views_are_equal<const int>(f1,f2);
  } else if (dt=="double") {
    return impl::views_are_equal<const double>(f1,f2);
  } else if (dt=="float") {
    return impl::views_are_equal<const float>(f1,f2);
  } else {
    EKAT_ERROR_MSG ("Error! Unsupported field data type '" + dt + "'.\n");
  }
}

template<typename Engine, typename PDF>
void randomize (const Field& f, Engine& engine, PDF&& pdf)
{
  EKAT_REQUIRE_MSG(f.is_allocated(),
      "Error! Cannot randomize the values of a field not yet allocated.\n");

  // Deduce scalar type from pdf
  using ST = decltype(pdf(engine));

  // Check compatibility between PDF and field data type
  const auto& dt = f.get_header().get_identifier().data_type();
  EKAT_REQUIRE_MSG (
      (std::is_same<ST,int>::value && dt=="int") ||
      (std::is_same<ST,float>::value && dt=="float") ||
      (std::is_same<ST,double>::value && dt=="double"),
      "Error! Field data type incompatible with input PDF.\n");

  impl::randomize<ST>(f,engine,pdf);
}

template<typename ST>
ST frobenius_norm(const Field& f)
{
  // Check compatibility between ST and field data type
  const auto& dt = f.get_header().get_identifier().data_type();
  EKAT_REQUIRE_MSG (dt=="float" || dt=="double",
      "Error! Frobenius norm only allowed for floating-point field value types.\n");

  EKAT_REQUIRE_MSG (
      (std::is_same<ST,float>::value && dt=="float") ||
      (std::is_same<ST,double>::value && dt=="double"),
      "Error! Field data type incompatible with template argument.\n");

  return impl::frobenius_norm<ST>(f);
}

template<typename ST>
ST field_sum(const Field& f)
{
  // Check compatibility between ST and field data type
  const auto& dt = f.get_header().get_identifier().data_type();

  EKAT_REQUIRE_MSG (
      (std::is_same<ST,int>::value && dt=="int") ||
      (std::is_same<ST,float>::value && dt=="float") ||
      (std::is_same<ST,double>::value && dt=="double"),
      "Error! Field data type incompatible with template argument.\n");

  return impl::field_sum<ST>(f);
}

template<typename ST>
ST field_max(const Field& f)
{
  // Check compatibility between ST and field data type
  const auto& dt = f.get_header().get_identifier().data_type();

  EKAT_REQUIRE_MSG (
      (std::is_same<ST,int>::value && dt=="int") ||
      (std::is_same<ST,float>::value && dt=="float") ||
      (std::is_same<ST,double>::value && dt=="double"),
      "Error! Field data type incompatible with template argument.\n");

  return impl::field_max<ST>(f);
}

template<typename ST>
ST field_min(const Field& f)
{
  // Check compatibility between ST and field data type
  const auto& dt = f.get_header().get_identifier().data_type();

  EKAT_REQUIRE_MSG (
      (std::is_same<ST,int>::value && dt=="int") ||
      (std::is_same<ST,float>::value && dt=="float") ||
      (std::is_same<ST,double>::value && dt=="double"),
      "Error! Field data type incompatible with template argument.\n");

  return impl::field_min<ST>(f);
}

} // namespace scream

#endif // SCREAM_FIELD_UTILS_HPP
