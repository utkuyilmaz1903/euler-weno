#pragma once

#include <string_view>

namespace euler {

// Project version as set in the top-level CMakeLists.txt, e.g. "0.1.0".
[[nodiscard]] std::string_view version() noexcept;

} // namespace euler
