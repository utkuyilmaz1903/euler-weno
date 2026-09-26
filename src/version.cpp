#include "euler/version.hpp"

namespace euler {

std::string_view version() noexcept {
    return EULER_VERSION_STRING; // injected by src/CMakeLists.txt
}

} // namespace euler
