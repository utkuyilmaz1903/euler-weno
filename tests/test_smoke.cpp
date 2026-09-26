#include "euler/version.hpp"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("version string is not empty", "[smoke]") {
    REQUIRE_FALSE(euler::version().empty());
}
