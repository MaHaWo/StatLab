#ifndef UTOPIA_MODELS_AMEEMULTI_TEST_UTILS_HH
#define UTOPIA_MODELS_AMEEMULTI_TEST_UTILS_HH

// custom assert macros go here
#include "../include/utils.hh"
#include <cmath>
#include <iomanip>

namespace Statistics
{
namespace TestUtils
{
#define ASSERT_EQ(lhs, rhs)                                                  \
    if (!Statistics::Utils::IsEqual()(lhs, rhs))                             \
    {                                                                        \
        std::cerr << std::setprecision(19)                                   \
                  << "Asserted equality wrong at line: " << __LINE__;        \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;               \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl; \
        std::exit(-1);                                                       \
    }

#define ASSERT_EQ_CUSTOM(lhs, rhs, tol)                                            \
    if (!Statistics::Utils::IsEqual()(lhs, rhs, tol))                              \
    {                                                                              \
        std::cerr << std::setprecision(19)                                         \
                  << "Asserted equality wrong at line: " << __LINE__ << std::endl; \
        std::cerr << std::setprecision(19) << " lhs: " << lhs << std::endl;        \
        std::cerr << std::setprecision(19) << " rhs: " << rhs << std::endl;        \
        std::cerr << std::setprecision(19)                                         \
                  << " abs. difference: " << std::abs(lhs - rhs) << std::endl;     \
        std::exit(-1);                                                             \
    }

#define EXPECT_EQ(lhs, rhs)                                                  \
    if (!Statistics::Utils::IsEqual()(lhs, rhs))                             \
    {                                                                        \
        std::cerr << std::setprecision(19)                                   \
                  << "Exepcted equality wrong at line: " << __LINE__;        \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;               \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl; \
    }

#define ASSERT_NEQ(lhs, rhs)                                                 \
    if (Statistics::Utils::IsEqual()(lhs, rhs))                              \
    {                                                                        \
        std::cerr << std::setprecision(19)                                   \
                  << "Asserted inequality wrong at line: " << __LINE__;      \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;               \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl; \
        std::exit(-1);                                                       \
    }

#define EXPECT_NEQ(lhs, rhs)                                                 \
    if (Statistics::Utils::IsEqual()(lhs, rhs))                              \
    {                                                                        \
        std::cerr << std::setprecision(19)                                   \
                  << "Exepcted inequality wrong at line: " << __LINE__;      \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;               \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl; \
    }

#define ASSERT_GREATER(lhs, rhs)                                                  \
    if (!(Statistics::Utils::IsGreater()(lhs, rhs)))                              \
    {                                                                             \
        std::cerr << std::setprecision(19)                                        \
                  << "Asserted relation 'lhs > rhs' wrong at line: " << __LINE__; \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;                    \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl;      \
        std::exit(-1);                                                            \
    }

#define EXPECT_GREATER(lhs, rhs)                                                  \
    if (!(Statistics::Utils::IsGreater()(lhs, rhs)))                              \
    {                                                                             \
        std::cerr << "Exepcted relation 'lhs > rhs' wrong at line: " << __LINE__; \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;                    \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl;      \
    }

#define ASSERT_GEQ(lhs, rhs)                                                       \
    if (!(Statistics::Utils::IsGreater()(lhs, rhs) or                              \
          Statistics::Utils::IsEqual()(lhs, rhs)))                                 \
    {                                                                              \
        std::cerr << "Asserted relation 'lhs <= rhs' wrong at line: " << __LINE__; \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;                     \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl;       \
        std::exit(-1);                                                             \
    }

#define EXPECT_GEQ(lhs, rhs)                                                       \
    if (!(Statistics::Utils::IsGreater()(lhs, rhs) or                              \
          Statistics::Utils::IsEqual()(lhs, rhs)))                                 \
    {                                                                              \
        std::cerr << "Exepcted relation 'lhs <= rhs' wrong at line: " << __LINE__; \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;                     \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl;       \
    }

#define ASSERT_LESS(lhs, rhs)                                                     \
    if (!(Statistics::Utils::IsLess()(lhs, rhs)))                                 \
    {                                                                             \
        std::cerr << std::setprecision(19)                                        \
                  << "Asserted relation 'lhs < rhs' wrong at line: " << __LINE__; \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;                    \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl;      \
        std::exit(-1);                                                            \
    }

#define EXPECT_LESS(lhs, rhs)                                                     \
    if (!(Statistics::Utils::IsLess()(lhs, rhs)))                                 \
    {                                                                             \
        std::cerr << std::setprecision(19)                                        \
                  << "Exepcted relation 'lhs < rhs' wrong at line: " << __LINE__; \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;                    \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl;      \
    }

#define ASSERT_LEQ(lhs, rhs)                                                                \
    if (!(Statistics::Utils::IsLess()(lhs, rhs) or Statistics::Utils::IsEqual()(lhs, rhs))) \
    {                                                                                       \
        std::cerr << std::setprecision(19)                                                  \
                  << "Asserted relation 'lhs <= rhs' wrong at line: " << __LINE__;          \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;                              \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl;                \
        std::exit(-1);                                                                      \
    }

#define EXPECT_LEQ(lhs, rhs)                                                                \
    if (!(Statistics::Utils::IsLess()(lhs, rhs) or Statistics::Utils::IsEqual()(lhs, rhs))) \
    {                                                                                       \
        std::cerr << std::setprecision(19)                                                  \
                  << "Exepcted relation 'lhs <= rhs' wrong at line: " << __LINE__;          \
        std::cerr << std::setprecision(19) << " lhs: " << lhs;                              \
        std::cerr << std::setprecision(19) << ", rhs: " << rhs << std::endl;                \
    }

} // namespace TestUtils
} // namespace Statistics
#endif
