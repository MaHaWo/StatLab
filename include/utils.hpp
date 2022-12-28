#ifndef STATLAB_UTILS_HPP
#define STATLAB_UTILS_HPP
namespace StatLab {
namespace Utils {
/**
 * @brief Indicate how NaN or Inf shall be treated:
 *  - propagate: default, this propagates nan or inf to the result
 *  - skip: skip them
 *  - error: throw invalid_agument exception when nan or inf is encountered
 */
enum class NaNPolicy
{
    propagate,
    skip,
    error
};


}
}
#endif