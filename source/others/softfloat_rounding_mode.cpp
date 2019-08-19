#include "softfloat/functions.h"

namespace {
static thread_local softfloat_round_mode softfloat_roundingMode = softfloat_round_near_even;
}

softfloat_round_mode
softfloat_get_roundingMode(void)
{
    return softfloat_roundingMode;
}

void
softfloat_set_roundingMode(softfloat_round_mode mode)
{
    softfloat_roundingMode = mode;
}
