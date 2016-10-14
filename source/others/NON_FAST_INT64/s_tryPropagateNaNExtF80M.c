
/*============================================================================

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions, and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions, and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

 3. Neither the name of the University nor the names of its contributors may
    be used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS", AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, ARE
DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=============================================================================*/

#include <stdbool.h>

#include "internals.h"
#include "specialize.h"

bool
 softfloat_tryPropagateNaNExtF80M(
    /** @bug use extFloat80_t */
    const struct extFloat80M *aSPtr,
    /** @bug use extFloat80_t */
    const struct extFloat80M *bSPtr,
    /** @bug use extFloat80_t */
    struct extFloat80M *zSPtr
 )
{
    uint16_t ui64;
    uint64_t ui0;

    ui64 = aSPtr->signExp;
    ui0  = aSPtr->signif;
    if ( isNaNExtF80UI( ui64, ui0 ) ) goto propagateNaN;
    ui64 = bSPtr->signExp;
    ui0  = bSPtr->signif;
    if ( isNaNExtF80UI( ui64, ui0 ) ) goto propagateNaN;
    return false;
 propagateNaN:
    softfloat_propagateNaNExtF80M( aSPtr, bSPtr, zSPtr );
    return true;

}

