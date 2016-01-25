/*
 * Copyright 2011-2013 Blender Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* Optimized CPU kernel entry points. This file is compiled with AVX
 * optimization flags and nearly all functions inlined, while kernel.cpp
 * is compiled without for other CPU's. */
 
/* SSE optimization disabled for now on 32 bit, see bug #36316 */
#if !(defined(__GNUC__) && (defined(i386) || defined(_M_IX86)))
#define __KERNEL_SSE2__
#define __KERNEL_SSE3__
#define __KERNEL_SSSE3__
#define __KERNEL_SSE41__
#define __KERNEL_AVX__
#endif
 
#include "util_optimization.h"
 
#ifdef WITH_CYCLES_OPTIMIZED_KERNEL_AVX
#  include "kernel.h"
#  define KERNEL_ARCH cpu_avx
#  include "kernel_cpu_impl.h"

#endif  /* WITH_CYCLES_OPTIMIZED_KERNEL_AVX */
