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

#ifndef __KERNEL_H__
#define __KERNEL_H__

/* CPU Kernel Interface */

#include "util_types.h"

CCL_NAMESPACE_BEGIN

#define KERNEL_NAME_JOIN(x, y, z) x ## _ ## y ## _ ## z
#define KERNEL_NAME_EVAL(arch, name)  KERNEL_NAME_JOIN(kernel, arch, name)
#define KERNEL_FUNCTION_FULL_NAME(name) KERNEL_NAME_EVAL(KERNEL_ARCH, name)

struct KernelGlobals;

KernelGlobals *kernel_globals_create();
void kernel_globals_free(KernelGlobals *kg);

void *kernel_osl_memory(KernelGlobals *kg);
bool kernel_osl_use(KernelGlobals *kg);

void kernel_const_copy(KernelGlobals *kg, const char *name, void *host, size_t size);
void kernel_tex_copy(KernelGlobals *kg,
                     const char *name,
                     device_ptr mem,
                     size_t width,
                     size_t height,
                     size_t depth,
                     InterpolationType interpolation=INTERPOLATION_LINEAR,
                     ExtensionType extension = EXTENSION_REPEAT);

#define KERNEL_ARCH cpu
#include "kernels/cpu/kernel_cpu.h"

#ifdef WITH_CYCLES_OPTIMIZED_KERNEL_SSE2
#  define KERNEL_ARCH cpu_sse2
#  include "kernels/cpu/kernel_cpu.h"
#endif  /* WITH_CYCLES_OPTIMIZED_KERNEL_SSE2 */

#ifdef WITH_CYCLES_OPTIMIZED_KERNEL_SSE3
#  define KERNEL_ARCH cpu_sse3
#  include "kernels/cpu/kernel_cpu.h"
#endif  /* WITH_CYCLES_OPTIMIZED_KERNEL_SSE2 */

#ifdef WITH_CYCLES_OPTIMIZED_KERNEL_SSE41
#  define KERNEL_ARCH cpu_sse41
#  include "kernels/cpu/kernel_cpu.h"
#endif  /* WITH_CYCLES_OPTIMIZED_KERNEL_SSE41 */

#ifdef WITH_CYCLES_OPTIMIZED_KERNEL_AVX
#  define KERNEL_ARCH cpu_avx
#  include "kernels/cpu/kernel_cpu.h"
#endif  /* WITH_CYCLES_OPTIMIZED_KERNEL_AVX */

#ifdef WITH_CYCLES_OPTIMIZED_KERNEL_AVX2
#  define KERNEL_ARCH cpu_avx2
#  include "kernels/cpu/kernel_cpu.h"
#endif  /* WITH_CYCLES_OPTIMIZED_KERNEL_AVX2 */

CCL_NAMESPACE_END

#endif /* __KERNEL_H__ */

