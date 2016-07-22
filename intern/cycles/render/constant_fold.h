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

#ifndef __CONSTANT_FOLD_H__
#define __CONSTANT_FOLD_H__

#include "util_types.h"

CCL_NAMESPACE_BEGIN

class ShaderGraph;
class ShaderInput;
class ShaderNode;
class ShaderOutput;

class ConstantFolder {
public:
	ShaderGraph *const graph;
	ShaderNode *const node;
	ShaderOutput *const output;

	ConstantFolder(ShaderGraph *graph, ShaderNode *node, ShaderOutput *output);

	bool all_inputs_constant() const;

	/* Constant folding helpers, always return true for convenience. */
	void make_constant(float value) const;
	void make_constant(float3 value) const;
	void make_constant_clamp(float value, bool clamp) const;
	void make_constant_clamp(float3 value, bool clamp) const;

	/* Bypass node, relinking to another output socket. */
	void bypass(ShaderOutput *output) const;

	/* For closure nodes, discard node entirely or bypass to one of its inputs. */
	void discard() const;
	void bypass_or_discard(ShaderInput *input) const;

	/* Bypass or make constant, unless we can't due to clamp being true. */
	bool try_bypass_or_make_constant(ShaderInput *input, float3 input_value, bool clamp) const;
};

CCL_NAMESPACE_END

#endif /* __CONSTANT_FOLD_H__ */

