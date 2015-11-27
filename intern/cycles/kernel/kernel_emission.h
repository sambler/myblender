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

CCL_NAMESPACE_BEGIN

/* Direction Emission */
ccl_device_noinline float3 direct_emissive_eval(KernelGlobals *kg,
	LightSample *ls, float3 I, differential3 dI, float t, float time, int bounce, int transparent_bounce
#ifdef __SPLIT_KERNEL__
	,ShaderData *sd_input
#endif
)
{
	/* setup shading at emitter */
#ifdef __SPLIT_KERNEL__
	ShaderData *sd = sd_input;
#else
	ShaderData sd_object;
	ShaderData *sd = &sd_object;
#endif
	float3 eval;

#ifdef __BACKGROUND_MIS__
	if(ls->type == LIGHT_BACKGROUND) {
		Ray ray;
		ray.D = ls->D;
		ray.P = ls->P;
		ray.t = 1.0f;
#ifdef __OBJECT_MOTION__
		ray.time = time;
#endif
		ray.dP = differential3_zero();
		ray.dD = dI;

		shader_setup_from_background(kg, sd, &ray, bounce+1, transparent_bounce);
		eval = shader_eval_background(kg, sd, 0, SHADER_CONTEXT_EMISSION);
	}
	else
#endif
	{
		shader_setup_from_sample(kg, sd, ls->P, ls->Ng, I, ls->shader, ls->object, ls->prim, ls->u, ls->v, t, time, bounce+1, transparent_bounce);

		ls->Ng = ccl_fetch(sd, Ng);

		/* no path flag, we're evaluating this for all closures. that's weak but
		 * we'd have to do multiple evaluations otherwise */
		shader_eval_surface(kg, sd, 0.0f, 0, SHADER_CONTEXT_EMISSION);

		/* evaluate emissive closure */
		if(ccl_fetch(sd, flag) & SD_EMISSION)
			eval = shader_emissive_eval(kg, sd);
		else
			eval = make_float3(0.0f, 0.0f, 0.0f);
	}
	
	eval *= ls->eval_fac;

	return eval;
}

ccl_device_noinline bool direct_emission(KernelGlobals *kg, ShaderData *sd,
	LightSample *ls, Ray *ray, BsdfEval *eval, bool *is_lamp,
	int bounce, int transparent_bounce
#ifdef __SPLIT_KERNEL__
	, ShaderData *sd_DL
#endif
	)
{
	if(ls->pdf == 0.0f)
		return false;

	/* todo: implement */
	differential3 dD = differential3_zero();

	/* evaluate closure */

	float3 light_eval = direct_emissive_eval(kg, ls, -ls->D, dD, ls->t, ccl_fetch(sd, time),
	                                         bounce,
	                                         transparent_bounce
#ifdef __SPLIT_KERNEL__
	                                         ,sd_DL
#endif
	                                         );

	if(is_zero(light_eval))
		return false;

	/* evaluate BSDF at shading point */

#ifdef __VOLUME__
	if(ccl_fetch(sd, prim) != PRIM_NONE)
		shader_bsdf_eval(kg, sd, ls->D, eval, ls->pdf, ls->shader & SHADER_USE_MIS);
	else {
		float bsdf_pdf;
		shader_volume_phase_eval(kg, sd, ls->D, eval, &bsdf_pdf);
		if(ls->shader & SHADER_USE_MIS) {
			/* Multiple importance sampling. */
			float mis_weight = power_heuristic(ls->pdf, bsdf_pdf);
			light_eval *= mis_weight;
		}
	}
#else
	shader_bsdf_eval(kg, sd, ls->D, eval, ls->pdf, ls->shader & SHADER_USE_MIS);
#endif

	bsdf_eval_mul(eval, light_eval/ls->pdf);

#ifdef __PASSES__
	/* use visibility flag to skip lights */
	if(ls->shader & SHADER_EXCLUDE_ANY) {
		if(ls->shader & SHADER_EXCLUDE_DIFFUSE)
			eval->diffuse = make_float3(0.0f, 0.0f, 0.0f);
		if(ls->shader & SHADER_EXCLUDE_GLOSSY)
			eval->glossy = make_float3(0.0f, 0.0f, 0.0f);
		if(ls->shader & SHADER_EXCLUDE_TRANSMIT)
			eval->transmission = make_float3(0.0f, 0.0f, 0.0f);
		if(ls->shader & SHADER_EXCLUDE_SCATTER)
			eval->scatter = make_float3(0.0f, 0.0f, 0.0f);
	}
#endif

	if(bsdf_eval_is_zero(eval))
		return false;

	if(ls->shader & SHADER_CAST_SHADOW) {
		/* setup ray */
		bool transmit = (dot(ccl_fetch(sd, Ng), ls->D) < 0.0f);
		ray->P = ray_offset(ccl_fetch(sd, P), (transmit)? -ccl_fetch(sd, Ng): ccl_fetch(sd, Ng));

		if(ls->t == FLT_MAX) {
			/* distant light */
			ray->D = ls->D;
			ray->t = ls->t;
		}
		else {
			/* other lights, avoid self-intersection */
			ray->D = ray_offset(ls->P, ls->Ng) - ray->P;
			ray->D = normalize_len(ray->D, &ray->t);
		}

		ray->dP = ccl_fetch(sd, dP);
		ray->dD = differential3_zero();
	}
	else {
		/* signal to not cast shadow ray */
		ray->t = 0.0f;
	}

	/* return if it's a lamp for shadow pass */
	*is_lamp = (ls->prim == PRIM_NONE && ls->type != LIGHT_BACKGROUND);

	return true;
}

/* Indirect Primitive Emission */

ccl_device_noinline float3 indirect_primitive_emission(KernelGlobals *kg, ShaderData *sd, float t, int path_flag, float bsdf_pdf)
{
	/* evaluate emissive closure */
	float3 L = shader_emissive_eval(kg, sd);

#ifdef __HAIR__
	if(!(path_flag & PATH_RAY_MIS_SKIP) && (ccl_fetch(sd, flag) & SD_USE_MIS) && (ccl_fetch(sd, type) & PRIMITIVE_ALL_TRIANGLE))
#else
	if(!(path_flag & PATH_RAY_MIS_SKIP) && (ccl_fetch(sd, flag) & SD_USE_MIS))
#endif
	{
		/* multiple importance sampling, get triangle light pdf,
		 * and compute weight with respect to BSDF pdf */
		float pdf = triangle_light_pdf(kg, ccl_fetch(sd, Ng), ccl_fetch(sd, I), t);
		float mis_weight = power_heuristic(bsdf_pdf, pdf);

		return L*mis_weight;
	}

	return L;
}

/* Indirect Lamp Emission */

ccl_device_noinline bool indirect_lamp_emission(KernelGlobals *kg, PathState *state, Ray *ray, float3 *emission
#ifdef __SPLIT_KERNEL__
                                                ,ShaderData *sd
#endif
                                                )
{
	bool hit_lamp = false;

	*emission = make_float3(0.0f, 0.0f, 0.0f);

	for(int lamp = 0; lamp < kernel_data.integrator.num_all_lights; lamp++) {
		LightSample ls;

		if(!lamp_light_eval(kg, lamp, ray->P, ray->D, ray->t, &ls))
			continue;

#ifdef __PASSES__
		/* use visibility flag to skip lights */
		if(ls.shader & SHADER_EXCLUDE_ANY) {
			if(((ls.shader & SHADER_EXCLUDE_DIFFUSE) && (state->flag & PATH_RAY_DIFFUSE)) ||
			   ((ls.shader & SHADER_EXCLUDE_GLOSSY) &&
			    ((state->flag & (PATH_RAY_GLOSSY|PATH_RAY_REFLECT)) == (PATH_RAY_GLOSSY|PATH_RAY_REFLECT))) ||
			   ((ls.shader & SHADER_EXCLUDE_TRANSMIT) && (state->flag & PATH_RAY_TRANSMIT)) ||
			   ((ls.shader & SHADER_EXCLUDE_SCATTER) && (state->flag & PATH_RAY_VOLUME_SCATTER)))
				continue;
		}
#endif

		float3 L = direct_emissive_eval(kg, &ls, -ray->D, ray->dD, ls.t, ray->time,
		                                state->bounce,
		                                state->transparent_bounce
#ifdef __SPLIT_KERNEL__
		                                ,sd
#endif
		                                );

#ifdef __VOLUME__
		if(state->volume_stack[0].shader != SHADER_NONE) {
			/* shadow attenuation */
			Ray volume_ray = *ray;
			volume_ray.t = ls.t;
			float3 volume_tp = make_float3(1.0f, 1.0f, 1.0f);
			kernel_volume_shadow(kg, state, &volume_ray, &volume_tp);
			L *= volume_tp;
		}
#endif

		if(!(state->flag & PATH_RAY_MIS_SKIP)) {
			/* multiple importance sampling, get regular light pdf,
			 * and compute weight with respect to BSDF pdf */
			float mis_weight = power_heuristic(state->ray_pdf, ls.pdf);
			L *= mis_weight;
		}

		*emission += L;
		hit_lamp = true;
	}

	return hit_lamp;
}

/* Indirect Background */

ccl_device_noinline float3 indirect_background(KernelGlobals *kg, ccl_addr_space PathState *state, ccl_addr_space Ray *ray
#ifdef __SPLIT_KERNEL__
                                               ,ShaderData *sd_global
#endif
                                               )
{
#ifdef __BACKGROUND__
	int shader = kernel_data.background.surface_shader;

	/* use visibility flag to skip lights */
	if(shader & SHADER_EXCLUDE_ANY) {
		if(((shader & SHADER_EXCLUDE_DIFFUSE) && (state->flag & PATH_RAY_DIFFUSE)) ||
		   ((shader & SHADER_EXCLUDE_GLOSSY) &&
		    ((state->flag & (PATH_RAY_GLOSSY|PATH_RAY_REFLECT)) == (PATH_RAY_GLOSSY|PATH_RAY_REFLECT))) ||
		   ((shader & SHADER_EXCLUDE_TRANSMIT) && (state->flag & PATH_RAY_TRANSMIT)) ||
		   ((shader & SHADER_EXCLUDE_CAMERA) && (state->flag & PATH_RAY_CAMERA)) ||
		   ((shader & SHADER_EXCLUDE_SCATTER) && (state->flag & PATH_RAY_VOLUME_SCATTER)))
			return make_float3(0.0f, 0.0f, 0.0f);
	}

#ifdef __SPLIT_KERNEL__
	/* evaluate background closure */
	Ray priv_ray = *ray;
	shader_setup_from_background(kg, sd_global, &priv_ray, state->bounce+1, state->transparent_bounce);
	float3 L = shader_eval_background(kg, sd_global, state->flag, SHADER_CONTEXT_EMISSION);
#else
	ShaderData sd;
	shader_setup_from_background(kg, &sd, ray, state->bounce+1, state->transparent_bounce);

	float3 L = shader_eval_background(kg, &sd, state->flag, SHADER_CONTEXT_EMISSION);
#endif

#ifdef __BACKGROUND_MIS__
	/* check if background light exists or if we should skip pdf */
	int res = kernel_data.integrator.pdf_background_res;

	if(!(state->flag & PATH_RAY_MIS_SKIP) && res) {
		/* multiple importance sampling, get background light pdf for ray
		 * direction, and compute weight with respect to BSDF pdf */
		float pdf = background_light_pdf(kg, ray->P, ray->D);
		float mis_weight = power_heuristic(state->ray_pdf, pdf);

		return L*mis_weight;
	}
#endif

	return L;
#else
	return make_float3(0.8f, 0.8f, 0.8f);
#endif
}

CCL_NAMESPACE_END

