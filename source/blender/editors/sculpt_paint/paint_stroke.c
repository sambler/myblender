/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2009 by Nicholas Bishop
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): Jason Wilkins, Tom Musgrove.
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/editors/sculpt_paint/paint_stroke.c
 *  \ingroup edsculpt
 */


#include "MEM_guardedalloc.h"

#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BLI_rand.h"
#include "BLI_listbase.h"

#include "DNA_object_types.h"
#include "DNA_scene_types.h"
#include "DNA_brush_types.h"
#include "DNA_curve_types.h"

#include "RNA_access.h"

#include "BKE_context.h"
#include "BKE_paint.h"
#include "BKE_brush.h"
#include "BKE_curve.h"
#include "BKE_colortools.h"
#include "BKE_image.h"

#include "WM_api.h"
#include "WM_types.h"

#include "BIF_gl.h"
#include "BIF_glutil.h"

#include "ED_screen.h"
#include "ED_view3d.h"

#include "IMB_imbuf_types.h"

#include "paint_intern.h"

#include <float.h>
#include <math.h>

typedef struct PaintSample {
	float mouse[2];
	float pressure;
} PaintSample;

typedef struct PaintStroke {
	void *mode_data;
	void *stroke_cursor;
	wmTimer *timer;

	/* Cached values */
	ViewContext vc;
	bglMats mats;
	Brush *brush;
	UnifiedPaintSettings *ups;

	/* used for lines and curves */
	ListBase line;

	/* Paint stroke can use up to PAINT_MAX_INPUT_SAMPLES prior inputs
	 * to smooth the stroke */
	PaintSample samples[PAINT_MAX_INPUT_SAMPLES];
	int num_samples;
	int cur_sample;

	float last_mouse_position[2];
	/* space distance covered so far */
	float stroke_distance;

	/* Set whether any stroke step has yet occurred
	 * e.g. in sculpt mode, stroke doesn't start until cursor
	 * passes over the mesh */
	bool stroke_started;
	/* event that started stroke, for modal() return */
	int event_type;
	/* check if stroke variables have been initialized */
	bool stroke_init;
	/* check if various brush mapping variables have been initialized */
	bool brush_init;
	float initial_mouse[2];
	/* cached_pressure stores initial pressure for size pressure influence mainly */
	float cached_size_pressure;
	/* last pressure will store last pressure value for use in interpolation for space strokes */
	float last_pressure;
	int stroke_mode;

	float zoom_2d;
	int pen_flip;

	StrokeGetLocation get_location;
	StrokeTestStart test_start;
	StrokeUpdateStep update_step;
	StrokeRedraw redraw;
	StrokeDone done;
} PaintStroke;

/*** Cursors ***/
static void paint_draw_smooth_cursor(bContext *C, int x, int y, void *customdata)
{
	Paint *paint = BKE_paint_get_active_from_context(C);
	Brush *brush = BKE_paint_brush(paint);
	PaintStroke *stroke = customdata;

	if (stroke && brush) {
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glColor4ubv(paint->paint_cursor_col);
		sdrawline(x, y, (int)stroke->last_mouse_position[0],
		          (int)stroke->last_mouse_position[1]);
		glDisable(GL_BLEND);
		glDisable(GL_LINE_SMOOTH);
	}
}

static void paint_draw_line_cursor(bContext *C, int x, int y, void *customdata)
{
	Paint *paint = BKE_paint_get_active_from_context(C);
	PaintStroke *stroke = customdata;

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);

	glEnable(GL_LINE_STIPPLE);
	glLineStipple(3, 0xAAAA);

	glColor4ub(0, 0, 0, paint->paint_cursor_col[3]);
	glLineWidth(3.0);
	sdrawline((int)stroke->last_mouse_position[0], (int)stroke->last_mouse_position[1],
	        x, y);

	glColor4ub(255, 255, 255, paint->paint_cursor_col[3]);
	glLineWidth(1.0);
	sdrawline((int)stroke->last_mouse_position[0], (int)stroke->last_mouse_position[1],
	        x, y);

	glDisable(GL_LINE_STIPPLE);

	glDisable(GL_BLEND);
	glDisable(GL_LINE_SMOOTH);
}

static bool paint_tool_require_location(Brush *brush, PaintMode mode)
{
	switch (mode) {
		case PAINT_SCULPT:
			if (ELEM(brush->sculpt_tool, SCULPT_TOOL_GRAB, SCULPT_TOOL_ROTATE,
			                             SCULPT_TOOL_SNAKE_HOOK, SCULPT_TOOL_THUMB))
			{
				return false;
			}
			else {
				return true;
			}
		default:
			break;
	}

	return true;
}

/* Initialize the stroke cache variants from operator properties */
static bool paint_brush_update(bContext *C,
                               Brush *brush,
                               PaintMode mode,
                               struct PaintStroke *stroke,
                               const float mouse_init[2],
                               float mouse[2], float pressure,
                               float location[3])
{
	Scene *scene = CTX_data_scene(C);
	UnifiedPaintSettings *ups = stroke->ups;
	bool location_sampled = false;
	bool location_success = false;
	/* XXX: Use pressure value from first brush step for brushes which don't
	 *      support strokes (grab, thumb). They depends on initial state and
	 *      brush coord/pressure/etc.
	 *      It's more an events design issue, which doesn't split coordinate/pressure/angle
	 *      changing events. We should avoid this after events system re-design */
	if (!stroke->brush_init) {
		copy_v2_v2(stroke->initial_mouse, mouse);
		copy_v2_v2(ups->last_rake, mouse);
		copy_v2_v2(ups->tex_mouse, mouse);
		copy_v2_v2(ups->mask_tex_mouse, mouse);
		stroke->cached_size_pressure = pressure;

		ups->do_linear_conversion = false;
		ups->colorspace = NULL;

		/* check here if color sampling the main brush should do color conversion. This is done here
		 * to avoid locking up to get the image buffer during sampling */
		if (brush->mtex.tex && brush->mtex.tex->type == TEX_IMAGE && brush->mtex.tex->ima) {
			ImBuf *tex_ibuf = BKE_image_pool_acquire_ibuf(brush->mtex.tex->ima, &brush->mtex.tex->iuser, NULL);
			if (tex_ibuf && tex_ibuf->rect_float == NULL) {
				ups->do_linear_conversion = true;
				ups->colorspace = tex_ibuf->rect_colorspace;
			}
			BKE_image_pool_release_ibuf(brush->mtex.tex->ima, tex_ibuf, NULL);
		}

		stroke->brush_init = true;
	}

	if (paint_supports_dynamic_size(brush, mode)) {
		copy_v2_v2(ups->tex_mouse, mouse);
		copy_v2_v2(ups->mask_tex_mouse, mouse);
		stroke->cached_size_pressure = pressure;
	}

	/* Truly temporary data that isn't stored in properties */

	ups->stroke_active = true;
	ups->size_pressure_value = stroke->cached_size_pressure;

	ups->pixel_radius = BKE_brush_size_get(scene, brush);

	if (BKE_brush_use_size_pressure(scene, brush) && paint_supports_dynamic_size(brush, mode)) {
		ups->pixel_radius *= stroke->cached_size_pressure;
	}

	if (paint_supports_dynamic_tex_coords(brush, mode)) {
		if (((brush->mtex.brush_map_mode == MTEX_MAP_MODE_VIEW) ||
		     (brush->mtex.brush_map_mode == MTEX_MAP_MODE_AREA) ||
		     (brush->mtex.brush_map_mode == MTEX_MAP_MODE_RANDOM)) &&
		    !(brush->flag & BRUSH_RAKE))
		{
			if (brush->flag & BRUSH_RANDOM_ROTATION)
				ups->brush_rotation = 2.0f * (float)M_PI * BLI_frand();
			else
				ups->brush_rotation = 0.0f;
		}

		if (brush->mtex.brush_map_mode == MTEX_MAP_MODE_RANDOM)
			BKE_brush_randomize_texture_coordinates(ups, false);
		else {
			copy_v2_v2(ups->tex_mouse, mouse);
		}

		/* take care of mask texture, if any */
		if (brush->mask_mtex.tex) {
			if (brush->mask_mtex.brush_map_mode == MTEX_MAP_MODE_RANDOM)
				BKE_brush_randomize_texture_coordinates(ups, true);
			else {
				copy_v2_v2(ups->mask_tex_mouse, mouse);
			}
		}
	}


	if (brush->flag & BRUSH_ANCHORED) {
		bool hit = false;
		float halfway[2];

		const float dx = mouse[0] - stroke->initial_mouse[0];
		const float dy = mouse[1] - stroke->initial_mouse[1];

		ups->anchored_size = ups->pixel_radius = sqrt(dx * dx + dy * dy);

		ups->brush_rotation = atan2(dx, dy) + M_PI;

		if (brush->flag & BRUSH_EDGE_TO_EDGE) {
			halfway[0] = dx * 0.5f + stroke->initial_mouse[0];
			halfway[1] = dy * 0.5f + stroke->initial_mouse[1];

			if (stroke->get_location) {
				if (stroke->get_location(C, location, halfway)) {
					hit = true;
					location_sampled = true;
					location_success = true;
				}
				else if (!paint_tool_require_location(brush, mode)) {
					hit = true;
				}
			}
			else {
				hit = true;
			}
		}
		if (hit) {
			copy_v2_v2(ups->anchored_initial_mouse, halfway);
			copy_v2_v2(ups->tex_mouse, halfway);
			copy_v2_v2(ups->mask_tex_mouse, halfway);
			copy_v2_v2(mouse, halfway);
			ups->anchored_size /= 2.0f;
			ups->pixel_radius  /= 2.0f;
			stroke->stroke_distance = ups->pixel_radius;
		}
		else {
			copy_v2_v2(ups->anchored_initial_mouse, stroke->initial_mouse);
			copy_v2_v2(mouse, stroke->initial_mouse);
			stroke->stroke_distance = ups->pixel_radius;
		}
		ups->pixel_radius /= stroke->zoom_2d;
		ups->draw_anchored = true;
	}
	else if (brush->flag & BRUSH_RAKE) {
		/* here we are using the initial mouse coordinate because we do not want the rake
		 * result to depend on jittering */
		if (!stroke->brush_init)
			copy_v2_v2(ups->last_rake, mouse_init);
		else
			paint_calculate_rake_rotation(ups, mouse_init);
	}

	if (!location_sampled) {
		if (stroke->get_location) {
			if (stroke->get_location(C, location, mouse))
				location_success = true;
			else if (!paint_tool_require_location(brush, mode))
				location_success = true;
		}
		else {
			zero_v3(location);
			location_success = true;
		}
	}

	return location_success;
}

static bool paint_stroke_use_jitter(PaintMode mode, Brush *brush, bool invert)
{
	bool use_jitter = (brush->flag & BRUSH_ABSOLUTE_JITTER) ?
		(brush->jitter_absolute != 0) : (brush->jitter != 0);

	/* jitter-ed brush gives weird and unpredictable result for this
	 * kinds of stroke, so manually disable jitter usage (sergey) */
	use_jitter &= (brush->flag & (BRUSH_DRAG_DOT | BRUSH_ANCHORED)) == 0;
	use_jitter &= (!ELEM(mode, PAINT_TEXTURE_2D, PAINT_TEXTURE_PROJECTIVE) ||
	               !(invert && brush->imagepaint_tool == PAINT_TOOL_CLONE));


	return use_jitter;
}

/* Put the location of the next stroke dot into the stroke RNA and apply it to the mesh */
static void paint_brush_stroke_add_step(bContext *C, wmOperator *op, const float mouse_in[2], float pressure)
{
	Scene *scene = CTX_data_scene(C);
	Paint *paint = BKE_paint_get_active_from_context(C);
	PaintMode mode = BKE_paintmode_get_active_from_context(C);
	Brush *brush = BKE_paint_brush(paint);
	PaintStroke *stroke = op->customdata;
	UnifiedPaintSettings *ups = stroke->ups;
	float mouse_out[2];
	PointerRNA itemptr;
	float location[3];

/* the following code is adapted from texture paint. It may not be needed but leaving here
 * just in case for reference (code in texpaint removed as part of refactoring).
 * It's strange that only texpaint had these guards. */
#if 0
	/* special exception here for too high pressure values on first touch in
	 * windows for some tablets, then we just skip first touch ..  */
	if (tablet && (pressure >= 0.99f) && ((pop->s.brush->flag & BRUSH_SPACING_PRESSURE) || BKE_brush_use_alpha_pressure(scene, pop->s.brush) || BKE_brush_use_size_pressure(scene, pop->s.brush)))
		return;

	/* This can be removed once fixed properly in
	 * BKE_brush_painter_paint(BrushPainter *painter, BrushFunc func, float *pos, double time, float pressure, void *user)
	 * at zero pressure we should do nothing 1/2^12 is 0.0002 which is the sensitivity of the most sensitive pen tablet available */
	if (tablet && (pressure < 0.0002f) && ((pop->s.brush->flag & BRUSH_SPACING_PRESSURE) || BKE_brush_use_alpha_pressure(scene, pop->s.brush) || BKE_brush_use_size_pressure(scene, pop->s.brush)))
		return;
#endif

	/* copy last position -before- jittering, or space fill code
	 * will create too many dabs */
	copy_v2_v2(stroke->last_mouse_position, mouse_in);
	stroke->last_pressure = pressure;

	if (paint_stroke_use_jitter(mode, brush, stroke->stroke_mode == BRUSH_STROKE_INVERT)) {
		float delta[2];
		float factor = stroke->zoom_2d;

		if (brush->flag & BRUSH_JITTER_PRESSURE)
			factor *= pressure;

		BKE_brush_jitter_pos(scene, brush, mouse_in, mouse_out);

		/* XXX: meh, this is round about because
		 * BKE_brush_jitter_pos isn't written in the best way to
		 * be reused here */
		if (factor != 1.0f) {
			sub_v2_v2v2(delta, mouse_out, mouse_in);
			mul_v2_fl(delta, factor);
			add_v2_v2v2(mouse_out, mouse_in, delta);
		}
	}
	else {
		copy_v2_v2(mouse_out, mouse_in);
	}

	if (!paint_brush_update(C, brush, mode, stroke, mouse_in, mouse_out, pressure, location)) {
		return;
	}

	/* Add to stroke */
	RNA_collection_add(op->ptr, "stroke", &itemptr);
	RNA_float_set(&itemptr, "size", ups->pixel_radius);
	RNA_float_set_array(&itemptr, "location", location);
	RNA_float_set_array(&itemptr, "mouse", mouse_out);
	RNA_boolean_set(&itemptr, "pen_flip", stroke->pen_flip);
	RNA_float_set(&itemptr, "pressure", pressure);

	stroke->update_step(C, stroke, &itemptr);

	/* don't record this for now, it takes up a lot of memory when doing long
	 * strokes with small brush size, and operators have register disabled */
	RNA_collection_clear(op->ptr, "stroke");
}

/* Returns zero if no sculpt changes should be made, non-zero otherwise */
static int paint_smooth_stroke(PaintStroke *stroke, float output[2], float *outpressure,
                               const PaintSample *sample, PaintMode mode)
{
	if (paint_supports_smooth_stroke(stroke->brush, mode)) {
		float radius = stroke->brush->smooth_stroke_radius * stroke->zoom_2d;
		float u = stroke->brush->smooth_stroke_factor, v = 1.0f - u;
		float dx = stroke->last_mouse_position[0] - sample->mouse[0];
		float dy = stroke->last_mouse_position[1] - sample->mouse[1];

		/* If the mouse is moving within the radius of the last move,
		 * don't update the mouse position. This allows sharp turns. */
		if (dx * dx + dy * dy <  radius * radius)
			return 0;

		output[0] = sample->mouse[0] * v + stroke->last_mouse_position[0] * u;
		output[1] = sample->mouse[1] * v + stroke->last_mouse_position[1] * u;
		*outpressure = sample->pressure * v + stroke->last_pressure * u;
	}
	else {
		output[0] = sample->mouse[0];
		output[1] = sample->mouse[1];
		*outpressure = sample->pressure;
	}

	return 1;
}

static float paint_space_stroke_spacing(const Scene *scene, PaintStroke *stroke, float size_pressure, float spacing_pressure)
{
	/* brushes can have a minimum size of 1.0 but with pressure it can be smaller then a pixel
	 * causing very high step sizes, hanging blender [#32381] */
	const float size_clamp = max_ff(1.0f, BKE_brush_size_get(scene, stroke->brush) * size_pressure);
	float spacing = stroke->brush->spacing;

	/* apply spacing pressure */
	if (stroke->brush->flag & BRUSH_SPACING_PRESSURE)
		spacing = spacing * (1.5f - spacing_pressure);

	/* stroke system is used for 2d paint too, so we need to account for
	 * the fact that brush can be scaled there. */
	spacing *= stroke->zoom_2d;

	return max_ff(1.0, size_clamp * spacing / 50.0f);
}



static float paint_stroke_overlapped_curve(Brush *br, float x, float spacing)
{
	int i;
	const int n = 100 / spacing;
	const float h = spacing / 50.0f;
	const float x0 = x - 1;

	float sum;

	sum = 0;
	for (i = 0; i < n; i++) {
		float xx;

		xx = fabs(x0 + i * h);

		if (xx < 1.0f)
			sum += BKE_brush_curve_strength(br, xx, 1);
	}

	return sum;
}

static float paint_stroke_integrate_overlap(Brush *br, float factor)
{
	int i;
	int m;
	float g;
	float max;

	float spacing = br->spacing * factor;

	if (!(br->flag & BRUSH_SPACE_ATTEN && (br->spacing < 100)))
		return 1.0;

	m = 10;
	g = 1.0f / m;
	max = 0;
	for (i = 0; i < m; i++) {
		float overlap = paint_stroke_overlapped_curve(br, i * g, spacing);

		if (overlap > max)
			max = overlap;
	}

	return 1.0f / max;
}

static float paint_space_stroke_spacing_variable(const Scene *scene, PaintStroke *stroke, float pressure, float dpressure, float length)
{
	if (BKE_brush_use_size_pressure(scene, stroke->brush)) {
		/* use pressure to modify size. set spacing so that at 100%, the circles
		 * are aligned nicely with no overlap. for this the spacing needs to be
		 * the average of the previous and next size. */
		float s = paint_space_stroke_spacing(scene, stroke, 1.0f, pressure);
		float q = s * dpressure / (2.0f * length);
		float pressure_fac = (1.0f + q) / (1.0f - q);

		float last_size_pressure = stroke->last_pressure;
		float new_size_pressure = stroke->last_pressure * pressure_fac;

		/* average spacing */
		float last_spacing = paint_space_stroke_spacing(scene, stroke, last_size_pressure, pressure);
		float new_spacing = paint_space_stroke_spacing(scene, stroke, new_size_pressure, pressure);

		return 0.5f * (last_spacing + new_spacing);
	}
	else {
		/* no size pressure */
		return paint_space_stroke_spacing(scene, stroke, 1.0f, pressure);
	}
}

/* For brushes with stroke spacing enabled, moves mouse in steps
 * towards the final mouse location. */
static int paint_space_stroke(bContext *C, wmOperator *op, const float final_mouse[2], float final_pressure)
{
	const Scene *scene = CTX_data_scene(C);
	PaintStroke *stroke = op->customdata;
	UnifiedPaintSettings *ups = stroke->ups;
	int cnt = 0;

	float pressure, dpressure;
	float mouse[2], dmouse[2];
	float length;
	float no_pressure_spacing = paint_space_stroke_spacing(scene, stroke, 1.0f, 1.0f);

	sub_v2_v2v2(dmouse, final_mouse, stroke->last_mouse_position);

	pressure = stroke->last_pressure;
	dpressure = final_pressure - stroke->last_pressure;

	length = normalize_v2(dmouse);

	while (length > 0.0f) {
		float spacing = paint_space_stroke_spacing_variable(scene, stroke, pressure, dpressure, length);

		if (length >= spacing) {
			mouse[0] = stroke->last_mouse_position[0] + dmouse[0] * spacing;
			mouse[1] = stroke->last_mouse_position[1] + dmouse[1] * spacing;
			pressure = stroke->last_pressure + (spacing / length) * dpressure;

			ups->overlap_factor = paint_stroke_integrate_overlap(stroke->brush, spacing / no_pressure_spacing);

			stroke->stroke_distance += spacing / stroke->zoom_2d;
			paint_brush_stroke_add_step(C, op, mouse, pressure);

			length -= spacing;
			pressure = stroke->last_pressure;
			dpressure = final_pressure - stroke->last_pressure;

			cnt++;
		}
		else {
			break;
		}
	}

	return cnt;
}

/**** Public API ****/

PaintStroke *paint_stroke_new(bContext *C,
                              wmOperator *op,
                              StrokeGetLocation get_location,
                              StrokeTestStart test_start,
                              StrokeUpdateStep update_step,
                              StrokeRedraw redraw,
                              StrokeDone done, int event_type)
{
	PaintStroke *stroke = MEM_callocN(sizeof(PaintStroke), "PaintStroke");
	ToolSettings *toolsettings = CTX_data_tool_settings(C);
	UnifiedPaintSettings *ups = &toolsettings->unified_paint_settings;
	Brush *br = stroke->brush = BKE_paint_brush(BKE_paint_get_active_from_context(C));
	float zoomx, zoomy;

	view3d_set_viewcontext(C, &stroke->vc);
	if (stroke->vc.v3d)
		view3d_get_transformation(stroke->vc.ar, stroke->vc.rv3d, stroke->vc.obact, &stroke->mats);

	stroke->get_location = get_location;
	stroke->test_start = test_start;
	stroke->update_step = update_step;
	stroke->redraw = redraw;
	stroke->done = done;
	stroke->event_type = event_type; /* for modal, return event */
	stroke->ups = ups;
	stroke->stroke_mode = RNA_enum_get(op->ptr, "mode");

	get_imapaint_zoom(C, &zoomx, &zoomy);
	stroke->zoom_2d = max_ff(zoomx, zoomy);

	if ((br->flag & BRUSH_CURVE) &&
	    RNA_struct_property_is_set(op->ptr, "mode"))
	{
		RNA_enum_set(op->ptr, "mode", BRUSH_STROKE_NORMAL);
	}
	/* initialize here */
	ups->overlap_factor = 1.0;
	ups->stroke_active = true;

	/* initialize here to avoid initialization conflict with threaded strokes */
	curvemapping_initialize(br->curve);
	
	BKE_paint_set_overlay_override(br->overlay_flags);

	return stroke;
}

void paint_stroke_data_free(struct wmOperator *op)
{
	BKE_paint_set_overlay_override(0);
	MEM_SAFE_FREE(op->customdata);
}

static void stroke_done(struct bContext *C, struct wmOperator *op)
{
	struct PaintStroke *stroke = op->customdata;
	UnifiedPaintSettings *ups = stroke->ups;

	ups->draw_anchored = false;
	ups->stroke_active = false;

	/* reset rotation here to avoid doing so in cursor display */
	if (!(stroke->brush->flag & BRUSH_RAKE))
		ups->brush_rotation = 0.0f;

	if (stroke->stroke_started) {
		if (stroke->redraw)
			stroke->redraw(C, stroke, true);

		if (stroke->done)
			stroke->done(C, stroke);
	}

	if (stroke->timer) {
		WM_event_remove_timer(
			CTX_wm_manager(C),
			CTX_wm_window(C),
			stroke->timer);
	}

	if (stroke->stroke_cursor)
		WM_paint_cursor_end(CTX_wm_manager(C), stroke->stroke_cursor);

	BLI_freelistN(&stroke->line);

	paint_stroke_data_free(op);
}

/* Returns zero if the stroke dots should not be spaced, non-zero otherwise */
bool paint_space_stroke_enabled(Brush *br, PaintMode mode)
{
	return (br->flag & BRUSH_SPACE) && paint_supports_dynamic_size(br, mode);
}

static bool sculpt_is_grab_tool(Brush *br)
{
	return ELEM(br->sculpt_tool,
	             SCULPT_TOOL_GRAB,
	             SCULPT_TOOL_THUMB,
	             SCULPT_TOOL_ROTATE,
	             SCULPT_TOOL_SNAKE_HOOK);
}

/* return true if the brush size can change during paint (normally used for pressure) */
bool paint_supports_dynamic_size(Brush *br, PaintMode mode)
{
	if (br->flag & BRUSH_ANCHORED)
		return false;

	switch (mode) {
		case PAINT_SCULPT:
			if (sculpt_is_grab_tool(br))
				return false;
			break;

		case PAINT_TEXTURE_2D: /* fall through */
		case PAINT_TEXTURE_PROJECTIVE:
			if ((br->imagepaint_tool == PAINT_TOOL_FILL) &&
			    (br->flag & BRUSH_USE_GRADIENT))
			{
				return false;
			}
			break;

		default:
			break;
	}
	return true;
}

bool paint_supports_smooth_stroke(Brush *br, PaintMode mode)
{
	if (!(br->flag & BRUSH_SMOOTH_STROKE) ||
	    (br->flag & (BRUSH_ANCHORED | BRUSH_DRAG_DOT | BRUSH_LINE)))
	{
		return false;
	}

	switch (mode) {
		case PAINT_SCULPT:
			if (sculpt_is_grab_tool(br))
				return false;
			break;
		default:
			break;
	}
	return true;
}

bool paint_supports_texture(PaintMode mode)
{
	/* ommit: PAINT_WEIGHT, PAINT_SCULPT_UV, PAINT_INVALID */
	return ELEM(mode, PAINT_SCULPT, PAINT_VERTEX, PAINT_TEXTURE_PROJECTIVE, PAINT_TEXTURE_2D);
}

/* return true if the brush size can change during paint (normally used for pressure) */
bool paint_supports_dynamic_tex_coords(Brush *br, PaintMode mode)
{
	if (br->flag & BRUSH_ANCHORED)
		return false;

	switch (mode) {
		case PAINT_SCULPT:
			if (sculpt_is_grab_tool(br))
				return false;
			break;
		default:
			break;
	}
	return true;
}

#define PAINT_STROKE_MODAL_CANCEL 1

/* called in paint_ops.c, on each regeneration of keymaps  */
struct wmKeyMap *paint_stroke_modal_keymap(struct wmKeyConfig *keyconf)
{
	static struct EnumPropertyItem modal_items[] = {
		{PAINT_STROKE_MODAL_CANCEL, "CANCEL", 0,
		"Cancel",
		"Cancel and undo a stroke in progress"},

		{ 0 }
	};

	static const char *name = "Paint Stroke Modal";

	struct wmKeyMap *keymap = WM_modalkeymap_get(keyconf, name);

	/* this function is called for each spacetype, only needs to add map once */
	if (!keymap) {
		keymap = WM_modalkeymap_add(keyconf, name, modal_items);

		/* items for modal map */
		WM_modalkeymap_add_item(
			keymap, ESCKEY, KM_PRESS, KM_ANY, 0, PAINT_STROKE_MODAL_CANCEL);
	}

	return keymap;
}

static void paint_stroke_add_sample(const Paint *paint,
                                    PaintStroke *stroke,
                                    float x, float y, float pressure)
{
	PaintSample *sample = &stroke->samples[stroke->cur_sample];
	int max_samples = MIN2(PAINT_MAX_INPUT_SAMPLES,
	                       MAX2(paint->num_input_samples, 1));

	sample->mouse[0] = x;
	sample->mouse[1] = y;
	sample->pressure = pressure;

	stroke->cur_sample++;
	if (stroke->cur_sample >= max_samples)
		stroke->cur_sample = 0;
	if (stroke->num_samples < max_samples)
		stroke->num_samples++;
}

static void paint_stroke_sample_average(const PaintStroke *stroke,
                                        PaintSample *average)
{
	int i;
	
	memset(average, 0, sizeof(*average));

	BLI_assert(stroke->num_samples > 0);
	
	for (i = 0; i < stroke->num_samples; i++) {
		add_v2_v2(average->mouse, stroke->samples[i].mouse);
		average->pressure += stroke->samples[i].pressure;
	}

	mul_v2_fl(average->mouse, 1.0f / stroke->num_samples);
	average->pressure /= stroke->num_samples;

	/*printf("avg=(%f, %f), num=%d\n", average->mouse[0], average->mouse[1], stroke->num_samples);*/
}

/**
 * Slightly different version of spacing for line/curve strokes,
 * makes sure the dabs stay on the line path.
 */
static void paint_line_strokes_spacing(
        bContext *C, wmOperator *op, PaintStroke *stroke, float spacing, float *length_residue,
        const float old_pos[2], const float new_pos[2])
{
	UnifiedPaintSettings *ups = stroke->ups;

	float mouse[2], dmouse[2];
	float length;

	sub_v2_v2v2(dmouse, new_pos, old_pos);
	copy_v2_v2(stroke->last_mouse_position, old_pos);

	length = normalize_v2(dmouse);

	BLI_assert(length >= 0.0f);

	if (length == 0.0f)
		return;

	while (length > 0.0f) {
		float spacing_final = spacing - *length_residue;
		length += *length_residue;
		*length_residue = 0.0;

		if (length >= spacing) {
			mouse[0] = stroke->last_mouse_position[0] + dmouse[0] * spacing_final;
			mouse[1] = stroke->last_mouse_position[1] + dmouse[1] * spacing_final;

			ups->overlap_factor = paint_stroke_integrate_overlap(stroke->brush, 1.0);

			stroke->stroke_distance += spacing / stroke->zoom_2d;
			paint_brush_stroke_add_step(C, op, mouse, 1.0);

			length -= spacing;
			spacing_final = spacing;
		}
		else {
			break;
		}
	}

	*length_residue = length;
}


static void paint_stroke_line_end(bContext *C, wmOperator *op, PaintStroke *stroke, float mouse[2])
{
	Brush *br = stroke->brush;
	if (stroke->stroke_started && (br->flag & BRUSH_LINE)) {
		stroke->ups->overlap_factor = paint_stroke_integrate_overlap(br, 1.0);

		paint_brush_stroke_add_step(C, op, stroke->last_mouse_position, 1.0);
		paint_space_stroke(C, op, mouse, 1.0);
	}
}

static bool paint_stroke_curve_end(bContext *C, wmOperator *op, PaintStroke *stroke)
{
	Brush *br = stroke->brush;
	if (br->flag & BRUSH_CURVE) {
		const Scene *scene = CTX_data_scene(C);
		const float spacing = paint_space_stroke_spacing(scene, stroke, 1.0f, 1.0f);
		PaintCurve *pc = br->paint_curve;
		PaintCurvePoint *pcp;
		float length_residue = 0.0f;
		int i;

		if (!pc)
			return true;

		pcp = pc->points;
		stroke->ups->overlap_factor = paint_stroke_integrate_overlap(br, 1.0);

		for (i = 0; i < pc->tot_points - 1; i++, pcp++) {
			int j;
			float data[(PAINT_CURVE_NUM_SEGMENTS + 1) * 2];
			PaintCurvePoint *pcp_next = pcp + 1;

			for (j = 0; j < 2; j++)
				BKE_curve_forward_diff_bezier(
				        pcp->bez.vec[1][j],
				        pcp->bez.vec[2][j],
				        pcp_next->bez.vec[0][j],
				        pcp_next->bez.vec[1][j],
				        data + j, PAINT_CURVE_NUM_SEGMENTS, sizeof(float[2]));


			for (j = 0; j < PAINT_CURVE_NUM_SEGMENTS; j++) {
				if (!stroke->stroke_started) {
					stroke->last_pressure = 1.0;
					copy_v2_v2(stroke->last_mouse_position, data + 2 * j);
					stroke->stroke_started = stroke->test_start(C, op, stroke->last_mouse_position);

					if (stroke->stroke_started) {
						paint_brush_stroke_add_step(C, op, data + 2 * j, 1.0);
						paint_line_strokes_spacing(C, op, stroke, spacing, &length_residue, data + 2 * j, data + 2 * (j + 1));
					}
				}
				else {
					paint_line_strokes_spacing(C, op, stroke, spacing, &length_residue, data + 2 * j, data + 2 * (j + 1));
				}
			}
		}

		stroke_done(C, op);
		return true;
	}

	return false;
}


int paint_stroke_modal(bContext *C, wmOperator *op, const wmEvent *event)
{
	Paint *p = BKE_paint_get_active_from_context(C);
	PaintMode mode = BKE_paintmode_get_active_from_context(C);
	PaintStroke *stroke = op->customdata;
	Brush *br = stroke->brush;
	PaintSample sample_average;
	float mouse[2];
	bool first_dab = false;
	bool first_modal = false;
	bool redraw = false;
	float pressure;

	/* see if tablet affects event. Line, anchored and drag dot strokes do not support pressure */
	pressure = (br->flag & (BRUSH_LINE | BRUSH_ANCHORED | BRUSH_DRAG_DOT)) ? 1.0f : WM_event_tablet_data(event, &stroke->pen_flip, NULL);

	paint_stroke_add_sample(p, stroke, event->mval[0], event->mval[1], pressure);
	paint_stroke_sample_average(stroke, &sample_average);

	/* let NDOF motion pass through to the 3D view so we can paint and rotate simultaneously!
	 * this isn't perfect... even when an extra MOUSEMOVE is spoofed, the stroke discards it
	 * since the 2D deltas are zero -- code in this file needs to be updated to use the
	 * post-NDOF_MOTION MOUSEMOVE */
	if (event->type == NDOF_MOTION)
		return OPERATOR_PASS_THROUGH;

	/* one time initialization */
	if (!stroke->stroke_init) {
		if (paint_stroke_curve_end(C, op, stroke))
			return OPERATOR_FINISHED;

		if (paint_supports_smooth_stroke(br, mode))
			stroke->stroke_cursor =
			    WM_paint_cursor_activate(CTX_wm_manager(C), paint_poll, paint_draw_smooth_cursor, stroke);

		stroke->stroke_init = true;
		first_modal = true;
	}

	/* one time stroke initialization */
	if (!stroke->stroke_started) {
		stroke->last_pressure = sample_average.pressure;
		copy_v2_v2(stroke->last_mouse_position, sample_average.mouse);
		stroke->stroke_started = stroke->test_start(C, op, sample_average.mouse);
		BLI_assert((stroke->stroke_started & ~1) == 0);  /* 0/1 */

		if (stroke->stroke_started) {
			if (br->flag & BRUSH_AIRBRUSH)
				stroke->timer = WM_event_add_timer(CTX_wm_manager(C), CTX_wm_window(C), TIMER, stroke->brush->rate);

			if (br->flag & BRUSH_LINE) {
				stroke->stroke_cursor =
					WM_paint_cursor_activate(CTX_wm_manager(C), paint_poll, paint_draw_line_cursor, stroke);
			}

			first_dab = true;
		}
	}

	/* Cancel */
	if (event->type == EVT_MODAL_MAP && event->val == PAINT_STROKE_MODAL_CANCEL) {
		if (op->type->cancel) {
			op->type->cancel(C, op);
		}
		else {
			paint_stroke_cancel(C, op);
		}
		return OPERATOR_CANCELLED;
	}

	if (event->type == stroke->event_type && !first_modal) {
		if (event->val == KM_RELEASE) {
			paint_stroke_line_end (C, op, stroke, sample_average.mouse);
			stroke_done(C, op);
			return OPERATOR_FINISHED;
		}
	}
	else if (ELEM(event->type, RETKEY, SPACEKEY)) {
		paint_stroke_line_end(C, op, stroke, sample_average.mouse);
		stroke_done(C, op);
		return OPERATOR_FINISHED;
	}
	else if ((br->flag & BRUSH_LINE) && stroke->stroke_started &&
	         (first_modal || (ELEM(event->type, MOUSEMOVE, INBETWEEN_MOUSEMOVE))))
	{
		if (br->flag & BRUSH_RAKE) {
			copy_v2_v2(stroke->ups->last_rake, stroke->last_mouse_position);
			paint_calculate_rake_rotation(stroke->ups,  sample_average.mouse);
		}
	}
	else if (first_modal ||
	         /* regular dabs */
	         (!(br->flag & (BRUSH_AIRBRUSH)) && (ELEM(event->type, MOUSEMOVE, INBETWEEN_MOUSEMOVE))) ||
	         /* airbrush */
	         ((br->flag & BRUSH_AIRBRUSH) && event->type == TIMER && event->customdata == stroke->timer))
	{
		if (paint_smooth_stroke(stroke, mouse, &pressure, &sample_average, mode)) {
			if (stroke->stroke_started) {
				if (paint_space_stroke_enabled(br, mode)) {
					if (paint_space_stroke(C, op, mouse, pressure))
						redraw = true;
				}
				else {
					float dmouse[2];
					sub_v2_v2v2(dmouse, mouse, stroke->last_mouse_position);
					stroke->stroke_distance += len_v2(dmouse);
					paint_brush_stroke_add_step(C, op, mouse, pressure);
					redraw = true;
				}
			}
		}
	}

	/* we want the stroke to have the first daub at the start location
	 * instead of waiting till we have moved the space distance */
	if (first_dab &&
	    paint_space_stroke_enabled(br, mode) &&
	    !(br->flag & BRUSH_SMOOTH_STROKE))
	{
		stroke->ups->overlap_factor = paint_stroke_integrate_overlap(br, 1.0);
		paint_brush_stroke_add_step(C, op, sample_average.mouse, sample_average.pressure);
		redraw = true;
	}

	/* do updates for redraw. if event is inbetween mousemove there are more
	 * coming, so postpone potentially slow redraw updates until all are done */
	if (event->type != INBETWEEN_MOUSEMOVE) {
		wmWindow *window = CTX_wm_window(C);
		ARegion *ar = CTX_wm_region(C);

		/* At the very least, invalidate the cursor */
		if (ar && (p->flags & PAINT_SHOW_BRUSH))
			WM_paint_cursor_tag_redraw(window, ar);

		if (redraw && stroke->redraw)
			stroke->redraw(C, stroke, false);
	}

	return OPERATOR_RUNNING_MODAL;
}

int paint_stroke_exec(bContext *C, wmOperator *op)
{
	PaintStroke *stroke = op->customdata;

	/* only when executed for the first time */
	if (stroke->stroke_started == 0) {
		/* XXX stroke->last_mouse_position is unset, this may cause problems */
		stroke->test_start(C, op, NULL);
		stroke->stroke_started = 1;
	}

	RNA_BEGIN (op->ptr, itemptr, "stroke")
	{
		stroke->update_step(C, stroke, &itemptr);
	}
	RNA_END;

	stroke_done(C, op);

	return OPERATOR_FINISHED;
}

void paint_stroke_cancel(bContext *C, wmOperator *op)
{
	stroke_done(C, op);
}

ViewContext *paint_stroke_view_context(PaintStroke *stroke)
{
	return &stroke->vc;
}

void *paint_stroke_mode_data(struct PaintStroke *stroke)
{
	return stroke->mode_data;
}

bool paint_stroke_flipped(struct PaintStroke *stroke)
{
	return stroke->pen_flip;
}

float paint_stroke_distance_get(struct PaintStroke *stroke)
{
	return stroke->stroke_distance;
}

void paint_stroke_set_mode_data(PaintStroke *stroke, void *mode_data)
{
	stroke->mode_data = mode_data;
}

int paint_poll(bContext *C)
{
	Paint *p = BKE_paint_get_active_from_context(C);
	Object *ob = CTX_data_active_object(C);
	ScrArea *sa = CTX_wm_area(C);
	ARegion *ar = CTX_wm_region(C);

	return p && ob && BKE_paint_brush(p) &&
	       (sa && ELEM(sa->spacetype, SPACE_VIEW3D, SPACE_IMAGE)) &&
	       (ar && ar->regiontype == RGN_TYPE_WINDOW);
}
