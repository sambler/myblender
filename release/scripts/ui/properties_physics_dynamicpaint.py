# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>
import bpy


class PhysicButtonsPanel():
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "physics"

    @classmethod
    def poll(cls, context):
        ob = context.object
        rd = context.scene.render
        return (ob and ob.type == 'MESH') and (not rd.use_game_engine)


class PHYSICS_PT_mesh_paint(PhysicButtonsPanel, bpy.types.Panel):
    bl_label = "Dynamic Paint"

    def draw(self, context):
        layout = self.layout

        md = context.dynamic_paint
        ob = context.object

        split = layout.split()

        if md:
            # remove modifier + settings
            split.context_pointer_set("modifier", md)
            split.operator("object.modifier_remove", text="Remove")

            row = split.row(align=True)
            row.prop(md, "show_render", text="")
            row.prop(md, "show_viewport", text="")

        else:
            # add modifier
            split.operator("object.modifier_add", text="Add").type = 'DYNAMIC_PAINT'
            split.label()

        if md:
            layout.prop(md, "dynamicpaint_type", expand=True)

            if md.dynamicpaint_type == 'CANVAS':

                canvas = md.canvas_settings


                layout.label(text="Settings:")
                layout.prop(canvas, "resolution", text="Resolution")
                layout.prop(canvas, "antialias", text="Anti-aliasing")

                layout.label(text="Frames:")
                layout.prop(canvas, "substeps", text="Sub-Steps")
                split = layout.split()
                col = split.column()
                col.prop(canvas, "start_frame", text="Start Frame")
                col = split.column()
                col.prop(canvas, "end_frame", text="End Frame")

                layout.label(text="Bake:")
                layout.operator("dpaint.bake", text="Bake Dynamic Paint", icon='MOD_FLUIDSIM')


            elif md.dynamicpaint_type == 'PAINT':
                paint = md.paint_settings
                
                layout.row().label(text="Paint:")
                split = layout.split()
                col = split.column()
                col.prop(paint, "do_paint", text="Affect Paint:")
                sub = col.column()
                sub.active = paint.do_paint
                sub.prop(paint, "use_material", text="Use Material")
                sub.prop(paint, "abs_alpha", text="Absolute Alpha")
                sub.prop(paint, "paint_erase", text="Erase Paint")
                subs = split.column()
                psub = subs.column()
                psub.active = (not paint.use_material) and paint.do_paint
                psub.prop(paint, "paint_color", text="Paint Color")
                psub.prop(paint, "paint_alpha", text="Alpha")
                psub= subs.column()
                psub.active = paint.do_paint
                psub.prop(paint, "paint_wetness", text="Wetness")
                
                layout.row().label(text="Displace:")
                split = layout.split()
                sub = split.column()
                sub.prop(paint, "do_displace", text="Affect Displace:")
                # try if curve display bug is fixed
                # sub.template_curve_mapping(paint, "disp_curve")


class PHYSICS_PT_mp_paintmap(PhysicButtonsPanel, bpy.types.Panel):
    bl_label = "Dynamic Paint: Output"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        md = context.dynamic_paint
        return md and (md.dynamicpaint_type == 'CANVAS')

    def draw(self, context):
        layout = self.layout

        canvas = context.dynamic_paint.canvas_settings

        layout.row().label(text="Paintmaps:")
        split = layout.split()
        col = split.column()
        col.prop(canvas, "out_paint", text="Output Paint")
        sub = split.column()
        sub.active = canvas.out_paint
        sub.prop(canvas, "premultiply", text="Premultiply alpha")

        col = layout.row().column()
        col.active = canvas.out_paint
        col.prop(canvas, "paint_output_path", text="")

        layout.row().label(text="Wetmaps:")
        split = layout.split()
        col = split.column()
        col.prop(canvas, "out_wet", text="Ouput Wetmaps:")
        sub = col.column()
        sub.active = canvas.out_wet
        sub.prop(canvas, "wet_output_path", text="")


        layout.row().label(text="Displacement:")
        split = layout.split()
        col = split.column()
        col.prop(canvas, "out_disp", text="Output Displacement")
        sub = col.column()
        sub.active = canvas.out_disp

        sub.prop(canvas, "displace_output_path", text="")
        sub.prop(canvas, "displacement", text="Strength")
        split = layout.row().split()
        psub = split.column()
        psub.active = canvas.out_disp
        psub.prop(canvas, "disp_type", text="Type")
        psub = split.column()
        psub.active = canvas.out_disp
        psub.prop(canvas, "disp_format", text="Format")


class PHYSICS_PT_mp_wetmap(PhysicButtonsPanel, bpy.types.Panel):
    bl_label = "Dynamic Paint: Advanced"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        md = context.dynamic_paint
        return md and (md.dynamicpaint_type == 'CANVAS')

    def draw(self, context):
        layout = self.layout

        canvas = context.dynamic_paint.canvas_settings
        ob = context.object

        col = layout.column()
        col.label(text="Paint Dry Speed:")
        col.prop(canvas, "dry_log", text="Slow")
        col.prop(canvas, "dry_speed", text="Dry Time")

        col = layout.column()
        col.prop(canvas, "dissolve_paint", text="Paint Fade:")
        sub = col.column()
        sub.active = canvas.dissolve_paint
        sub.prop(canvas, "dissolve_speed", text="Time")

        col = layout.column()
        col.prop(canvas, "flatten_disp", text="Flatten Displace:")
        sub = col.column()
        sub.active = canvas.flatten_disp
        sub.prop(canvas, "flatten_speed", text="Time")
		
        layout.row().label(text="UV:")
        layout.prop_search(canvas, "uv_layer", ob.data, "uv_textures")

        layout.row().label(text="Paint Effects:")
        layout.prop(canvas, "effect_ui", text="Effect")

        if canvas.effect_ui == "SPREAD":
            col = layout.row().column()
            col.prop(canvas, "do_spread", text="Do Spread Effect")
            sub = col.column()
            sub.active = canvas.do_spread
            sub.prop(canvas, "spread_speed", text="Spread Speed")

        if canvas.effect_ui == "DRIP":
            col = layout.row().column()
            col.prop(canvas, "do_drip", text="Do Drip Effect")
            sub = col.column()
            sub.active = canvas.do_drip
            sub.prop(canvas, "drip_speed", text="Drip Speed")

        if canvas.effect_ui == "SHRINK":
            col = layout.row().column()
            col.prop(canvas, "do_shrink", text="Do Shrink Effect")
            sub = col.column()
            sub.active = canvas.do_shrink
            sub.prop(canvas, "shrink_speed", text="Shrink Speed")


class PHYSICS_PT_mp_displacemap(PhysicButtonsPanel, bpy.types.Panel):
    bl_label = "Dynamic Paint: Advanced"
    bl_options = {'DEFAULT_CLOSED'}

    @classmethod
    def poll(cls, context):
        md = context.dynamic_paint
        return md and (md.dynamicpaint_type == 'PAINT')

    def draw(self, context):
        layout = self.layout

        paint = context.dynamic_paint.paint_settings
        ob = context.object
		
        split = layout.split()
        col = split.column()
        col.label(text="Source:")
        col.prop(paint, "paint_source", expand=False)

        if paint.paint_source == "PSYS":
            col.label(text="Particle System:")
            col.prop_search(paint, "psys", ob, "particle_systems", text="")
            if paint.psys:
                col.label(text="Particle effect:")
                sub = col.column()
                sub.active = not paint.use_part_radius
                sub.prop(paint, "solid_radius", text="Solid Radius")
                col.prop(paint, "use_part_radius", text="Use Particle's Radius")
                col.prop(paint, "smooth_radius", text="Smooth radius")

        if paint.paint_source == "DISTANCE" or paint.paint_source == "VOLDIST":
            col.prop(paint, "paint_distance", text="Paint Distance")
            split = layout.row().split()
            sub = split.column()
            sub.prop(paint, "prox_facealigned", text="Face Aligned")
            sub = split.column()
            sub.prop(paint, "prox_falloff", text="Falloff")
            if paint.prox_falloff == "RAMP":
                col = layout.row().column()
                col.label(text="Fallout Curve:")
                col.prop(paint, "prox_ramp_alpha", text="Only Use Alpha")
                col.template_color_ramp(paint, "paint_ramp", expand=True)


            # edge displace (dummy) panels, waiting due to curve display bug
            #sub = col.column()
            #sub.active = paint.do_displace
            #sub.prop(paint, "displace_distance", text="Displace Distance")
            #sub.prop(paint, "edge_displace", text="Edge Displace:")
            #sub = sub.column()
            #sub.active = paint.edge_displace
            #sub.prop(paint, "prox_displace_strength", text="Amount")
            #sub.template_curve_mapping(paint, "disp_curve")

def register():
    pass


def unregister():
    pass

if __name__ == "__main__":
    register()
