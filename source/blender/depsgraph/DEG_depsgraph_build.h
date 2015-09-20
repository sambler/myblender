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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2013 Blender Foundation.
 * All rights reserved.
 *
 * Original Author: Joshua Leung
 * Contributor(s): Lukas Toenne
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/depsgraph/DEG_depsgraph_build.h
 *  \ingroup depsgraph
 *
 * Public API for Depsgraph
 */

#ifndef __DEG_DEPSGRAPH_BUILD_H__
#define __DEG_DEPSGRAPH_BUILD_H__

/* ************************************************* */

/* Dependency Graph */
struct Depsgraph;

/* ------------------------------------------------ */

struct Main;
struct Scene;

struct PointerRNA;
struct PropertyRNA;

#ifdef __cplusplus
extern "C" {
#endif

/* Graph Building -------------------------------- */

/* Build depsgraph for the given scene, and dump results in given graph container */
void DEG_graph_build_from_scene(struct Depsgraph *graph, struct Main *bmain, struct Scene *scene);

/* Tag relations from the given graph for update. */
void DEG_graph_tag_relations_update(struct Depsgraph *graph);

/* Tag all relations in the database for update.*/
void DEG_relations_tag_update(struct Main *bmain);

/* Create new graph if didn't exist yet,
 * or update relations if graph was tagged for update.
 */
void DEG_scene_relations_update(struct Main *bmain, struct Scene *scene);

/* Rebuild dependency graph only for a given scene. */
void DEG_scene_relations_rebuild(struct Main *bmain,
                                 struct Scene *scene);

/* Delete scene graph. */
void DEG_scene_graph_free(struct Scene *scene);

/* Add Dependencies  ----------------------------- */

/* Handle for components to define their dependencies from callbacks.
 * This is generated by the depsgraph and passed to dependency callbacks
 * as a symbolic reference to the current DepsNode.
 * All relations will be defined in reference to that node.
 */
struct DepsNodeHandle;

struct Object;

typedef enum eDepsSceneComponentType {
	DEG_SCENE_COMP_PARAMETERS,     /* Parameters Component - Default when nothing else fits (i.e. just SDNA property setting) */
	DEG_SCENE_COMP_ANIMATION,      /* Animation Component */                 // XXX: merge in with parameters?
	DEG_SCENE_COMP_SEQUENCER,      /* Sequencer Component (Scene Only) */
} eDepsSceneComponentType;

typedef enum eDepsObjectComponentType {
	DEG_OB_COMP_PARAMETERS,        /* Parameters Component - Default when nothing else fits (i.e. just SDNA property setting) */
	DEG_OB_COMP_PROXY,             /* Generic "Proxy-Inherit" Component */   // XXX: Also for instancing of subgraphs?
	DEG_OB_COMP_ANIMATION,         /* Animation Component */                 // XXX: merge in with parameters?
	DEG_OB_COMP_TRANSFORM,         /* Transform Component (Parenting/Constraints) */
	DEG_OB_COMP_GEOMETRY,          /* Geometry Component (DerivedMesh/Displist) */
	
	/* Evaluation-Related Outer Types (with Subdata) */
	DEG_OB_COMP_EVAL_POSE,         /* Pose Component - Owner/Container of Bones Eval */
	DEG_OB_COMP_BONE,              /* Bone Component - Child/Subcomponent of Pose */
	
	DEG_OB_COMP_EVAL_PARTICLES,    /* Particle Systems Component */
	DEG_OB_COMP_SHADING,           /* Material Shading Component */
} eDepsObjectComponentType;

void DEG_add_scene_relation(struct DepsNodeHandle *node, struct Scene *scene, eDepsSceneComponentType component, const char *description);
void DEG_add_object_relation(struct DepsNodeHandle *node, struct Object *ob, eDepsObjectComponentType component, const char *description);
void DEG_add_bone_relation(struct DepsNodeHandle *handle, struct Object *ob, const char *bone_name, eDepsObjectComponentType component, const char *description);

/* TODO(sergey): Remove once all geometry update is granular. */
void DEG_add_special_eval_flag(struct Depsgraph *graph, struct ID *id, short flag);

/* ************************************************ */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  /* __DEG_DEPSGRAPH_BUILD_H__ */
