/**
 * $Id: LightExporter.h 32355 2010-10-06 20:40:16Z gsrb3d $
 *
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
 * Contributor(s): Chingiz Dyussenov, Arystanbek Dyussenov, Jan Diederich, Tod Liverseed,
 *                 Nathan Letwory
 *
 * ***** END GPL LICENSE BLOCK *****
 */

#include <map>

#include "COLLADASWEffectProfile.h"

#include "EffectExporter.h"
#include "MaterialExporter.h"

#include "DNA_mesh_types.h"
#include "DNA_texture_types.h"

#include "BKE_customdata.h"

#include "collada_internal.h"
#include "collada_utils.h"

// OB_MESH is assumed
static std::string getActiveUVLayerName(Object *ob)
{
	Mesh *me = (Mesh*)ob->data;

	int num_layers = CustomData_number_of_layers(&me->fdata, CD_MTFACE);
	if (num_layers)
		return std::string(bc_CustomData_get_active_layer_name(&me->fdata, CD_MTFACE));
		
	return "";
}


EffectsExporter::EffectsExporter(COLLADASW::StreamWriter *sw) : COLLADASW::LibraryEffects(sw){}
void EffectsExporter::exportEffects(Scene *sce)
{
	openLibrary();
	MaterialFunctor mf;
	mf.forEachMaterialInScene<EffectsExporter>(sce, *this);

	closeLibrary();
}

void EffectsExporter::operator()(Material *ma, Object *ob)
{
	// create a list of indices to textures of type TEX_IMAGE
	std::vector<int> tex_indices;
	createTextureIndices(ma, tex_indices);

	openEffect(translate_id(id_name(ma)) + "-effect");
	
	COLLADASW::EffectProfile ep(mSW);
	ep.setProfileType(COLLADASW::EffectProfile::COMMON);
	ep.openProfile();
	// set shader type - one of three blinn, phong or lambert
	if (ma->spec_shader == MA_SPEC_BLINN) {
		ep.setShaderType(COLLADASW::EffectProfile::BLINN);
		// shininess
		ep.setShininess(ma->har);
	}
	else if (ma->spec_shader == MA_SPEC_PHONG) {
		ep.setShaderType(COLLADASW::EffectProfile::PHONG);
		// shininess
		ep.setShininess(ma->har);
	}
	else {
		// XXX write warning "Current shader type is not supported" 
		ep.setShaderType(COLLADASW::EffectProfile::LAMBERT);
	}
	// index of refraction
	if (ma->mode & MA_RAYTRANSP) {
		ep.setIndexOfRefraction(ma->ang);
	}
	else {
		ep.setIndexOfRefraction(1.0f);
	}

	COLLADASW::ColorOrTexture cot;

	// transparency
	if (ma->mode & MA_TRANSP) {
		// Tod: because we are in A_ONE mode transparency is calculated like this:
		ep.setTransparency(ma->alpha);
		// cot = getcol(1.0f, 1.0f, 1.0f, 1.0f);
		// ep.setTransparent(cot);
	}

	// emission
	cot=getcol(ma->emit, ma->emit, ma->emit, 1.0f);
	ep.setEmission(cot);

	// diffuse multiplied by diffuse intensity
	cot = getcol(ma->r * ma->ref, ma->g * ma->ref, ma->b * ma->ref, 1.0f);
	ep.setDiffuse(cot);

	// ambient
	cot = getcol(ma->ambr, ma->ambg, ma->ambb, 1.0f);
	ep.setAmbient(cot);

	// reflective, reflectivity
	if (ma->mode & MA_RAYMIRROR) {
		cot = getcol(ma->mirr, ma->mirg, ma->mirb, 1.0f);
		ep.setReflective(cot);
		ep.setReflectivity(ma->ray_mirror);
	}
	// else {
	// 	cot = getcol(ma->specr, ma->specg, ma->specb, 1.0f);
	// 	ep.setReflective(cot);
	// 	ep.setReflectivity(ma->spec);
	// }

	// specular
	if (ep.getShaderType() != COLLADASW::EffectProfile::LAMBERT) {
		cot = getcol(ma->specr * ma->spec, ma->specg * ma->spec, ma->specb * ma->spec, 1.0f);
		ep.setSpecular(cot);
	}	

	// XXX make this more readable if possible

	// create <sampler> and <surface> for each image
	COLLADASW::Sampler samplers[MAX_MTEX];
	//COLLADASW::Surface surfaces[MAX_MTEX];
	//void *samp_surf[MAX_MTEX][2];
	void *samp_surf[MAX_MTEX][1];
	
	// image to index to samp_surf map
	// samp_surf[index] stores 2 pointers, sampler and surface
	std::map<std::string, int> im_samp_map;

	unsigned int a, b;
	for (a = 0, b = 0; a < tex_indices.size(); a++) {
		MTex *t = ma->mtex[tex_indices[a]];
		Image *ima = t->tex->ima;
		
		// Image not set for texture
		if(!ima) continue;
		
		std::string key(id_name(ima));
		key = translate_id(key);

		// create only one <sampler>/<surface> pair for each unique image
		if (im_samp_map.find(key) == im_samp_map.end()) {
			// //<newparam> <surface> <init_from>
			// COLLADASW::Surface surface(COLLADASW::Surface::SURFACE_TYPE_2D,
			// 						   key + COLLADASW::Surface::SURFACE_SID_SUFFIX);
			// COLLADASW::SurfaceInitOption sio(COLLADASW::SurfaceInitOption::INIT_FROM);
			// sio.setImageReference(key);
			// surface.setInitOption(sio);

			// COLLADASW::NewParamSurface surface(mSW);
			// surface->setParamType(COLLADASW::CSW_SURFACE_TYPE_2D);
			
			//<newparam> <sampler> <source>
			COLLADASW::Sampler sampler(COLLADASW::Sampler::SAMPLER_TYPE_2D,
									   key + COLLADASW::Sampler::SAMPLER_SID_SUFFIX,
									   key + COLLADASW::Sampler::SURFACE_SID_SUFFIX);
			sampler.setImageId(key);
			// copy values to arrays since they will live longer
			samplers[a] = sampler;
			//surfaces[a] = surface;
			
			// store pointers so they can be used later when we create <texture>s
			samp_surf[b][0] = &samplers[a];
			//samp_surf[b][1] = &surfaces[a];
			
			im_samp_map[key] = b;
			b++;
		}
	}

	// used as fallback when MTex->uvname is "" (this is pretty common)
	// it is indeed the correct value to use in that case
	std::string active_uv(getActiveUVLayerName(ob));

	// write textures
	// XXX very slow
	for (a = 0; a < tex_indices.size(); a++) {
		MTex *t = ma->mtex[tex_indices[a]];
		Image *ima = t->tex->ima;
		
		// Image not set for texture
		if(!ima) continue;

		// we assume map input is always TEXCO_UV

		std::string key(id_name(ima));
		key = translate_id(key);
		int i = im_samp_map[key];
		COLLADASW::Sampler *sampler = (COLLADASW::Sampler*)samp_surf[i][0];
		//COLLADASW::Surface *surface = (COLLADASW::Surface*)samp_surf[i][1];

		std::string uvname = strlen(t->uvname) ? t->uvname : active_uv;

		// color
		if (t->mapto & MAP_COL) {
			ep.setDiffuse(createTexture(ima, uvname, sampler));
		}
		// ambient
		if (t->mapto & MAP_AMB) {
			ep.setAmbient(createTexture(ima, uvname, sampler));
		}
		// specular
		if (t->mapto & MAP_SPEC) {
			ep.setSpecular(createTexture(ima, uvname, sampler));
		}
		// emission
		if (t->mapto & MAP_EMIT) {
			ep.setEmission(createTexture(ima, uvname, sampler));
		}
		// reflective
		if (t->mapto & MAP_REF) {
			ep.setReflective(createTexture(ima, uvname, sampler));
		}
		// alpha
		if (t->mapto & MAP_ALPHA) {
			ep.setTransparent(createTexture(ima, uvname, sampler));
		}
		// extension:
		// Normal map --> Must be stored with <extra> tag as different technique, 
		// since COLLADA doesn't support normal maps, even in current COLLADA 1.5.
		if (t->mapto & MAP_NORM) {
			COLLADASW::Texture texture(key);
			texture.setTexcoord(uvname);
			texture.setSampler(*sampler);
			// technique FCOLLADA, with the <bump> tag, is most likely the best understood,
			// most widespread de-facto standard.
			texture.setProfileName("FCOLLADA");
			texture.setChildElementName("bump");
			ep.addExtraTechniqueColorOrTexture(COLLADASW::ColorOrTexture(texture));
		}
	}
	// performs the actual writing
	ep.addProfileElements();
	bool twoSided = false;
	if (ob->type == OB_MESH && ob->data) {
		Mesh *me = (Mesh*)ob->data;
		if (me->flag & ME_TWOSIDED)
			twoSided = true;
	}
	if (twoSided)
		ep.addExtraTechniqueParameter("GOOGLEEARTH", "show_double_sided", 1);
	ep.addExtraTechniques(mSW);

	ep.closeProfile();
	if (twoSided)
		mSW->appendTextBlock("<extra><technique profile=\"MAX3D\"><double_sided>1</double_sided></technique></extra>");
	closeEffect();	
}

COLLADASW::ColorOrTexture EffectsExporter::createTexture(Image *ima,
										std::string& uv_layer_name,
										COLLADASW::Sampler *sampler
										/*COLLADASW::Surface *surface*/)
{
	
	COLLADASW::Texture texture(translate_id(id_name(ima)));
	texture.setTexcoord(uv_layer_name);
	//texture.setSurface(*surface);
	texture.setSampler(*sampler);
	
	COLLADASW::ColorOrTexture cot(texture);
	return cot;
}

COLLADASW::ColorOrTexture EffectsExporter::getcol(float r, float g, float b, float a)
{
	COLLADASW::Color color(r,g,b,a);
	COLLADASW::ColorOrTexture cot(color);
	return cot;
}

//returns the array of mtex indices which have image 
//need this for exporting textures
void EffectsExporter::createTextureIndices(Material *ma, std::vector<int> &indices)
{
	indices.clear();

	for (int a = 0; a < MAX_MTEX; a++) {
		if (ma->mtex[a] &&
			ma->mtex[a]->tex &&
			ma->mtex[a]->tex->type == TEX_IMAGE &&
			ma->mtex[a]->texco == TEXCO_UV){
			indices.push_back(a);
		}
	}
}
