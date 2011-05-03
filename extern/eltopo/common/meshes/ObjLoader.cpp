/*    
      glm.c

      Nate Robins, 1997, 2000

      nate@pobox.com, http://www.pobox.com/~nate

 

      Wavefront OBJ model file format reader/writer/manipulator.



      Includes routines for generating smooth normals with

      preservation of edges, welding redundant vertices & texture

      coordinate generation (spheremap and planar projections) + more.

  
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "ObjLoader.hpp"
#include <iostream>
//#include <QGLViewer/vec.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif
#include "vec.h"
//#include <opengl_image.h>

using namespace std;

#define T(x) (model->triangles[(x)])


/* _GLMnode: general purpose node */
typedef struct _GLMnode {
    GLuint         index;
    GLboolean      averaged;
    struct _GLMnode* next;
} GLMnode;


/* glmMax: returns the maximum of two floats */
static GLfloat
glmMax(GLfloat a, GLfloat b) 
{
    if (b > a)
        return b;
    return a;
}

/* glmAbs: returns the absolute value of a float */
static GLfloat
glmAbs(GLfloat f)
{
    if (f < 0)
        return -f;
    return f;
}

/* glmDot: compute the dot product of two vectors
 *
 * u - array of 3 GLfloats (GLfloat u[3])
 * v - array of 3 GLfloats (GLfloat v[3])
 */
static GLfloat
glmDot(GLfloat* u, GLfloat* v)
{
    assert(u); assert(v);
    
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

/* glmCross: compute the cross product of two vectors
 *
 * u - array of 3 GLfloats (GLfloat u[3])
 * v - array of 3 GLfloats (GLfloat v[3])
 * n - array of 3 GLfloats (GLfloat n[3]) to return the cross product in
 */
static GLvoid
glmCross(GLfloat* u, GLfloat* v, GLfloat* n)
{
    assert(u); assert(v); assert(n);
    
    n[0] = u[1]*v[2] - u[2]*v[1];
    n[1] = u[2]*v[0] - u[0]*v[2];
    n[2] = u[0]*v[1] - u[1]*v[0];
}

/* glmNormalize: normalize a vector
 *
 * v - array of 3 GLfloats (GLfloat v[3]) to be normalized
 */
static GLvoid
glmNormalize(GLfloat* v)
{
    GLfloat l;
    
    assert(v);
    
    l = (GLfloat)sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] /= l;
    v[1] /= l;
    v[2] /= l;
}

/* glmEqual: compares two vectors and returns GL_TRUE if they are
 * equal (within a certain threshold) or GL_FALSE if not. An epsilon
 * that works fairly well is 0.000001.
 *
 * u - array of 3 GLfloats (GLfloat u[3])
 * v - array of 3 GLfloats (GLfloat v[3]) 
 */
static GLboolean
glmEqual(GLfloat* u, GLfloat* v, GLfloat epsilon)
{
    if (glmAbs(u[0] - v[0]) < epsilon &&
        glmAbs(u[1] - v[1]) < epsilon &&
        glmAbs(u[2] - v[2]) < epsilon) 
    {
        return GL_TRUE;
    }
    return GL_FALSE;
}

/* glmWeldVectors: eliminate (weld) vectors that are within an
 * epsilon of each other.
 *
 * vectors     - array of GLfloat[3]'s to be welded
 * numvectors - number of GLfloat[3]'s in vectors
 * epsilon     - maximum difference between vectors 
 *
 */
GLfloat*
glmWeldVectors(GLfloat* vectors, GLuint* numvectors, GLfloat epsilon)
{
    GLfloat* copies;
    GLuint   copied;
    GLuint   i, j;
    
    copies = (GLfloat*)malloc(sizeof(GLfloat) * 3 * (*numvectors + 1));
    memcpy(copies, vectors, (sizeof(GLfloat) * 3 * (*numvectors + 1)));
    
    copied = 1;
    for (i = 1; i <= *numvectors; i++) {
        for (j = 1; j <= copied; j++) {
            if (glmEqual(&vectors[3 * i], &copies[3 * j], epsilon)) {
                goto duplicate;
            }
        }
        
        /* must not be any duplicates -- add to the copies array */
        copies[3 * copied + 0] = vectors[3 * i + 0];
        copies[3 * copied + 1] = vectors[3 * i + 1];
        copies[3 * copied + 2] = vectors[3 * i + 2];
        j = copied;             /* pass this along for below */
        copied++;
        
duplicate:
/* set the first component of this vector to point at the correct
        index into the new copies array */
        vectors[3 * i + 0] = (GLfloat)j;
    }
    
    *numvectors = copied-1;
    return copies;
}

/* glmFindGroup: Find a group in the model */
GLMgroup*
glmFindGroup(GLMmodel* model, char* name)
{
    GLMgroup* group;
    
    assert(model);
    
    group = model->groups;
    while(group) {
        if (!strcmp(name, group->name))
            break;
        group = group->next;
    }
    
    return group;
}

/* glmAddGroup: Add a group to the model */
GLMgroup*
glmAddGroup(GLMmodel* model, char* name)
{
    GLMgroup* group;
    
    group = glmFindGroup(model, name);
    if (!group) {
        group = (GLMgroup*)malloc(sizeof(GLMgroup));
        group->name = strdup(name);
        group->material = 0;
        group->numtriangles = 0;
        group->triangles = NULL;
        group->next = model->groups;
        model->groups = group;
        model->numgroups++;
    }
    
    return group;
}

/* glmFindGroup: Find a material in the model */
GLuint
glmFindMaterial(GLMmodel* model, char* name)
{
    GLuint i;
    
    /* XXX doing a linear search on a string key'd list is pretty lame,
    but it works and is fast enough for now. */
    for (i = 0; i < model->nummaterials; i++) {
        if (!strcmp(model->materials[i].name, name))
            goto found;
    }
    
    /* didn't find the name, so print a warning and return the default
    material (0). */
    printf("glmFindMaterial():  can't find material \"%s\".\n", name);
    i = 0;
    
found:
    return i;
}


/* glmDirName: return the directory given a path
 *
 * path - filesystem path
 *
 * NOTE: the return value should be free'd.
 */
static char*
glmDirName(char* path)
{
    char* dir;
    char* s;
    
    dir = strdup(path);
    
    s = strrchr(dir, '/');
    if (s)
        s[1] = '\0';
    else
        dir[0] = '\0';
    
    return dir;
}

/*   written by  M a r c  */
#define MAX_TEXS 64

int glmInitTex(char *name)
{
  static int idTexs[MAX_TEXS], nbTex=0;
  static char nameTexs[MAX_TEXS][256];
  for(int i=0;i<nbTex;i++)
    {
      // search if texture already loaded
      if (!strcmp(nameTexs[i],name)){
	return idTexs[i];
      }
    }
  
  int texWidth, texHeight;

  GLubyte *image=glmReadPPM (name, &texWidth, &texHeight);

  GLuint idOglTex;
  if (image==(GLubyte*)NULL){
    fprintf(stderr,"Sorry, map '%s' not found.\n", name);
    return -1;
  }
  
  fprintf(stderr,"glmInitTex map '%s' loaded (%dx%d).\n", name, texWidth, texHeight);
  glGenTextures(1, &idOglTex );
  glBindTexture(GL_TEXTURE_2D, idOglTex);	idTexs[nbTex]=(int)idOglTex;	nbTex++;
  glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  //glTexEnvi (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  glTexEnvi (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
  gluBuild2DMipmaps (GL_TEXTURE_2D,    // target: 2D texture
		     3,                // 3 components (R, G, and B)
		     texWidth,         // image width
		     texHeight,        // image height
		     GL_RGB,           // image format
		     GL_UNSIGNED_BYTE, // format of data within image file
		     (const void*)image);		// image file
  fprintf(stderr,"glmInitTex finito. tex %d is %s\n", (int)idOglTex, name);
  return (int)idOglTex;
}

void glmUseTex(int id){
	// fprintf(stderr,"glmUseTex %d\n",id);
	if (id>=0){
		glBindTexture(GL_TEXTURE_2D, (GLuint)id);
		glEnable(GL_TEXTURE_2D);
	}
}

/* glmReadMTL: read a wavefront material library file
 *
 * model - properly initialized GLMmodel structure
 * name  - name of the material library
 */
static GLvoid
glmReadMTL(GLMmodel* model, char* name, GLuint texId)
{
  FILE* file;
  char* dir;
  char* filename;
  char    buf[128];
  GLuint nummaterials, i;
    
  dir = glmDirName(model->pathname);
  filename = (char*)malloc(sizeof(char) * (strlen(dir) + strlen(name) + 1));
  strcpy(filename, dir);
  strcat(filename, name);
  free(dir);
    
  file = fopen(filename, "r");
  if (!file) {
    fprintf(stderr, "glmReadMTL() failed: can't open material file \"%s\".\n",
            filename);
    exit(1);
  }
  free(filename);
    
  /* count the number of materials in the file */
  nummaterials = 1;
  while(fscanf(file, "%s", buf) != EOF) {
    switch(buf[0]) {
    case '#':               /* comment */
      /* eat up rest of line */
      fgets(buf, sizeof(buf), file);
      break;
    case 'n':               /* newmtl */
      fgets(buf, sizeof(buf), file);
      nummaterials++;
      sscanf(buf, "%s %s", buf, buf);
      break;
    default:
      /* eat up rest of line */
      fgets(buf, sizeof(buf), file);
      break;
    }
  }
    
  rewind(file);
    
  model->materials = (GLMmaterial*)malloc(sizeof(GLMmaterial) * nummaterials);
  model->nummaterials = nummaterials;
    
  /* set the default material */
  for (i = 0; i < nummaterials; i++) {
    model->materials[i].name = NULL;
    model->materials[i].shininess = 65.0f;
    model->materials[i].diffuse[0] = 0.8f;
    model->materials[i].diffuse[1] = 0.8f;
    model->materials[i].diffuse[2] = 0.8f;
    model->materials[i].diffuse[3] = 1.0f;
    model->materials[i].ambient[0] = 0.2f;
    model->materials[i].ambient[1] = 0.2f;
    model->materials[i].ambient[2] = 0.2f;
    model->materials[i].ambient[3] = 1.0f;
    model->materials[i].specular[0] = 0.0f;
    model->materials[i].specular[1] = 0.0f;
    model->materials[i].specular[2] = 0.0f;
    model->materials[i].specular[3] = 1.0f;
  }
  model->materials[0].name = strdup("default");
    
  /* now, read in the data */
  nummaterials = 0;
  while(fscanf(file, "%s", buf) != EOF) {
    switch(buf[0]) {
    case '#':               /* comment */
      /* eat up rest of line */
      fgets(buf, sizeof(buf), file);
      break;
    case 'n':               /* newmtl */
      fgets(buf, sizeof(buf), file);
      sscanf(buf, "%s %s", buf, buf);
      nummaterials++;
      model->materials[nummaterials].name = strdup(buf);
      break;
    case 'N':
      fscanf(file, "%f", &model->materials[nummaterials].shininess);
      /* wavefront shininess is from [0, 1000], so scale for OpenGL */
      model->materials[nummaterials].shininess /= 1000.0f;
      model->materials[nummaterials].shininess *= 128.0f;
      break;
      /*  M a r c  */
    case 'M':
    case 'm':
      fgets(buf, sizeof(buf), file);
      sscanf(buf, "%s %s", buf, buf);
      //model->materials[nummaterials].idTex=glmInitTex(buf);
      model->materials[nummaterials].idTex=texId;
      break;
    case 'K':
      switch(buf[1]) {
      case 'd':
	fscanf(file, "%f %f %f",
	       &model->materials[nummaterials].diffuse[0],
	       &model->materials[nummaterials].diffuse[1],
	       &model->materials[nummaterials].diffuse[2]);
	break;
      case 's':
	fscanf(file, "%f %f %f",
	       &model->materials[nummaterials].specular[0],
	       &model->materials[nummaterials].specular[1],
	       &model->materials[nummaterials].specular[2]);
	break;
      case 'a':
	fscanf(file, "%f %f %f",
	       &model->materials[nummaterials].ambient[0],
	       &model->materials[nummaterials].ambient[1],
	       &model->materials[nummaterials].ambient[2]);
	break;
      default:
	/* eat up rest of line */
	fgets(buf, sizeof(buf), file);
	break;
      }
      break;
    default:
      /* eat up rest of line */
      fgets(buf, sizeof(buf), file);
      break;
    }
  }
					 
  fclose(file);
}

/* glmWriteMTL: write a wavefront material library file
 *
 * model   - properly initialized GLMmodel structure
 * modelpath  - pathname of the model being written
 * mtllibname - name of the material library to be written
 */
static GLvoid
glmWriteMTL(GLMmodel* model, char* modelpath, char* mtllibname)
{
    FILE* file;
    char* dir;
    char* filename;
    GLMmaterial* material;
    GLuint i;
    
    dir = glmDirName(modelpath);
    filename = (char*)malloc(sizeof(char) * (strlen(dir)+strlen(mtllibname)));
    strcpy(filename, dir);
    strcat(filename, mtllibname);
    free(dir);
    
    /* open the file */
    file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "glmWriteMTL() failed: can't open file \"%s\".\n",
            filename);
        exit(1);
    }
    free(filename);
    
    /* spit out a header */
    fprintf(file, "#  \n");
    fprintf(file, "#  Wavefront MTL generated by GLM library\n");
    fprintf(file, "#  \n");
    fprintf(file, "#  GLM library\n");
    fprintf(file, "#  Nate Robins\n");
    fprintf(file, "#  ndr@pobox.com\n");
    fprintf(file, "#  http://www.pobox.com/~ndr\n");
    fprintf(file, "#  \n\n");
    
    for (i = 0; i < model->nummaterials; i++) {
        material = &model->materials[i];
        fprintf(file, "newmtl %s\n", material->name);
        fprintf(file, "Ka %f %f %f\n", 
            material->ambient[0], material->ambient[1], material->ambient[2]);
        fprintf(file, "Kd %f %f %f\n", 
            material->diffuse[0], material->diffuse[1], material->diffuse[2]);
        fprintf(file, "Ks %f %f %f\n", 
            material->specular[0],material->specular[1],material->specular[2]);
        fprintf(file, "Ns %f\n", material->shininess / 128.0 * 1000.0);
        fprintf(file, "\n");
    }
}


/* glmFirstPass: first pass at a Wavefront OBJ file that gets all the
 * statistics of the model (such as #vertices, #normals, etc)
 *
 * model - properly initialized GLMmodel structure
 * file  - (fopen'd) file descriptor 
 */
static GLvoid
glmFirstPass(GLMmodel* model, FILE* file, GLuint texID) 
{
    GLuint  numvertices;        /* number of vertices in model */
    GLuint  numnormals;         /* number of normals in model */
    GLuint  numtexcoords;       /* number of texcoords in model */
    GLuint  numtriangles;       /* number of triangles in model */
    GLMgroup* group;            /* current group */
    unsigned    v, n, t;
    char        buf[128];
    
    /* make a default group */
    group = glmAddGroup(model, "default");
    
    numvertices = numnormals = numtexcoords = numtriangles = 0;
    while(fscanf(file, "%s", buf) != EOF) {
        switch(buf[0]) {
        case '#':               /* comment */
            /* eat up rest of line */
            fgets(buf, sizeof(buf), file);
            break;
        case 'v':               /* v, vn, vt */
            switch(buf[1]) {
            case '\0':          /* vertex */
                /* eat up rest of line */
                fgets(buf, sizeof(buf), file);
                numvertices++;
                break;
            case 'n':           /* normal */
                /* eat up rest of line */
                fgets(buf, sizeof(buf), file);
                numnormals++;
                break;
            case 't':           /* texcoord */
                /* eat up rest of line */
                fgets(buf, sizeof(buf), file);
                numtexcoords++;
                break;
            default:
                printf("glmFirstPass(): Unknown token \"%s\".\n", buf);
                exit(1);
                break;
            }
            break;
            case 'm':
                fgets(buf, sizeof(buf), file);
                sscanf(buf, "%s %s", buf, buf);
                model->mtllibname = strdup(buf);
                glmReadMTL(model, buf, texID);
                break;
            case 'u':
                /* eat up rest of line */
                fgets(buf, sizeof(buf), file);
                break;
            case 'g':               /* group */
                /* eat up rest of line */
                fgets(buf, sizeof(buf), file);
#if SINGLE_STRING_GROUP_NAMES
                sscanf(buf, "%s", buf);
#else
                buf[strlen(buf)-1] = '\0';  /* nuke '\n' */
#endif
                group = glmAddGroup(model, buf);
                break;
            case 'f':               /* face */
                v = n = t = 0;
                fscanf(file, "%s", buf);
                /* can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d */
                if (strstr(buf, "//")) {
                    /* v//n */
                    sscanf(buf, "%u//%u", &v, &n);
                    fscanf(file, "%u//%u", &v, &n);
                    fscanf(file, "%u//%u", &v, &n);
                    numtriangles++;
                    group->numtriangles++;
                    while(fscanf(file, "%u//%u", &v, &n) > 0) {
                        numtriangles++;
                        group->numtriangles++;
                    }
                } else if (sscanf(buf, "%u/%u/%u", &v, &t, &n) == 3) {
                    /* v/t/n */
                    fscanf(file, "%u/%u/%u", &v, &t, &n);
                    fscanf(file, "%u/%u/%u", &v, &t, &n);
                    numtriangles++;
                    group->numtriangles++;
                    while(fscanf(file, "%u/%u/%u", &v, &t, &n) > 0) {
                        numtriangles++;
                        group->numtriangles++;
                    }
                } else if (sscanf(buf, "%u/%u", &v, &t) == 2) {
                    /* v/t */
                    fscanf(file, "%u/%u", &v, &t);
                    fscanf(file, "%u/%u", &v, &t);
                    numtriangles++;
                    group->numtriangles++;
                    while(fscanf(file, "%u/%u", &v, &t) > 0) {
                        numtriangles++;
                        group->numtriangles++;
                    }
                } else {
                    /* v */
                    fscanf(file, "%u", &v);
                    fscanf(file, "%u", &v);
                    numtriangles++;
                    group->numtriangles++;
                    while(fscanf(file, "%u", &v) > 0) {
                        numtriangles++;
                        group->numtriangles++;
                    }
                }
                break;
                
            default:
                /* eat up rest of line */
                fgets(buf, sizeof(buf), file);
                break;
        }
  }
  
  /* set the stats in the model structure */
  model->numvertices  = numvertices;
  model->numnormals   = numnormals;
  model->numtexcoords = numtexcoords;
  model->numtriangles = numtriangles;
  
  /* allocate memory for the triangles in each group */
  group = model->groups;
  while(group) {
      group->triangles = (GLuint*)malloc(sizeof(GLuint) * group->numtriangles);
      group->numtriangles = 0;
      group = group->next;
  }
}

/* glmSecondPass: second pass at a Wavefront OBJ file that gets all
 * the data.
 *
 * model - properly initialized GLMmodel structure
 * file  - (fopen'd) file descriptor 
 */
static GLvoid
glmSecondPass(GLMmodel* model, FILE* file) 
{
    GLuint  numvertices;        /* number of vertices in model */
    GLuint  numnormals;         /* number of normals in model */
    GLuint  numtexcoords;       /* number of texcoords in model */
    GLuint  numtriangles;       /* number of triangles in model */
    GLfloat*    vertices;           /* array of vertices  */
    GLfloat*    normals;            /* array of normals */
    GLfloat*    texcoords;          /* array of texture coordinates */
    GLMgroup* group;            /* current group pointer */
    GLuint  material;           /* current material */
    GLuint  v, n, t;
    char        buf[128];
    
    /* set the pointer shortcuts */
    vertices       = model->vertices;
    normals    = model->normals;
    texcoords    = model->texcoords;
    group      = model->groups;
    
    /* on the second pass through the file, read all the data into the
    allocated arrays */
    numvertices = numnormals = numtexcoords = 1;
    numtriangles = 0;
    material = 0;
    while(fscanf(file, "%s", buf) != EOF) {
        switch(buf[0]) {
        case '#':               /* comment */
            /* eat up rest of line */
            fgets(buf, sizeof(buf), file);
            break;
        case 'v':               /* v, vn, vt */
            switch(buf[1]) {
            case '\0':          /* vertex */
                fscanf(file, "%f %f %f", 
                    &vertices[3 * numvertices + 0], 
                    &vertices[3 * numvertices + 1], 
                    &vertices[3 * numvertices + 2]);
                numvertices++;
                break;
            case 'n':           /* normal */
                fscanf(file, "%f %f %f", 
                    &normals[3 * numnormals + 0],
                    &normals[3 * numnormals + 1], 
                    &normals[3 * numnormals + 2]);
                numnormals++;
                break;
            case 't':           /* texcoord */
                fscanf(file, "%f %f", 
                    &texcoords[2 * numtexcoords + 0],
                    &texcoords[2 * numtexcoords + 1]);
                numtexcoords++;
                break;
            }
            break;
            case 'u':
                fgets(buf, sizeof(buf), file);
                sscanf(buf, "%s %s", buf, buf);
                group->material = material = glmFindMaterial(model, buf);
                break;
            case 'g':               /* group */
                /* eat up rest of line */
                fgets(buf, sizeof(buf), file);
#if SINGLE_STRING_GROUP_NAMES
                sscanf(buf, "%s", buf);
#else
                buf[strlen(buf)-1] = '\0';  /* nuke '\n' */
#endif
                group = glmFindGroup(model, buf);
                group->material = material;
                break;
            case 'f':               /* face */
                v = n = t = 0;
                fscanf(file, "%s", buf);
                /* can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d */
                if (strstr(buf, "//")) {
                    /* v//n */
                    sscanf(buf, "%d//%d", (int*)&v, (int*)&n);
                    T(numtriangles).vindices[0] = v;
                    T(numtriangles).nindices[0] = n;
                    fscanf(file, "%d//%d", (int*)&v, (int*)&n);
                    T(numtriangles).vindices[1] = v;
                    T(numtriangles).nindices[1] = n;
                    fscanf(file, "%d//%d", (int*)&v, (int*)&n);
                    T(numtriangles).vindices[2] = v;
                    T(numtriangles).nindices[2] = n;
                    group->triangles[group->numtriangles++] = numtriangles;
                    numtriangles++;
                    while(fscanf(file, "%d//%d", (int*)&v, (int*)&n) > 0) {
                        T(numtriangles).vindices[0] = T(numtriangles-1).vindices[0];
                        T(numtriangles).nindices[0] = T(numtriangles-1).nindices[0];
                        T(numtriangles).vindices[1] = T(numtriangles-1).vindices[2];
                        T(numtriangles).nindices[1] = T(numtriangles-1).nindices[2];
                        T(numtriangles).vindices[2] = v;
                        T(numtriangles).nindices[2] = n;
                        group->triangles[group->numtriangles++] = numtriangles;
                        numtriangles++;
                    }
                } else if (sscanf(buf, "%d/%d/%d", (int*)&v, (int*)&t, (int*)&n) == 3) {
                    /* v/t/n */
                    T(numtriangles).vindices[0] = v;
                    T(numtriangles).tindices[0] = t;
                    T(numtriangles).nindices[0] = n;
                    fscanf(file, "%d/%d/%d", (int*)&v, (int*)&t, (int*)&n);
                    T(numtriangles).vindices[1] = v;
                    T(numtriangles).tindices[1] = t;
                    T(numtriangles).nindices[1] = n;
                    fscanf(file, "%d/%d/%d", (int*)&v, (int*)&t, (int*)&n);
                    T(numtriangles).vindices[2] = v;
                    T(numtriangles).tindices[2] = t;
                    T(numtriangles).nindices[2] = n;
                    group->triangles[group->numtriangles++] = numtriangles;
                    numtriangles++;
                    while(fscanf(file, "%d/%d/%d", (int*)&v, (int*)&t, (int*)&n) > 0) {
                        T(numtriangles).vindices[0] = T(numtriangles-1).vindices[0];
                        T(numtriangles).tindices[0] = T(numtriangles-1).tindices[0];
                        T(numtriangles).nindices[0] = T(numtriangles-1).nindices[0];
                        T(numtriangles).vindices[1] = T(numtriangles-1).vindices[2];
                        T(numtriangles).tindices[1] = T(numtriangles-1).tindices[2];
                        T(numtriangles).nindices[1] = T(numtriangles-1).nindices[2];
                        T(numtriangles).vindices[2] = v;
                        T(numtriangles).tindices[2] = t;
                        T(numtriangles).nindices[2] = n;
                        group->triangles[group->numtriangles++] = numtriangles;
                        numtriangles++;
                    }
                } else if (sscanf(buf, "%d/%d", (int*)&v, (int*)&t) == 2) {
                    /* v/t */
                    T(numtriangles).vindices[0] = v;
                    T(numtriangles).tindices[0] = t;
                    fscanf(file, "%d/%d", (int*)&v, (int*)&t);
                    T(numtriangles).vindices[1] = v;
                    T(numtriangles).tindices[1] = t;
                    fscanf(file, "%d/%d", (int*)&v, (int*)&t);
                    T(numtriangles).vindices[2] = v;
                    T(numtriangles).tindices[2] = t;
                    group->triangles[group->numtriangles++] = numtriangles;
                    numtriangles++;
                    while(fscanf(file, "%d/%d", (int*)&v, (int*)&t) > 0) {
                        T(numtriangles).vindices[0] = T(numtriangles-1).vindices[0];
                        T(numtriangles).tindices[0] = T(numtriangles-1).tindices[0];
                        T(numtriangles).vindices[1] = T(numtriangles-1).vindices[2];
                        T(numtriangles).tindices[1] = T(numtriangles-1).tindices[2];
                        T(numtriangles).vindices[2] = v;
                        T(numtriangles).tindices[2] = t;
                        group->triangles[group->numtriangles++] = numtriangles;
                        numtriangles++;
                    }
                } else {
                    /* v */
                    sscanf(buf, "%d", (int*)&v);
                    T(numtriangles).vindices[0] = v;
                    fscanf(file, "%d", (int*)&v);
                    T(numtriangles).vindices[1] = v;
                    fscanf(file, "%d", (int*)&v);
                    T(numtriangles).vindices[2] = v;
                    group->triangles[group->numtriangles++] = numtriangles;
                    numtriangles++;
                    while(fscanf(file, "%d", (int*)&v) > 0) {
                        T(numtriangles).vindices[0] = T(numtriangles-1).vindices[0];
                        T(numtriangles).vindices[1] = T(numtriangles-1).vindices[2];
                        T(numtriangles).vindices[2] = v;
                        group->triangles[group->numtriangles++] = numtriangles;
                        numtriangles++;
                    }
                }
                break;
                
            default:
                /* eat up rest of line */
                fgets(buf, sizeof(buf), file);
                break;
    }
  }
  
#if 1
  /* announce the memory requirements */
  printf(" Memory: %lu bytes\n",
      numvertices  * 3*sizeof(GLfloat) +
      numnormals   * 3*sizeof(GLfloat) * (numnormals ? 1 : 0) +
      numtexcoords * 3*sizeof(GLfloat) * (numtexcoords ? 1 : 0) +
      numtriangles * sizeof(GLMtriangle));
#endif
}


/* public functions */


/* glmUnitize: "unitize" a model by translating it to the origin and
 * scaling it to fit in a unit cube around the origin.   Returns the
 * scalefactor used.
 *
 * model - properly initialized GLMmodel structure 
 */
GLfloat
glmUnitize(GLMmodel* model)
{
    GLuint  i;
    GLfloat maxx, minx, maxy, miny, maxz, minz;
    GLfloat cx, cy, cz, w, h, d;
    GLfloat scale;
    
    assert(model);
    assert(model->vertices);
    
    /* get the max/mins */
    maxx = minx = model->vertices[3 + 0];
    maxy = miny = model->vertices[3 + 1];
    maxz = minz = model->vertices[3 + 2];
    for (i = 1; i <= model->numvertices; i++) {
        if (maxx < model->vertices[3 * i + 0])
            maxx = model->vertices[3 * i + 0];
        if (minx > model->vertices[3 * i + 0])
            minx = model->vertices[3 * i + 0];
        
        if (maxy < model->vertices[3 * i + 1])
            maxy = model->vertices[3 * i + 1];
        if (miny > model->vertices[3 * i + 1])
            miny = model->vertices[3 * i + 1];
        
        if (maxz < model->vertices[3 * i + 2])
            maxz = model->vertices[3 * i + 2];
        if (minz > model->vertices[3 * i + 2])
            minz = model->vertices[3 * i + 2];
    }
    
    /* calculate model width, height, and depth */
    w = glmAbs(maxx) + glmAbs(minx);
    h = glmAbs(maxy) + glmAbs(miny);
    d = glmAbs(maxz) + glmAbs(minz);
    
    /* calculate center of the model */
    cx = (maxx + minx) / 2.0f;
    cy = (maxy + miny) / 2.0f;
    cz = (maxz + minz) / 2.0f;
    
    /* calculate unitizing scale factor */
    scale = 2.0f / glmMax(glmMax(w, h), d);
    
    /* translate around center then scale */
    for (i = 1; i <= model->numvertices; i++) {
        model->vertices[3 * i + 0] -= cx;
        model->vertices[3 * i + 1] -= cy;
        model->vertices[3 * i + 2] -= cz;
        model->vertices[3 * i + 0] *= scale;
        model->vertices[3 * i + 1] *= scale;
        model->vertices[3 * i + 2] *= scale;
    }
    
    return scale;
}

/* glmDimensions: Calculates the dimensions (width, height, depth) of
 * a model.
 *
 * model   - initialized GLMmodel structure
 * dimensions - array of 3 GLfloats (GLfloat dimensions[3])
 */
GLvoid
glmDimensions(GLMmodel* model, GLfloat* dimensions)
{
    GLuint i;
    GLfloat maxx, minx, maxy, miny, maxz, minz;
    
    assert(model);
    assert(model->vertices);
    assert(dimensions);
    
    /* get the max/mins */
    maxx = minx = model->vertices[3 + 0];
    maxy = miny = model->vertices[3 + 1];
    maxz = minz = model->vertices[3 + 2];
    for (i = 1; i <= model->numvertices; i++) {
        if (maxx < model->vertices[3 * i + 0])
            maxx = model->vertices[3 * i + 0];
        if (minx > model->vertices[3 * i + 0])
            minx = model->vertices[3 * i + 0];
        
        if (maxy < model->vertices[3 * i + 1])
            maxy = model->vertices[3 * i + 1];
        if (miny > model->vertices[3 * i + 1])
            miny = model->vertices[3 * i + 1];
        
        if (maxz < model->vertices[3 * i + 2])
            maxz = model->vertices[3 * i + 2];
        if (minz > model->vertices[3 * i + 2])
            minz = model->vertices[3 * i + 2];
    }
    
    /* calculate model width, height, and depth */
    dimensions[0] = glmAbs(maxx) + glmAbs(minx);
    dimensions[1] = glmAbs(maxy) + glmAbs(miny);
    dimensions[2] = glmAbs(maxz) + glmAbs(minz);
}

/* glmScale: Scales a model by a given amount.
 * 
 * model - properly initialized GLMmodel structure
 * scale - scalefactor (0.5 = half as large, 2.0 = twice as large)
 */
GLvoid
glmScale(GLMmodel* model, GLfloat scale)
{
    GLuint i;
    
    for (i = 1; i <= model->numvertices; i++) {
        model->vertices[3 * i + 0] *= scale;
        model->vertices[3 * i + 1] *= scale;
        model->vertices[3 * i + 2] *= scale;
    }
}

/* glmTranslate: Translates along direction d a model by a given amount.
 * 
 * model - properly initialized GLMmodel structure
 * transl - translationfactor 
 */
GLvoid
glmTranslate(GLMmodel* model, GLfloat transl, int d)
{
    GLuint i;
    
    for (i = 1; i <= model->numvertices; i++)
	{
		switch (d)
		{
			case 0 :
				// translation selon x
        		model->vertices[3 * i + 0] += transl;
        		break;
			case 1 :
				// translation selon y
				model->vertices[3 * i + 1] += transl;
        		break;
			case 2 :
				// translation selon z
				model->vertices[3 * i + 2] += transl;
				break;
		}
    }
}

/* glmReverseWinding: Reverse the polygon winding for all polygons in
 * this model.   Default winding is counter-clockwise.  Also changes
 * the direction of the normals.
 * 
 * model - properly initialized GLMmodel structure 
 */
GLvoid
glmReverseWinding(GLMmodel* model)
{
    GLuint i, swap;
    
    assert(model);
    
    for (i = 0; i < model->numtriangles; i++) {
        swap = T(i).vindices[0];
        T(i).vindices[0] = T(i).vindices[2];
        T(i).vindices[2] = swap;
        
        if (model->numnormals) {
            swap = T(i).nindices[0];
            T(i).nindices[0] = T(i).nindices[2];
            T(i).nindices[2] = swap;
        }
        
        if (model->numtexcoords) {
            swap = T(i).tindices[0];
            T(i).tindices[0] = T(i).tindices[2];
            T(i).tindices[2] = swap;
        }
    }
    
    /* reverse facet normals */
    for (i = 1; i <= model->numfacetnorms; i++) {
        model->facetnorms[3 * i + 0] = -model->facetnorms[3 * i + 0];
        model->facetnorms[3 * i + 1] = -model->facetnorms[3 * i + 1];
        model->facetnorms[3 * i + 2] = -model->facetnorms[3 * i + 2];
    }
    
    /* reverse vertex normals */
    for (i = 1; i <= model->numnormals; i++) {
        model->normals[3 * i + 0] = -model->normals[3 * i + 0];
        model->normals[3 * i + 1] = -model->normals[3 * i + 1];
        model->normals[3 * i + 2] = -model->normals[3 * i + 2];
    }
}

/* glmFacetNormals: Generates facet normals for a model (by taking the
 * cross product of the two vectors derived from the sides of each
 * triangle).  Assumes a counter-clockwise winding.
 *
 * model - initialized GLMmodel structure
 */
GLvoid
glmFacetNormals(GLMmodel* model)
{
    GLuint  i;
    GLfloat u[3];
    GLfloat v[3];
    
    assert(model);
    assert(model->vertices);
    
    /* clobber any old facetnormals */
    if (model->facetnorms)
        free(model->facetnorms);
    
    /* allocate memory for the new facet normals */
    model->numfacetnorms = model->numtriangles;
    model->facetnorms = (GLfloat*)malloc(sizeof(GLfloat) *
                       3 * (model->numfacetnorms + 1));
    
    for (i = 0; i < model->numtriangles; i++) {
        model->triangles[i].findex = i+1;
        
        u[0] = model->vertices[3 * T(i).vindices[1] + 0] -
            model->vertices[3 * T(i).vindices[0] + 0];
        u[1] = model->vertices[3 * T(i).vindices[1] + 1] -
            model->vertices[3 * T(i).vindices[0] + 1];
        u[2] = model->vertices[3 * T(i).vindices[1] + 2] -
            model->vertices[3 * T(i).vindices[0] + 2];
        
        v[0] = model->vertices[3 * T(i).vindices[2] + 0] -
            model->vertices[3 * T(i).vindices[0] + 0];
        v[1] = model->vertices[3 * T(i).vindices[2] + 1] -
            model->vertices[3 * T(i).vindices[0] + 1];
        v[2] = model->vertices[3 * T(i).vindices[2] + 2] -
            model->vertices[3 * T(i).vindices[0] + 2];
        
        glmCross(u, v, &model->facetnorms[3 * (i+1)]);
        glmNormalize(&model->facetnorms[3 * (i+1)]);
    }
}

/* glmVertexNormals: Generates smooth vertex normals for a model.
 * First builds a list of all the triangles each vertex is in.   Then
 * loops through each vertex in the the list averaging all the facet
 * normals of the triangles each vertex is in.   Finally, sets the
 * normal index in the triangle for the vertex to the generated smooth
 * normal.   If the dot product of a facet normal and the facet normal
 * associated with the first triangle in the list of triangles the
 * current vertex is in is greater than the cosine of the angle
 * parameter to the function, that facet normal is not added into the
 * average normal calculation and the corresponding vertex is given
 * the facet normal.  This tends to preserve hard edges.  The angle to
 * use depends on the model, but 90 degrees is usually a good start.
 *
 * model - initialized GLMmodel structure
 * angle - maximum angle (in degrees) to smooth across
 */
GLvoid
glmVertexNormals(GLMmodel* model, GLfloat angle)
{
    GLMnode*    node;
    GLMnode*    tail;
    GLMnode** members;
    GLfloat*    normals;
    GLuint  numnormals;
    GLfloat average[3];
    GLfloat dot, cos_angle;
    GLuint  i, avg;
    
    assert(model);
    assert(model->facetnorms);
    
    /* calculate the cosine of the angle (in degrees) */
    cos_angle = cos(angle * (float)M_PI / 180.0f);
    
    /* nuke any previous normals */
    if (model->normals)
        free(model->normals);
    
    /* allocate space for new normals */
    model->numnormals = model->numtriangles * 3; /* 3 normals per triangle */
    model->normals = (GLfloat*)malloc(sizeof(GLfloat)* 3* (model->numnormals+1));
    
    /* allocate a structure that will hold a linked list of triangle
    indices for each vertex */
    members = (GLMnode**)malloc(sizeof(GLMnode*) * (model->numvertices + 1));
    for (i = 1; i <= model->numvertices; i++)
        members[i] = NULL;
    
    /* for every triangle, create a node for each vertex in it */
    for (i = 0; i < model->numtriangles; i++) {
        node = (GLMnode*)malloc(sizeof(GLMnode));
        node->index = i;
        node->next  = members[T(i).vindices[0]];
        members[T(i).vindices[0]] = node;
        
        node = (GLMnode*)malloc(sizeof(GLMnode));
        node->index = i;
        node->next  = members[T(i).vindices[1]];
        members[T(i).vindices[1]] = node;
        
        node = (GLMnode*)malloc(sizeof(GLMnode));
        node->index = i;
        node->next  = members[T(i).vindices[2]];
        members[T(i).vindices[2]] = node;
    }
    
    /* calculate the average normal for each vertex */
    numnormals = 1;
    for (i = 1; i <= model->numvertices; i++) {
    /* calculate an average normal for this vertex by averaging the
        facet normal of every triangle this vertex is in */
        node = members[i];
        if (!node)
            fprintf(stderr, "glmVertexNormals(): vertex w/o a triangle\n");
        average[0] = 0.0f; average[1] = 0.0f; average[2] = 0.0f;
        avg = 0;
        while (node) {
        /* only average if the dot product of the angle between the two
        facet normals is greater than the cosine of the threshold
        angle -- or, said another way, the angle between the two
            facet normals is less than (or equal to) the threshold angle */
            dot = glmDot(&model->facetnorms[3 * T(node->index).findex],
                &model->facetnorms[3 * T(members[i]->index).findex]);
            if (dot > cos_angle) {
                node->averaged = GL_TRUE;
                average[0] += model->facetnorms[3 * T(node->index).findex + 0];
                average[1] += model->facetnorms[3 * T(node->index).findex + 1];
                average[2] += model->facetnorms[3 * T(node->index).findex + 2];
                avg = 1;            /* we averaged at least one normal! */
            } else {
                node->averaged = GL_FALSE;
            }
            node = node->next;
        }
        
        if (avg) {
            /* normalize the averaged normal */
            glmNormalize(average);
            
            /* add the normal to the vertex normals list */
            model->normals[3 * numnormals + 0] = average[0];
            model->normals[3 * numnormals + 1] = average[1];
            model->normals[3 * numnormals + 2] = average[2];
            avg = numnormals;
            numnormals++;
        }
        
        /* set the normal of this vertex in each triangle it is in */
        node = members[i];
        while (node) {
            if (node->averaged) {
                /* if this node was averaged, use the average normal */
                if (T(node->index).vindices[0] == i)
                    T(node->index).nindices[0] = avg;
                else if (T(node->index).vindices[1] == i)
                    T(node->index).nindices[1] = avg;
                else if (T(node->index).vindices[2] == i)
                    T(node->index).nindices[2] = avg;
            } else {
                /* if this node wasn't averaged, use the facet normal */
                model->normals[3 * numnormals + 0] = 
                    model->facetnorms[3 * T(node->index).findex + 0];
                model->normals[3 * numnormals + 1] = 
                    model->facetnorms[3 * T(node->index).findex + 1];
                model->normals[3 * numnormals + 2] = 
                    model->facetnorms[3 * T(node->index).findex + 2];
                if (T(node->index).vindices[0] == i)
                    T(node->index).nindices[0] = numnormals;
                else if (T(node->index).vindices[1] == i)
                    T(node->index).nindices[1] = numnormals;
                else if (T(node->index).vindices[2] == i)
                    T(node->index).nindices[2] = numnormals;
                numnormals++;
            }
            node = node->next;
        }
    }
    
    model->numnormals = numnormals - 1;
    
    /* free the member information */
    for (i = 1; i <= model->numvertices; i++) {
        node = members[i];
        while (node) {
            tail = node;
            node = node->next;
            free(tail);
        }
    }
    free(members);
    
    /* pack the normals array (we previously allocated the maximum
    number of normals that could possibly be created (numtriangles *
    3), so get rid of some of them (usually alot unless none of the
    facet normals were averaged)) */
    normals = model->normals;
    model->normals = (GLfloat*)malloc(sizeof(GLfloat)* 3* (model->numnormals+1));
    for (i = 1; i <= model->numnormals; i++) {
        model->normals[3 * i + 0] = normals[3 * i + 0];
        model->normals[3 * i + 1] = normals[3 * i + 1];
        model->normals[3 * i + 2] = normals[3 * i + 2];
    }
    free(normals);
}


/* glmLinearTexture: Generates texture coordinates according to a
 * linear projection of the texture map.  It generates these by
 * linearly mapping the vertices onto a square.
 *
 * model - pointer to initialized GLMmodel structure
 */
GLvoid
glmLinearTexture(GLMmodel* model)
{
    GLMgroup *group;
    GLfloat dimensions[3];
    GLfloat x, y, scalefactor;
    GLuint i;
    
    assert(model);
    
    if (model->texcoords)
        free(model->texcoords);
    model->numtexcoords = model->numvertices;
    model->texcoords=(GLfloat*)malloc(sizeof(GLfloat)*2*(model->numtexcoords+1));
    
    glmDimensions(model, dimensions);
    scalefactor = 2.0f / 
        glmAbs(glmMax(glmMax(dimensions[0], dimensions[1]), dimensions[2]));
    
    /* do the calculations */
    for(i = 1; i <= model->numvertices; i++) {
        x = model->vertices[3 * i + 0] * scalefactor;
        y = model->vertices[3 * i + 2] * scalefactor;
        model->texcoords[2 * i + 0] = (x + 1.0f) / 2.0f;
        model->texcoords[2 * i + 1] = (y + 1.0f) / 2.0f;
    }
    
    /* go through and put texture coordinate indices in all the triangles */
    group = model->groups;
    while(group) {
        for(i = 0; i < group->numtriangles; i++) {
            T(group->triangles[i]).tindices[0] = T(group->triangles[i]).vindices[0];
            T(group->triangles[i]).tindices[1] = T(group->triangles[i]).vindices[1];
            T(group->triangles[i]).tindices[2] = T(group->triangles[i]).vindices[2];
        }    
        group = group->next;
    }
    
#if 0
    printf("glmLinearTexture(): generated %d linear texture coordinates\n",
        model->numtexcoords);
#endif
}

/* glmSpheremapTexture: Generates texture coordinates according to a
 * spherical projection of the texture map.  Sometimes referred to as
 * spheremap, or reflection map texture coordinates.  It generates
 * these by using the normal to calculate where that vertex would map
 * onto a sphere.  Since it is impossible to map something flat
 * perfectly onto something spherical, there is distortion at the
 * poles.  This particular implementation causes the poles along the X
 * axis to be distorted.
 *
 * model - pointer to initialized GLMmodel structure
 */
GLvoid
glmSpheremapTexture(GLMmodel* model)
{
    GLMgroup* group;
    GLfloat theta, phi, rho, x, y, z, r;
    GLuint i;
    
    assert(model);
    assert(model->normals);
    
    if (model->texcoords)
        free(model->texcoords);
    model->numtexcoords = model->numnormals;
    model->texcoords=(GLfloat*)malloc(sizeof(GLfloat)*2*(model->numtexcoords+1));
    
    for (i = 1; i <= model->numnormals; i++) {
        z = model->normals[3 * i + 0];  /* re-arrange for pole distortion */
        y = model->normals[3 * i + 1];
        x = model->normals[3 * i + 2];
        r = sqrt((x * x) + (y * y));
        rho = sqrt((r * r) + (z * z));
        
        if(r == 0.0) {
            theta = 0.0f;
            phi = 0.0f;
        } else {
            if(z == 0.0)
                phi = 3.14159265f / 2.0f;
            else
                phi = acos(z / rho);
            
            if(y == 0.0)
                theta = 3.141592365f / 2.0f;
            else
                theta = asin(y / r) + (3.14159265f / 2.0f);
        }
        
        model->texcoords[2 * i + 0] = theta / 3.14159265f;
        model->texcoords[2 * i + 1] = phi / 3.14159265f;
    }
    
    /* go through and put texcoord indices in all the triangles */
    group = model->groups;
    while(group) {
        for (i = 0; i < group->numtriangles; i++) {
            T(group->triangles[i]).tindices[0] = T(group->triangles[i]).nindices[0];
            T(group->triangles[i]).tindices[1] = T(group->triangles[i]).nindices[1];
            T(group->triangles[i]).tindices[2] = T(group->triangles[i]).nindices[2];
        }
        group = group->next;
    }
}

/* glmDelete: Deletes a GLMmodel structure.
 *
 * model - initialized GLMmodel structure
 */
GLvoid
glmDelete(GLMmodel* model)
{
  GLMgroup* group;
  GLuint i;
    
  assert(model);
    
  if (model->pathname)     free(model->pathname);
  if (model->mtllibname) free(model->mtllibname);
  if (model->vertices)     free(model->vertices);
  if (model->normals)  free(model->normals);
  if (model->texcoords)  free(model->texcoords);
  if (model->facetnorms) free(model->facetnorms);
  if (model->triangles)  free(model->triangles);
  if (model->materials)
    {
      for (i = 0; i < model->nummaterials; i++)
	free(model->materials[i].name);
    }
  free(model->materials);
  while(model->groups)
    {
      group = model->groups;
      model->groups = model->groups->next;
      if (group->name)
	{
	  // cout<<"on efface "<<group->name<<endl;
	  free(group->name);
	}
      if (group->triangles)
	free(group->triangles);
      free(group);
    }
  free(model);
}

/* glmReadOBJ: Reads a model description from a Wavefront .OBJ file.
 * Returns a pointer to the created object which should be free'd with
 * glmDelete().
 *
 * filename - name of the file containing the Wavefront .OBJ format data.  
 */
GLMmodel* 
glmReadOBJ(char* filename, GLuint texId)
{
    GLMmodel* model;
    FILE*   file;
    
    /* open the file */
    file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "glmReadOBJ() failed: can't open data file \"%s\".\n",
            filename);
        exit(1);
    }
    
    /* allocate a new model */
    model = (GLMmodel*)malloc(sizeof(GLMmodel));
    model->pathname    = strdup(filename);
    model->mtllibname    = NULL;
    model->numvertices   = 0;
    model->vertices    = NULL;
    model->numnormals    = 0;
    model->normals     = NULL;
    model->numtexcoords  = 0;
    model->texcoords       = NULL;
    model->numfacetnorms = 0;
    model->facetnorms    = NULL;
    model->numtriangles  = 0;
    model->triangles       = NULL;
    model->nummaterials  = 0;
    model->materials       = NULL;
    model->numgroups       = 0;
    model->groups      = NULL;
    model->position[0]   = 0.0f;
    model->position[1]   = 0.0f;
    model->position[2]   = 0.0f;
    
    /* make a first pass through the file to get a count of the number
    of vertices, normals, texcoords & triangles */
    glmFirstPass(model, file, texId);
    
    /* allocate memory */
    model->vertices = (GLfloat*)malloc(sizeof(GLfloat) *
        3 * (model->numvertices + 1));
    model->triangles = (GLMtriangle*)malloc(sizeof(GLMtriangle) *
        model->numtriangles);
    if (model->numnormals) {
        model->normals = (GLfloat*)malloc(sizeof(GLfloat) *
            3 * (model->numnormals + 1));
    }
    if (model->numtexcoords) {
        model->texcoords = (GLfloat*)malloc(sizeof(GLfloat) *
            2 * (model->numtexcoords + 1));
    }
    
    /* rewind to beginning of file and read in the data this pass */
    rewind(file);
    
    glmSecondPass(model, file);
    
    /* close the file */
    fclose(file);
    
    return model;
}

/* glmWriteOBJ: Writes a model description in Wavefront .OBJ format to
 * a file.
 *
 * model - initialized GLMmodel structure
 * filename - name of the file to write the Wavefront .OBJ format data to
 * mode  - a bitwise or of values describing what is written to the file
 *             GLM_NONE     -  render with only vertices
 *             GLM_FLAT     -  render with facet normals
 *             GLM_SMOOTH   -  render with vertex normals
 *             GLM_TEXTURE  -  render with texture coords
 *             GLM_COLOR    -  render with colors (color material)
 *             GLM_MATERIAL -  render with materials
 *             GLM_COLOR and GLM_MATERIAL should not both be specified.  
 *             GLM_FLAT and GLM_SMOOTH should not both be specified.  
 */
GLvoid
glmWriteOBJ(GLMmodel* model, char* filename, GLuint mode)
{
    GLuint  i;
    FILE*   file;
    GLMgroup* group;
    
    assert(model);
    
    /* do a bit of warning */
    if (mode & GLM_FLAT && !model->facetnorms) {
        printf("glmWriteOBJ() warning: flat normal output requested "
            "with no facet normals defined.\n");
        mode &= ~GLM_FLAT;
    }
    if (mode & GLM_SMOOTH && !model->normals) {
        printf("glmWriteOBJ() warning: smooth normal output requested "
            "with no normals defined.\n");
        mode &= ~GLM_SMOOTH;
    }
    if (mode & GLM_TEXTURE && !model->texcoords) {
        printf("glmWriteOBJ() warning: texture coordinate output requested "
            "with no texture coordinates defined.\n");
        mode &= ~GLM_TEXTURE;
    }
    if (mode & GLM_FLAT && mode & GLM_SMOOTH) {
        printf("glmWriteOBJ() warning: flat normal output requested "
            "and smooth normal output requested (using smooth).\n");
        mode &= ~GLM_FLAT;
    }
    if (mode & GLM_COLOR && !model->materials) {
        printf("glmWriteOBJ() warning: color output requested "
            "with no colors (materials) defined.\n");
        mode &= ~GLM_COLOR;
    }
    if (mode & GLM_MATERIAL && !model->materials) {
        printf("glmWriteOBJ() warning: material output requested "
            "with no materials defined.\n");
        mode &= ~GLM_MATERIAL;
    }
    if (mode & GLM_COLOR && mode & GLM_MATERIAL) {
        printf("glmWriteOBJ() warning: color and material output requested "
            "outputting only materials.\n");
        mode &= ~GLM_COLOR;
    }
    
    
    /* open the file */
    file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "glmWriteOBJ() failed: can't open file \"%s\" to write.\n",
            filename);
        exit(1);
    }
    
    /* spit out a header */
    fprintf(file, "#  \n");
    fprintf(file, "#  Wavefront OBJ generated by GLM library\n");
    fprintf(file, "#  \n");
    fprintf(file, "#  GLM library\n");
    fprintf(file, "#  Nate Robins\n");
    fprintf(file, "#  ndr@pobox.com\n");
    fprintf(file, "#  http://www.pobox.com/~ndr\n");
    fprintf(file, "#  \n");
    
    if (mode & GLM_MATERIAL && model->mtllibname) {
        fprintf(file, "\nmtllib %s\n\n", model->mtllibname);
        glmWriteMTL(model, filename, model->mtllibname);
    }
    
    /* spit out the vertices */
    fprintf(file, "\n");
    fprintf(file, "# %d vertices\n", model->numvertices);
    for (i = 1; i <= model->numvertices; i++) {
        fprintf(file, "v %f %f %f\n", 
            model->vertices[3 * i + 0],
            model->vertices[3 * i + 1],
            model->vertices[3 * i + 2]);
    }
    
    /* spit out the smooth/flat normals */
    if (mode & GLM_SMOOTH) {
        fprintf(file, "\n");
        fprintf(file, "# %d normals\n", model->numnormals);
        for (i = 1; i <= model->numnormals; i++) {
            fprintf(file, "vn %f %f %f\n", 
                model->normals[3 * i + 0],
                model->normals[3 * i + 1],
                model->normals[3 * i + 2]);
        }
    } else if (mode & GLM_FLAT) {
        fprintf(file, "\n");
        fprintf(file, "# %d normals\n", model->numfacetnorms);
        for (i = 1; i <= model->numnormals; i++) {
            fprintf(file, "vn %f %f %f\n", 
                model->facetnorms[3 * i + 0],
                model->facetnorms[3 * i + 1],
                model->facetnorms[3 * i + 2]);
        }
    }
    
    /* spit out the texture coordinates */
    if (mode & GLM_TEXTURE) {
        fprintf(file, "\n");
        fprintf(file, "# %d texcoords\n", model->numtexcoords);
        for (i = 1; i <= model->numtexcoords; i++) {
            fprintf(file, "vt %f %f\n", 
                model->texcoords[2 * i + 0],
                model->texcoords[2 * i + 1]);
        }
    }
    
    fprintf(file, "\n");
    fprintf(file, "# %d groups\n", model->numgroups);
    fprintf(file, "# %d faces (triangles)\n", model->numtriangles);
    fprintf(file, "\n");
    
    group = model->groups;
    while(group) {
        fprintf(file, "g %s\n", group->name);
        if (mode & GLM_MATERIAL)
            fprintf(file, "usemtl %s\n", model->materials[group->material].name);
        for (i = 0; i < group->numtriangles; i++) {
            if (mode & GLM_SMOOTH && mode & GLM_TEXTURE) {
                fprintf(file, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
                    T(group->triangles[i]).vindices[0], 
                    T(group->triangles[i]).nindices[0], 
                    T(group->triangles[i]).tindices[0],
                    T(group->triangles[i]).vindices[1],
                    T(group->triangles[i]).nindices[1],
                    T(group->triangles[i]).tindices[1],
                    T(group->triangles[i]).vindices[2],
                    T(group->triangles[i]).nindices[2],
                    T(group->triangles[i]).tindices[2]);
            } else if (mode & GLM_FLAT && mode & GLM_TEXTURE) {
                fprintf(file, "f %d/%d %d/%d %d/%d\n",
                    T(group->triangles[i]).vindices[0],
                    T(group->triangles[i]).findex,
                    T(group->triangles[i]).vindices[1],
                    T(group->triangles[i]).findex,
                    T(group->triangles[i]).vindices[2],
                    T(group->triangles[i]).findex);
            } else if (mode & GLM_TEXTURE) {
                fprintf(file, "f %d/%d %d/%d %d/%d\n",
                    T(group->triangles[i]).vindices[0],
                    T(group->triangles[i]).tindices[0],
                    T(group->triangles[i]).vindices[1],
                    T(group->triangles[i]).tindices[1],
                    T(group->triangles[i]).vindices[2],
                    T(group->triangles[i]).tindices[2]);
            } else if (mode & GLM_SMOOTH) {
                fprintf(file, "f %d//%d %d//%d %d//%d\n",
                    T(group->triangles[i]).vindices[0],
                    T(group->triangles[i]).nindices[0],
                    T(group->triangles[i]).vindices[1],
                    T(group->triangles[i]).nindices[1],
                    T(group->triangles[i]).vindices[2], 
                    T(group->triangles[i]).nindices[2]);
            } else if (mode & GLM_FLAT) {
                fprintf(file, "f %d//%d %d//%d %d//%d\n",
                    T(group->triangles[i]).vindices[0], 
                    T(group->triangles[i]).findex,
                    T(group->triangles[i]).vindices[1],
                    T(group->triangles[i]).findex,
                    T(group->triangles[i]).vindices[2],
                    T(group->triangles[i]).findex);
            } else {
                fprintf(file, "f %d %d %d\n",
                    T(group->triangles[i]).vindices[0],
                    T(group->triangles[i]).vindices[1],
                    T(group->triangles[i]).vindices[2]);
            }
        }
        fprintf(file, "\n");
        group = group->next;
    }
    
    fclose(file);
}

/* glmDraw: Renders the model to the current OpenGL context using the
 * mode specified.
 *
 * model - initialized GLMmodel structure
 * mode  - a bitwise OR of values describing what is to be rendered.
 *             GLM_NONE     -  render with only vertices
 *             GLM_FLAT     -  render with facet normals
 *             GLM_SMOOTH   -  render with vertex normals
 *             GLM_TEXTURE  -  render with texture coords
 *             GLM_COLOR    -  render with colors (color material)
 *             GLM_MATERIAL -  render with materials
 *             GLM_COLOR and GLM_MATERIAL should not both be specified.  
 *             GLM_FLAT and GLM_SMOOTH should not both be specified.  
 */
GLvoid
//glmDraw(GLMmodel* model, GLuint mode, const bool names, const int selectTriangle)
glmDraw(GLMmodel* model, GLuint mode, const float basisColor[3], const float alpha, const bool names, std::map<int, TriangleIndex> map)
{
  static GLuint i;
  static GLMgroup* group;
  static GLMtriangle* triangle;
  static GLMmaterial* material;
    
  assert(model);
  assert(model->vertices);
    
  /* do a bit of warning */
  if (mode & GLM_FLAT && !model->facetnorms) {
    printf("glmDraw() warning: flat render mode requested "
	   "with no facet normals defined.\n");
    mode &= ~GLM_FLAT;
  }
  if (mode & GLM_SMOOTH && !model->normals) {
    printf("glmDraw() warning: smooth render mode requested "
	   "with no normals defined.\n");
    mode &= ~GLM_SMOOTH;
  }
  if (mode & GLM_TEXTURE && !model->texcoords) {
    printf("glmDraw() warning: texture render mode requested "
	   "with no texture coordinates defined.\n");
    mode &= ~GLM_TEXTURE;
  }
  if (mode & GLM_FLAT && mode & GLM_SMOOTH) {
    printf("glmDraw() warning: flat render mode requested "
	   "and smooth render mode requested (using smooth).\n");
    mode &= ~GLM_FLAT;
  }
  if (mode & GLM_COLOR && !model->materials) {
    printf("glmDraw() warning: color render mode requested "
	   "with no materials defined.\n");
    mode &= ~GLM_COLOR;
  }
  if (mode & GLM_MATERIAL && !model->materials) {
    printf("glmDraw() warning: material render mode requested "
	   "with no materials defined.\n");
    mode &= ~GLM_MATERIAL;
  }
  if (mode & GLM_COLOR && mode & GLM_MATERIAL) {
    printf("glmDraw() warning: color and material render mode requested "
	   "using only material mode.\n");
    mode &= ~GLM_COLOR;
  }
  if (mode & GLM_COLOR)
    glEnable(GL_COLOR_MATERIAL);
  else if (mode & GLM_MATERIAL)
    glDisable(GL_COLOR_MATERIAL);
    
  /* perhaps this loop should be unrolled into material, color, flat,
     smooth, etc. loops?  since most cpu's have good branch prediction
     schemes (and these branches will always go one way), probably
     wouldn't gain too much?  */
  
  group = model->groups;

  // FLO
  unsigned short nTrianglesPrev = 0;
  unsigned short int groupID = 0;
  // --

  if (!names)
    {
      while (group)
	{
	  // FLO
	  //cout<<"groupID = "<<groupID<<endl;
	  // --
      
	  if (mode & GLM_MATERIAL)
	    {
	      material = &model->materials[group->material];
	      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material->ambient);
	      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material->diffuse);
	      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material->specular);
	      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material->shininess);
	    }

	  /*  M a r c  */
	  if (mode & GLM_TEXTURE)
	    {
	      glmUseTex(material->idTex);
	    }
	  else glDisable(GL_TEXTURE_2D);
    
	  
	  if (mode & GLM_COLOR)
	    {
	      glColor3fv(material->diffuse);
	    }
        
	  glBegin(GL_TRIANGLES);
	  for (i = 0; i < group->numtriangles; i++)
	    {
	      // FLO
	  
	      // Global index of the triangle i of current group
	      //   ie. index of the triangle in the whole model
	      GLuint triangleID = (group->triangles)[i];
	  
	      // selected triangles are drawn in a special color
	      // IMPROVE the test !!
	      GLuint triangleIndex = (map[triangleID]).getIndex();
	      //  if (map[(int)triangleID] == triangleID)
	      if (triangleIndex == triangleID)
		{
		  //cout<<"triangleID = "<<triangleID<<endl;
		  //glColor3f(0.30, 0.20, 0.1);
		  glColor3f(0.15f, 0.1f, 0.05f);
		}
	      else
		// Access to an empty triangle
		{
		  assert(triangleIndex == EMPTY_TRIANGLE);
		  //glColor3fv(basisColor);
		  glColor4f(basisColor[0], basisColor[1], basisColor[2], alpha);
		}
	      // --

	      triangle = &T(group->triangles[i]);

	      // FLO

	      // NB : in 'vertices' array of the model, points are stores as blocks of 3 floats
	      //   [P0x, P0y, P0z, P1x, P1y, P1z, P2x, P2y, P2z, ...]
	  
	      // So to access to point Pk knowing index k, you have to compute vertices[3*k] (which
	      // gives the index of Pkx), and then get the block of 3 floats starting at this address

	      // Here, 3*triangle->vindices[k] is the starting index of the point Pk in 'vertices' array
	    
	      // --
	  
	      if (mode & GLM_FLAT)
		glNormal3fv(&model->facetnorms[3 * triangle->findex]);
            
	      if (mode & GLM_SMOOTH)
		glNormal3fv(&model->normals[3 * triangle->nindices[0]]);
	      if (mode & GLM_TEXTURE)
		glTexCoord2fv(&model->texcoords[2 * triangle->tindices[0]]);
	      
	      glVertex3fv(&model->vertices[3 * triangle->vindices[0]]);
            
	      if (mode & GLM_SMOOTH)
		glNormal3fv(&model->normals[3 * triangle->nindices[1]]);
	      if (mode & GLM_TEXTURE)
		glTexCoord2fv(&model->texcoords[2 * triangle->tindices[1]]);
	      
	      glVertex3fv(&model->vertices[3 * triangle->vindices[1]]);
            
	      if (mode & GLM_SMOOTH)
		glNormal3fv(&model->normals[3 * triangle->nindices[2]]);
	      if (mode & GLM_TEXTURE)
		glTexCoord2fv(&model->texcoords[2 * triangle->tindices[2]]);
	      
	      glVertex3fv(&model->vertices[3 * triangle->vindices[2]]);

	    }
      
	  glEnd();

	  // FLO
	  // Updating nTrianglesPrev
	  nTrianglesPrev = group->numtriangles;
	  //cout<<"numTriangles = "<<nTrianglesPrev<<endl;
	  // --
      
	  group = group->next;

	  // FLO
	  ++groupID;
	  // --
	}
    }

  // FLO
  else
    // names == true
    {
      group = model->groups;
      //group = group->next;
      //group = group->next;

      while (group)
	{
	  for (i = 0; i < group->numtriangles; i++)
	    {		  
	      // Global index of the triangle i of current group
	      //   ie. index of the triangle in the whole model
	      int triangleID = (group->triangles)[i];
	  	  
	      glPushName(triangleID);

	      glBegin(GL_TRIANGLES);

	      // triangle = &((model->triangle)[group->triangles[i]])
	      //   gives
	      triangle = &T(group->triangles[i]);

	      // NB : in 'vertices' array of the model, points are stores as blocks of 3 floats
	      //   [P0x, P0y, P0z, P1x, P1y, P1z, P2x, P2y, P2z, ...]
	  
	      // So to access to point Pk knowing index k, you have to compute vertices[3*k] (which
	      // gives the index of Pkx), and then get the block of 3 floats starting at this address

	      // Here, 3*triangle->vindices[k] is the starting index of the point Pk in 'vertices' array

	      glVertex3fv(&model->vertices[3*triangle->vindices[0]]);
	      glVertex3fv(&model->vertices[3*triangle->vindices[1]]);
	      glVertex3fv(&model->vertices[3*triangle->vindices[2]]);

	      glEnd();
	  
	      glPopName();
	    } // end for
	  group = group->next;
	} // end while
    } // end else
  // --
}


GLvoid
glmDrawSelectedTriangles(GLMmodel* model, GLuint mode, std::map<int, TriangleIndex> * selectedTriangles)
{

  //cout<<"Entre dans glmDrawSelectedTriangles"<<endl;

  //cout<<"map size = "<<selectedTriangles->size()<<endl;
  static GLMtriangle* triangle;
    
  assert(model);
  assert(model->vertices);
    
  /* do a bit of warning */
  if (mode & GLM_FLAT && !model->facetnorms) {
    printf("glmDraw() warning: flat render mode requested "
	   "with no facet normals defined.\n");
    mode &= ~GLM_FLAT;
  }
  if (mode & GLM_SMOOTH && !model->normals) {
    printf("glmDraw() warning: smooth render mode requested "
	   "with no normals defined.\n");
    mode &= ~GLM_SMOOTH;
  }
  if (mode & GLM_TEXTURE && !model->texcoords) {
    printf("glmDraw() warning: texture render mode requested "
	   "with no texture coordinates defined.\n");
    mode &= ~GLM_TEXTURE;
  }
  if (mode & GLM_FLAT && mode & GLM_SMOOTH) {
    printf("glmDraw() warning: flat render mode requested "
	   "and smooth render mode requested (using smooth).\n");
    mode &= ~GLM_FLAT;
  }
  if (mode & GLM_COLOR && !model->materials) {
    printf("glmDraw() warning: color render mode requested "
	   "with no materials defined.\n");
    mode &= ~GLM_COLOR;
  }
  if (mode & GLM_MATERIAL && !model->materials) {
    printf("glmDraw() warning: material render mode requested "
	   "with no materials defined.\n");
    mode &= ~GLM_MATERIAL;
  }
  if (mode & GLM_COLOR && mode & GLM_MATERIAL) {
    printf("glmDraw() warning: color and material render mode requested "
	   "using only material mode.\n");
    mode &= ~GLM_COLOR;
  }
  if (mode & GLM_COLOR)
    glEnable(GL_COLOR_MATERIAL);
  else if (mode & GLM_MATERIAL)
    glDisable(GL_COLOR_MATERIAL);
    
  /* perhaps this loop should be unrolled into material, color, flat,
     smooth, etc. loops?  since most cpu's have good branch prediction
     schemes (and these branches will always go one way), probably
     wouldn't gain too much?  */


  glColor3f(0.0f, 0.7f, 1.0f);

  std::map<int, TriangleIndex>::const_iterator mapIter;

  glBegin(GL_TRIANGLES);
  for (mapIter = selectedTriangles->begin(); mapIter != selectedTriangles->end(); ++mapIter)
    {      

      glColor3f(0.0f, 0.7f, 1.0f);
      // Global index of the triangle i of current group
      //   ie. index of the triangle in the whole model
      GLuint triangleID = (mapIter->second).getIndex();

      cout<<"triangleID = "<<triangleID<<endl;
      triangle = &T(triangleID);
	  
      if (mode & GLM_FLAT)
	glNormal3fv(&model->facetnorms[3 * triangle->findex]);
            
      if (mode & GLM_SMOOTH)
	glNormal3fv(&model->normals[3 * triangle->nindices[0]]);
      if (mode & GLM_TEXTURE)
	glTexCoord2fv(&model->texcoords[2 * triangle->tindices[0]]);

      glVertex3fv(&model->vertices[3 * triangle->vindices[0]]);
            
      if (mode & GLM_SMOOTH)
	glNormal3fv(&model->normals[3 * triangle->nindices[1]]);
      if (mode & GLM_TEXTURE)
	glTexCoord2fv(&model->texcoords[2 * triangle->tindices[1]]);

      glVertex3fv(&model->vertices[3 * triangle->vindices[1]]);
            
      if (mode & GLM_SMOOTH)
	glNormal3fv(&model->normals[3 * triangle->nindices[2]]);
      if (mode & GLM_TEXTURE)
	glTexCoord2fv(&model->texcoords[2 * triangle->tindices[2]]);

      glVertex3fv(&model->vertices[3 * triangle->vindices[2]]);	
            
    } // end for
  glEnd();
}


GLMmodel * extractSubModel(GLMmodel * model, const std::map<int, TriangleIndex>& selectedTrianglesMap, char * name)
{
  GLMmodel * extrModel;
  
  // allocate a new model
  extrModel = (GLMmodel*)malloc(sizeof(GLMmodel));
  extrModel->pathname    = NULL;
  extrModel->mtllibname    = NULL;
  extrModel->numvertices   = 0;
  extrModel->vertices    = NULL;
  extrModel->numnormals    = 0;
  extrModel->normals     = NULL;
  extrModel->numtexcoords  = 0;
  extrModel->texcoords       = NULL;
  extrModel->numfacetnorms = 0;
  extrModel->facetnorms    = NULL;
  extrModel->numtriangles  = 0;
  extrModel->triangles       = NULL;
  extrModel->nummaterials  = 0;
  extrModel->materials       = NULL;
  extrModel->numgroups       = 0;
  extrModel->groups      = NULL;
  extrModel->position[0]   = model->position[0];
  extrModel->position[1]   = model->position[1];
  extrModel->position[2]   = model->position[2];

  //cout<<"Position modèle = "<<extrModel->position[0]<<", "<<extrModel->position[1]<<", "<<extrModel->position[2]<<endl;
  

  
  // Gettint the number of vertices, normals, texcoords and triangles

  // VERTICES

  // Building a vertices map to get the number of vertices + their values
  //    map : index of the vertex in the WHOLE model --> (value of the vertex, index of the vertex in the EXTRACTED model)

  std::map<GLuint, std::pair<Vec3d, GLuint> > selectedVerticesMap;


  // Browsing the triangles map to find the number of vertices (getting rid of doublons)
  
  std::map<int, TriangleIndex>::const_iterator trianglesMapIter;

  for (trianglesMapIter = selectedTrianglesMap.begin(); trianglesMapIter != selectedTrianglesMap.end(); ++trianglesMapIter)
    {
      // index of the corresponding triangle in the whole model
      GLuint index = (trianglesMapIter->second).getIndex();

      // finding the corresponding triangle in the whole model
      GLMtriangle triangle = (model->triangles)[index];

      // indexes of the 3 vertices of the triangle in the whole model
      GLuint v1index, v2index, v3index;

      v1index = triangle.vindices[0];
      v2index = triangle.vindices[1];
      v3index = triangle.vindices[2];

      // values of the corresponding vertices
      Vec3d v1, v2, v3;

      for (int k = 0; k != 3; ++k)
	{
	  v1[k] = model->vertices[3*v1index + k];
	  v2[k] = model->vertices[3*v2index + k];
	  v3[k] = model->vertices[3*v3index + k];
	}

      // filling vertices map
      //   for the moment, with a fake new index (to be corrected during the filling step)
      selectedVerticesMap[v1index] = std::pair<Vec3d, GLuint>(v1, 0);
      selectedVerticesMap[v2index] = std::pair<Vec3d, GLuint>(v2, 0);
      selectedVerticesMap[v3index] = std::pair<Vec3d, GLuint>(v3, 0);
    }
  
  extrModel->numvertices = (unsigned int)selectedVerticesMap.size();

  //cout<<"numvertices = "<<extrModel->numvertices<<endl;
  
  // allocate memory
  // !! CAREFUL : indexation begins at 1 !
  extrModel->vertices = (GLfloat*)malloc(sizeof(GLfloat) *
					 3 * (extrModel->numvertices + 1));

  // filling
  std::map<GLuint, std::pair<Vec3d, GLuint> >::iterator verticesMapIter;
  unsigned short int j = 1;
  for (verticesMapIter = selectedVerticesMap.begin(); verticesMapIter != selectedVerticesMap.end(); ++verticesMapIter)
    {
      // index of the vertex in the whole model
      //GLuint index = verticesMapIter->first;
      // corresponding vertex
      Vec3d vertex = (verticesMapIter->second).first;
      
      // filling the vertices array of extrModel
      for (int k = 0; k != 3; ++k)
         extrModel->vertices[3*j + k] = (float)vertex[k];

      // correcting the vertices map with the proper new index
      (verticesMapIter->second).second = j;
      
      ++j;
    }

  
  // TRIANGLES
  extrModel->numtriangles = (unsigned int)selectedTrianglesMap.size();

  //cout<<"numtriangles = "<<extrModel->numtriangles<<endl;
  
  // allocation
  extrModel->triangles = (GLMtriangle*)malloc(sizeof(GLMtriangle) *
					      extrModel->numtriangles);
  // filling
  //std::map<int, GLuint>::const_iterator trianglesMapIter;
  //unsigned short j = 0;
  j = 0;
  for (trianglesMapIter = selectedTrianglesMap.begin(); trianglesMapIter != selectedTrianglesMap.end(); ++trianglesMapIter)
    {
      // index of the corresponding triangle in the whole model
      int index = trianglesMapIter->first;

      //cout<<"index = "<<index<<endl;
      
      // building the corresponding triangle
      GLMtriangle triangle;

      for (int k = 0; k != 3; ++k)
	{
	  GLuint oldVindk = (model->triangles[index]).vindices[k]; 
	  triangle.vindices[k] = (selectedVerticesMap[oldVindk]).second;
	  triangle.nindices[k] = (model->triangles[index]).nindices[k];
	  triangle.tindices[k] = (model->triangles[index]).tindices[k];
	}
      triangle.findex = (model->triangles[index]).findex;

      // filling the triangles array of extrModel
      extrModel->triangles[j] = triangle;
      
      ++j;
    }

  


  // GROUPS
  extrModel->numgroups = 1;

  // allocate memory
  extrModel->groups = (GLMgroup*)malloc(sizeof(GLMgroup) *
					extrModel->numgroups);

  (extrModel->groups)->name = (char*)malloc(sizeof(char) * strlen(name) + 1);

  (extrModel->groups)->triangles = (GLuint*)malloc(sizeof(GLuint) *
						   extrModel->numtriangles);

  // filling

  //(extrModel->groups)->name = NULL;//"scalp";
  strcpy((extrModel->groups)->name, name);
  
  (extrModel->groups)->numtriangles = extrModel->numtriangles;

  for (GLuint i = 0; i != extrModel->numtriangles; ++i)
    {
      ((extrModel->groups)->triangles)[i] = i;
    }
  
  (extrModel->groups)->material = 0;
  (extrModel->groups)->next = NULL;


  
  //glmFacetNormals(extrModel);
  
  /*
    if (extrModel->numnormals)
    {
    extrModel->normals = (GLfloat*)malloc(sizeof(GLfloat) *
    3 * (extrModel->numnormals + 1));
    }
  
    if (extrModel->numtexcoords)
    {
    extrModel->texcoords = (GLfloat*)malloc(sizeof(GLfloat) *
    2 * (extrModel->numtexcoords + 1));
    }
  */
  return extrModel;
}


//void glmBuildNeighborsMap(GLMmodel * model, std::map<GLuint, TriangleIDVector>& neighborsMap)
void glmBuildNeighborsMap(GLMmodel * model, std::map<GLuint, TriangleIDTriplet>& neighborsMap)
{
  std::map<Edge, TriangleIDVector> adjacencyMap;
  GLMtriangle triangle;
  int v0, v1, v2;

  // Clearing the neigborsMap
  neighborsMap.clear();
  
  // Building an adjacent map based on the edges
  
  for (GLuint i = 0; i != model->numtriangles; ++i)
    {
      // Current triangle
      triangle = model->triangles[i];

      // The 3 vertices of the triangle
      v0 = triangle.vindices[0];
      v1 = triangle.vindices[1];
      v2 = triangle.vindices[2];

      // Filling the adjacency map
      Edge e0(v0, v1);
      Edge e1(v1, v2);
      Edge e2(v2, v0);
      
      adjacencyMap[e0].push_back(TriangleIndex(i));
      adjacencyMap[e1].push_back(TriangleIndex(i));
      adjacencyMap[e2].push_back(TriangleIndex(i));
    }
  
  /*
  // Showing the adjacency map
  cout<<endl<<"Adjacency map"<<endl;
  cout<<"size = "<<adjacencyMap.size()<<endl;
  for (std::map<Edge, TriangleIDVector>::const_iterator it = adjacencyMap.begin(); it != adjacencyMap.end(); ++it)
  {
  cout<<"Edge "<<it->first;
  for (TriangleIDVector::const_iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2)
  cout<<it2->getIndex()<<" ";
  cout<<endl;
  }
  */
  
  // Building the neighbors map from the adjacent map
  for (GLuint i = 0; i != model->numtriangles; ++i)
    {
      //cout<<endl<<"i = "<<i<<" *************"<<endl;
      // Current triangle
      triangle = model->triangles[i];

      // The 3 vertices of the triangle
      v0 = triangle.vindices[0];
      v1 = triangle.vindices[1];
      v2 = triangle.vindices[2];

      // The vector that contains all the triangles adjacent to the current one
      //   (with no doublon)
      //TriangleIDVector concat;
      TriangleIDTriplet concat;
      
      // Filling this concac vector using the adjacency map

      TriangleIDVector::const_iterator iter;

      // First edge *****************************
      TriangleIDVector adj = adjacencyMap[Edge(v0, v1)];
      assert((adj.size() > 0) && (adj.size() <= 2));

      if (adj.size() == 1)
	{
	  assert((adj.begin())->getIndex() == i);
	  //concat.push_back(TriangleIndex(EMPTY_TRIANGLE));
	  concat.setFirstIndex(TriangleIndex(EMPTY_TRIANGLE));
	}
      else
	{
	  for (iter = adj.begin(); iter != adj.end(); ++iter)
	    {
	      if (iter->getIndex() != i)
		//concat.push_back(*iter);
		concat.setFirstIndex(*iter);
	    }
	}
      
      
      // Second edge *****************************
      adj = adjacencyMap[Edge(v1, v2)];
      assert((adj.size() > 0) && (adj.size() <= 2));

     
      if (adj.size() == 1)
	{
	  assert((adj.begin())->getIndex() == i);
	  //concat.push_back(TriangleIndex(EMPTY_TRIANGLE));
	  concat.setSecondIndex(TriangleIndex(EMPTY_TRIANGLE));
	}
      else
	{
	  for (iter = adj.begin(); iter != adj.end(); ++iter)
	    {
	      if (iter->getIndex() != i)
		//concat.push_back(*iter);
		concat.setSecondIndex(*iter);
	    }
	}
      
      
      // Third edge *****************************
      adj = adjacencyMap[Edge(v2, v0)];
      assert((adj.size() > 0) && (adj.size() <= 2));

      if (adj.size() == 1)
	{
	  assert((adj.begin())->getIndex() == i);
	  //concat.push_back(TriangleIndex(EMPTY_TRIANGLE));
	  concat.setThirdIndex(TriangleIndex(EMPTY_TRIANGLE));
	}
      else
	{
	  for (iter = adj.begin(); iter != adj.end(); ++iter)
	    {
	      if (iter->getIndex() != i)
		//concat.push_back(*iter);
		concat.setThirdIndex(*iter);
	    }
	}
      
      //cout<<"concat.size = "<<concat.size()<<endl;
      //assert(concat.size() == 3);

      // Filling the neighbors map
      /*
      for (iter = concat.begin(); iter != concat.end(); ++iter)
	{
	  neighborsMap[i].push_back(*iter);
	}
      assert(neighborsMap[i].size() == 3);
      */
      neighborsMap[i] = concat;
    }

  
  // Showing the neighbors map
  cout<<endl<<"Neighbors map"<<endl;
  cout<<"size = "<<(unsigned int)neighborsMap.size()<<endl;
  for (std::map<GLuint, TriangleIDTriplet>::const_iterator nit = neighborsMap.begin(); nit != neighborsMap.end(); ++nit)
    {
      cout<<"Triangle "<<nit->first<<" : ";
      cout<<nit->second;
      cout<<endl;
    }
  
}

/*
  void glmBuildEdgePointsMap(GLMmodel * model, std::map<int, qglviewer::Vec>& edgePointsMap)
  {
  std::map<Edge, TriangleIDVector> adjacencyMap;
  GLMtriangle triangle;
  int v0, v1, v2;

  // Clearing the edgesPointMap
  edgePointsMap.clear();
  
  // Building an adjacent map based on the edges
  
  for (GLuint i = 0; i != model->numtriangles; ++i)
  {
  // Current triangle
  triangle = model->triangles[i];

  // The 3 vertices of the triangle
  v0 = triangle.vindices[0];
  v1 = triangle.vindices[1];
  v2 = triangle.vindices[2];

  // Filling the adjacency map
  Edge e0(v0, v1);
  Edge e1(v1, v2);
  Edge e2(v2, v0);
      
  adjacencyMap[e0].push_back(TriangleIndex(i));
  adjacencyMap[e1].push_back(TriangleIndex(i));
  adjacencyMap[e2].push_back(TriangleIndex(i));
  }
  
  //
  //// Showing the adjacency map
  //cout<<endl<<"Adjacency map"<<endl;
  //cout<<"size = "<<adjacencyMap.size()<<endl;
  //for (std::map<Edge, TriangleIDVector>::const_iterator it = adjacencyMap.begin(); it != adjacencyMap.end(); ++it)
  //{
  //cout<<"Edge "<<it->first;
  //for (TriangleIDVector::const_iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2)
  //cout<<it2->getIndex()<<" ";
  //cout<<endl;
  //}
  
  // Building the edge points map from the adjacent map
  std::map<Edge, TriangleIDVector>::const_iterator it;
  for (it = adjacencyMap.begin(); it != adjacencyMap.end(); ++it)
  {
  Edge ed = it->first;
      
  TriangleIDVector triangleVec = it->second;
  int nbAdjacentTriangles = triangleVec.size();

  assert((nbAdjacentTriangles > 0) && (nbAdjacentTriangles <= 2));

  if (nbAdjacentTriangles == 1)
  {
  // Inserting the points of the corresponding edge in the Edge Points Map
  GLuint p1, p2;
  qglviewer::Vec pos1, pos2;
  p1 = ed.first;
  p2 = ed.second;

  assert((p1 <= model->numvertices) && (p2 <= model->numvertices));

  for (int k = 0; k != 3; ++k)
  {
  pos1[k] = model->vertices[3*p1 + k];
  pos2[k] = model->vertices[3*p2 + k];
  }

  edgePointsMap[p1] = pos1;
  edgePointsMap[p2] = pos2;
  }
  }
  }
*/


/* glmList: Generates and returns a display list for the model using
 * the mode specified.
 *
 * model - initialized GLMmodel structure
 * mode  - a bitwise OR of values describing what is to be rendered.
 *             GLM_NONE     -  render with only vertices
 *             GLM_FLAT     -  render with facet normals
 *             GLM_SMOOTH   -  render with vertex normals
 *             GLM_TEXTURE  -  render with texture coords
 *             GLM_COLOR    -  render with colors (color material)
 *             GLM_MATERIAL -  render with materials
 *             GLM_COLOR and GLM_MATERIAL should not both be specified.  
 * GLM_FLAT and GLM_SMOOTH should not both be specified.  

 */
GLuint
glmList(GLMmodel* model, GLuint mode)
{
  // AJOUT FLO
  float color[3] = {1.0f, 1.0f, 1.0f};
  //
  
  GLuint list;
    
  list = glGenLists(1);
  glNewList(list, GL_COMPILE);
  glmDraw(model, mode, color);
  glEndList();
    
  return list;
}

/* glmWeld: eliminate (weld) vectors that are within an epsilon of
 * each other.
 *
 * model   - initialized GLMmodel structure
 * epsilon     - maximum difference between vertices
 *               ( 0.00001 is a good start for a unitized model)
 *
 */
GLvoid
glmWeld(GLMmodel* model, GLfloat epsilon)
{
    GLfloat* vectors;
    GLfloat* copies;
    GLuint   numvectors;
    GLuint   i;
    
    /* vertices */
    numvectors = model->numvertices;
    vectors  = model->vertices;
    copies = glmWeldVectors(vectors, &numvectors, epsilon);
    
#if 1
    printf("glmWeld(): %d redundant vertices.\n", 
        model->numvertices - numvectors - 1);
#endif
    
    for (i = 0; i < model->numtriangles; i++) {
        T(i).vindices[0] = (GLuint)vectors[3 * T(i).vindices[0] + 0];
        T(i).vindices[1] = (GLuint)vectors[3 * T(i).vindices[1] + 0];
        T(i).vindices[2] = (GLuint)vectors[3 * T(i).vindices[2] + 0];
    }
    
    /* free space for old vertices */
    free(vectors);
    
    /* allocate space for the new vertices */
    model->numvertices = numvectors;
    model->vertices = (GLfloat*)malloc(sizeof(GLfloat) * 
        3 * (model->numvertices + 1));
    
    /* copy the optimized vertices into the actual vertex list */
    for (i = 1; i <= model->numvertices; i++) {
        model->vertices[3 * i + 0] = copies[3 * i + 0];
        model->vertices[3 * i + 1] = copies[3 * i + 1];
        model->vertices[3 * i + 2] = copies[3 * i + 2];
    }
    
    free(copies);
}

/* glmReadPPM: read a PPM raw (type P6) file.  The PPM file has a header

 * that should look something like:

 *

 *    P6

 *    # comment

 *    width height max_value

 *    rgbrgbrgb...

 *

 * where "P6" is the magic cookie which identifies the file type and

 * should be the only characters on the first line followed by a

 * carriage return.  Any line starting with a # mark will be treated

 * as a comment and discarded.   After the magic cookie, three integer

 * values are expected: width, height of the image and the maximum

 * value for a pixel (max_value must be < 256 for PPM raw files).  The

 * data section consists of width*height rgb triplets (one byte each)

 * in binary format (i.e., such as that written with fwrite() or

 * equivalent).

 *

 * The rgb data is returned as an array of unsigned chars (packed

 * rgb).  The malloc()'d memory should be free()'d by the caller.  If

 * an error occurs, an error message is sent to stderr and NULL is

 * returned.

 *

 * filename   - name of the .ppm file.

 * width      - will contain the width of the image on return.

 * height     - will contain the height of the image on return.

 *

 */

GLubyte* 

glmReadPPM(char* filename, int* width, int* height)

{

    FILE* fp;

    int i, w, h, d;

    unsigned char* image;

    char head[70];          /* max line <= 70 in PPM (per spec). */

    

    fp = fopen(filename, "rb");

    if (!fp) {

        perror(filename);

        return NULL;

    }

    

    /* grab first two chars of the file and make sure that it has the

       correct magic cookie for a raw PPM file. */

    fgets(head, 70, fp);

    if (strncmp(head, "P6", 2)) {

        fprintf(stderr, "%s: Not a raw PPM file\n", filename);

        return NULL;

    }

    

    /* grab the three elements in the header (width, height, maxval). */

    i = 0;

    while(i < 3) {

        fgets(head, 70, fp);

        if (head[0] == '#')     /* skip comments. */

            continue;

        if (i == 0)

            i += sscanf(head, "%d %d %d", &w, &h, &d);

        else if (i == 1)

            i += sscanf(head, "%d %d", &h, &d);

        else if (i == 2)

            i += sscanf(head, "%d", &d);

    }

    

    /* grab all the image data in one fell swoop. */

    image = (unsigned char*)malloc(sizeof(unsigned char)*w*h*3);

    fread(image, sizeof(unsigned char), w*h*3, fp);

    fclose(fp);

    

    *width = w;

    *height = h;

    return image;

}


#if 0
/* normals */
if (model->numnormals) {
    numvectors = model->numnormals;
    vectors  = model->normals;
    copies = glmOptimizeVectors(vectors, &numvectors);
    
    printf("glmOptimize(): %d redundant normals.\n", 
        model->numnormals - numvectors);
    
    for (i = 0; i < model->numtriangles; i++) {
        T(i).nindices[0] = (GLuint)vectors[3 * T(i).nindices[0] + 0];
        T(i).nindices[1] = (GLuint)vectors[3 * T(i).nindices[1] + 0];
        T(i).nindices[2] = (GLuint)vectors[3 * T(i).nindices[2] + 0];
    }
    
    /* free space for old normals */
    free(vectors);
    
    /* allocate space for the new normals */
    model->numnormals = numvectors;
    model->normals = (GLfloat*)malloc(sizeof(GLfloat) * 
        3 * (model->numnormals + 1));
    
    /* copy the optimized vertices into the actual vertex list */
    for (i = 1; i <= model->numnormals; i++) {
        model->normals[3 * i + 0] = copies[3 * i + 0];
        model->normals[3 * i + 1] = copies[3 * i + 1];
        model->normals[3 * i + 2] = copies[3 * i + 2];
    }
    
    free(copies);
}

/* texcoords */
if (model->numtexcoords) {
    numvectors = model->numtexcoords;
    vectors  = model->texcoords;
    copies = glmOptimizeVectors(vectors, &numvectors);
    
    printf("glmOptimize(): %d redundant texcoords.\n", 
        model->numtexcoords - numvectors);
    
    for (i = 0; i < model->numtriangles; i++) {
        for (j = 0; j < 3; j++) {
            T(i).tindices[j] = (GLuint)vectors[3 * T(i).tindices[j] + 0];
        }
    }
    
    /* free space for old texcoords */
    free(vectors);
    
    /* allocate space for the new texcoords */
    model->numtexcoords = numvectors;
    model->texcoords = (GLfloat*)malloc(sizeof(GLfloat) * 
        2 * (model->numtexcoords + 1));
    
    /* copy the optimized vertices into the actual vertex list */
    for (i = 1; i <= model->numtexcoords; i++) {
        model->texcoords[2 * i + 0] = copies[2 * i + 0];
        model->texcoords[2 * i + 1] = copies[2 * i + 1];
    }
    
    free(copies);
}
#endif

#if 0
/* look for unused vertices */
/* look for unused normals */
/* look for unused texcoords */
for (i = 1; i <= model->numvertices; i++) {
    for (j = 0; j < model->numtriangles; i++) {
        if (T(j).vindices[0] == i || 
            T(j).vindices[1] == i || 
            T(j).vindices[1] == i)
            break;
    }
}
#endif
