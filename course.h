/* --------------------------------------------------------------------
EXTREME TUXRACER

Copyright (C) 1999-2001 Jasmin F. Patry (Tuxracer)
Copyright (C) 2010 Extreme Tuxracer Team

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
---------------------------------------------------------------------*/

#ifndef COURSE_H
#define COURSE_H

#include "bh.h"
#include <vector>
#include <map>

#define FLOATVAL(i) (*(GLfloat*)(vnc_array+idx+(i)*sizeof(GLfloat)))
#ifdef USE_GLES1
#define BYTEVAL(i) (*(GLubyte*)(vnc_array+idx+10*sizeof(GLfloat) +\
    i*sizeof(GLubyte)))
#define STRIDE_GL_ARRAY (10 * sizeof(GLfloat) + 4 * sizeof(GLubyte))
#else
#define BYTEVAL(i) (*(GLubyte*)(vnc_array+idx+8*sizeof(GLfloat) +\
    i*sizeof(GLubyte)))
#define STRIDE_GL_ARRAY (8 * sizeof(GLfloat) + 4 * sizeof(GLubyte))
#endif
#define ELEV(x,y) (elevation[(x) + nx*(y)] )
#define NORM_INTERPOL 0.05
#define XCD(x) ((ETR_DOUBLE)(x) / (nx-1.0) * curr_course->width)
#define ZCD(y) (-(ETR_DOUBLE)(y) / (ny-1.0) * curr_course->length)
#define NMLPOINT(x,y) TVector3(XCD(x), ELEV(x,y), ZCD(y) )


#define MAX_COURSES 64
#define MAX_TERR_TYPES 64
#define MAX_OBJECT_TYPES 128
#define MAX_DESCRIPTION_LINES 8

class TTexture;

struct TCollidable {
	TVector3 pt;
	ETR_DOUBLE height;
	ETR_DOUBLE diam;
	size_t tree_type;
};

struct TItem {
	TVector3 pt;
	ETR_DOUBLE height;
	ETR_DOUBLE diam;
	size_t item_type;
	int collectable;
	bool drawable;
};

struct TCourse {
	string name;
	string dir;
	string author;
	string desc[MAX_DESCRIPTION_LINES];
	size_t num_lines;
	TTexture* preview;
	ETR_DOUBLE width;
	ETR_DOUBLE length;
	ETR_DOUBLE play_width;
	ETR_DOUBLE play_length;
	ETR_DOUBLE angle;
	ETR_DOUBLE scale;
	ETR_DOUBLE startx;
	ETR_DOUBLE starty;
	size_t env;
	size_t music_theme;
 	bool use_keyframe;
	ETR_DOUBLE finish_brake;
};

struct TTerrType {
	string textureFile;
	TTexture* texture;
	string sound;
	TColor3 col;

	ETR_DOUBLE friction;
	ETR_DOUBLE depth;
	int vol_type;
	int particles;
	int trackmarks;
	int starttex;
	int tracktex;
	int stoptex;
	bool shiny;
};

struct TObjectType {
	string		name;
	string		textureFile;
	TTexture*	texture;
	bool		collidable;
    int			collectable;
    bool		drawable;
    bool		reset_point;
    bool		use_normal;
    TVector3	normal;
    int			num_items;
    int			poly;
};

class CCourse {
private:
	const TCourse* curr_course;
	map<string, size_t> CourseIndex;
	map<string, size_t> ObjectIndex;
	string		CourseDir;
	bool		SaveItemsFlag;

	int			nx;
	int			ny;
	TVector2	start_pt;
	int			base_height_value;
	bool		mirrored;

	void		FreeTerrainTextures ();
	void		FreeObjectTextures ();
	void		CalcNormals ();
	void		MakeCourseNormals ();
	bool		LoadElevMap ();
	void		LoadItemList ();
	bool		LoadObjectMap ();
	bool		LoadTerrainMap ();
	int			GetTerrain (unsigned char pixel[]) const;

	void		MirrorCourseData ();
public:
	CCourse ();
	~CCourse();

	vector<TCourse>		CourseList;
	vector<TTerrType>	TerrList;
	vector<TObjectType>	ObjTypes;
	vector<TCollidable>	CollArr;
	vector<TItem>		NocollArr;
	vector<TPolyhedron>	PolyArr;

	char		*terrain;
	ETR_DOUBLE		*elevation;
	TVector3	*nmls;
	GLubyte		*vnc_array;

	void ResetCourse ();
 	size_t GetCourseIdx (const string& dir) const;
	bool LoadCourseList ();
	void FreeCourseList ();
	bool LoadCourse (size_t idx);
	bool LoadTerrainTypes ();
	bool LoadObjectTypes ();
	void MakeStandardPolyhedrons ();
	void GetGLArrays (GLubyte **vnc_array) const;
	void FillGlArrays();

	void GetDimensions (ETR_DOUBLE *w, ETR_DOUBLE *l) const;
	void GetPlayDimensions (ETR_DOUBLE *pw, ETR_DOUBLE *pl) const;
	void GetDivisions (int *nx, int *ny) const;
	ETR_DOUBLE GetCourseAngle () const;
	ETR_DOUBLE GetBaseHeight (ETR_DOUBLE distance) const;
	ETR_DOUBLE GetMaxHeight (ETR_DOUBLE distance) const;
	size_t GetEnv () const;
	const TVector2& GetStartPoint () const;
	const TPolyhedron& GetPoly (size_t type) const;
	void MirrorCourse ();

	void GetIndicesForPoint (ETR_DOUBLE x, ETR_DOUBLE z, int *x0, int *y0, int *x1, int *y1) const;
	void FindBarycentricCoords (ETR_DOUBLE x, ETR_DOUBLE z,
		TIndex2 *idx0, TIndex2 *idx1, TIndex2 *idx2, ETR_DOUBLE *u, ETR_DOUBLE *v) const;
	TVector3 FindCourseNormal (ETR_DOUBLE x, ETR_DOUBLE z) const;
	ETR_DOUBLE FindYCoord (ETR_DOUBLE x, ETR_DOUBLE z) const;
	void GetSurfaceType (ETR_DOUBLE x, ETR_DOUBLE z, ETR_DOUBLE weights[]) const;
	int GetTerrainIdx (ETR_DOUBLE x, ETR_DOUBLE z, ETR_DOUBLE level) const;
	TPlane GetLocalCoursePlane (TVector3 pt) const;
};

extern CCourse Course;

#endif
