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

#ifndef PHYSICS_H
#define PHYSICS_H

#include "bh.h"

#define MAX_PADDLING_SPEED (60.0 / 3.6)   /* original 60 */
#define PADDLE_FACT 1.0 /* original 1.0 */

#define EARTH_GRAV 9.81  /* ok, why not 10.0 ? */
#define JUMP_FORCE_DURATION 0.20
#define TUX_MASS 20
#define MIN_TUX_SPEED 1.4
#define INIT_TUX_SPEED 3.0
#define COLL_TOLERANCE 0.1

#define MAX_SURF_PEN 0.2
#define TUX_Y_CORR 0.36
#define IDEAL_ROLL_SPEED 6.0
#define IDEAL_ROLL_FRIC 0.35
#define WIND_FACTOR 1.5

#define MIN_FRICT_SPEED 2.8
#define MAX_FRICT_FORCE 800
#define MAX_TURN_ANGLE 45
#define MAX_TURN_PERP 400
#define MAX_TURN_PEN 0.15
#define PADDLING_DURATION 0.40
#define IDEAL_PADD_FRIC 0.35
#define MAX_PADD_FORCE 122.5
#define BRAKE_FORCE 200

#define MIN_TIME_STEP 0.01
#define MAX_TIME_STEP 0.10
#define MAX_STEP_DIST 0.20
#define MAX_POS_ERR 0.005
#define MAX_VEL_ERR	0.05

// constants for finish stage
#define FIN_AIR_GRAV 500
#define FIN_GRAV 500
#define FIN_AIR_BRAKE 20
#define FIN_BRAKE 12

struct TForce {
	TVector3 surfnml;
	TVector3 rollnml;
	TVector3 pos;
	TVector3 vel;
	TVector3 frictdir;

	ETR_DOUBLE frict_coeff;
	ETR_DOUBLE comp_depth;
	ETR_DOUBLE surfdistance;
    ETR_DOUBLE compression;
};

class CControl {
private:
	TForce ff;
	ETR_DOUBLE ode_time_step;
	ETR_DOUBLE finish_speed;

	bool     CheckTreeCollisions (const TVector3& pos, TVector3 *tree_loc, ETR_DOUBLE *tree_diam);
	void     AdjustTreeCollision (const TVector3& pos, TVector3 *vel);
	void     CheckItemCollection (const TVector3& pos);

	TVector3 CalcRollNormal (ETR_DOUBLE speed);
	TVector3 CalcAirForce ();
	TVector3 CalcSpringForce ();
	TVector3 CalcNormalForce ();
	TVector3 CalcJumpForce ();
	TVector3 CalcFrictionForce (ETR_DOUBLE speed, const TVector3& nmlforce);
	TVector3 CalcPaddleForce (ETR_DOUBLE speed);
	TVector3 CalcBrakeForce (ETR_DOUBLE speed);
	TVector3 CalcGravitationForce ();
	TVector3 CalcNetForce (const TVector3& pos, const TVector3& vel);
	TVector3 CalcFinishForce (const TVector3& pos, const TVector3& vel);

	void     AdjustVelocity (const TPlane& surf_plane);
	void     AdjustPosition (const TPlane& surf_plane, ETR_DOUBLE dist_from_surface);
	void     SetTuxPosition (ETR_DOUBLE speed);
	ETR_DOUBLE   AdjustTimeStep (ETR_DOUBLE h, TVector3 vel);
	void     SolveOdeSystem (ETR_DOUBLE timestep);
public:
	CControl ();

	// view:
	TViewMode viewmode;
    TVector3 viewpos;
    TVector3 plyr_pos;
    TVector3 viewdir;
    TVector3 viewup;
    TMatrix view_mat;
    bool view_init;
	// main:
	TVector3 cpos;
    TVector3 cvel;
	TVector3 last_pos;
    TVector3 cnet_force;
    TVector3 cdirection;
	TQuaternion corientation;
	ETR_DOUBLE way;

    bool orientation_initialized;
    TVector3 plane_nml;
	// steering:
    ETR_DOUBLE turn_fact;
    ETR_DOUBLE turn_animation;
	ETR_DOUBLE paddle_time;
    ETR_DOUBLE jump_amt;
    ETR_DOUBLE jump_start_time;
    bool   is_paddling;
    bool   is_braking;
    bool   begin_jump;
    bool   jumping;
    bool   jump_charging;
	// trick:
    bool   front_flip;
    bool   back_flip;
    bool   cairborne;
    bool   roll_left;
    bool   roll_right;
    ETR_DOUBLE roll_factor;
    ETR_DOUBLE flip_factor;
	// pseudo constants:
	ETR_DOUBLE minSpeed;
	ETR_DOUBLE minFrictspeed;

	void Init ();
	void UpdatePlayerPos (ETR_DOUBLE timestep);
};

// STRUCTURE OF PHYSICS.CPP

// ------------------- init -------------------------------------------
// void InitSimulation ()

// ----------------- collision ----------------------------------------
// bool  CheckTreeCollisions (TVector3 pos, TVector3 *tree_loc, ETR_DOUBLE *tree_diam)
// void  AdjustTreeCollision (TVector3 pos, TVector3 *vel)
// void  CheckItemCollection (TVector3 pos)

// ----------------- position and vlelocity ---------------------------
// ETR_DOUBLE AdjustVelocity (TVector3 *vel, TVector3 pos, TPlane surf_plane,
// 	   ETR_DOUBLE dist_from_surface)
// void   AdjustPosition (TVector3 *pos, TPlane surf_plane, ETR_DOUBLE dist_from_surface)
// void   SetTuxPosition (TVector3 new_pos)

// ----------------- forces -------------------------------------------
// TVector3 AdjustRollNormal (TVector3 vel, ETR_DOUBLE fric_coeff, TVector3 nml)
// TVector3 CalcAirForce (TVector3 player_vel)
// TVector3 CalcSpringForce (ETR_DOUBLE compression, TVector3 vel, TVector3 surf_nml)
// TVector3 CalcNetForce (TVector3 pos, TVector3 vel)

// ------------------ ODE solver and port function
// static ETR_DOUBLE adjust_time_step_size (ETR_DOUBLE h, TVector3 vel)
// void   SolveOdeSystem (ETR_DOUBLE timestep)
// void   UpdatePlayerPos (ETR_DOUBLE timestep)

#endif
