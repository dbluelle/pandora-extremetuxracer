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
#include "mathlib.h"

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
#define MAX_FRICT_FORCE 800.0f
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

#define MAX_ROLL_ANGLE 30
#define BRAKING_ROLL_ANGLE 55

// constants for finish stage
#define FIN_AIR_GRAV 500
#define FIN_GRAV 500
#define FIN_AIR_BRAKE 20
#define FIN_BRAKE 12

struct TForce {
	TVector3d surfnml;
	TVector3d rollnml;
	TVector3d pos;
	TVector3d vel;
	TVector3d frictdir;

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

	bool     CheckTreeCollisions (const TVector3d& pos, TVector3d *tree_loc, ETR_DOUBLE *tree_diam);
	void     AdjustTreeCollision (const TVector3d& pos, TVector3d *vel);
	void     CheckItemCollection (const TVector3d& pos);

	TVector3d CalcRollNormal (ETR_DOUBLE speed);
	TVector3d CalcAirForce ();
	TVector3d CalcSpringForce ();
	TVector3d CalcNormalForce ();
	TVector3d CalcJumpForce ();
	TVector3d CalcFrictionForce (ETR_DOUBLE speed, const TVector3d& nmlforce);
	TVector3d CalcPaddleForce (ETR_DOUBLE speed);
	TVector3d CalcBrakeForce (ETR_DOUBLE speed);
	TVector3d CalcGravitationForce ();
	TVector3d CalcNetForce (const TVector3d& pos, const TVector3d& vel);
	TVector3d CalcFinishForce (const TVector3d& pos, const TVector3d& vel);

	void     AdjustVelocity (const TPlane& surf_plane);
	void     AdjustPosition (const TPlane& surf_plane, ETR_DOUBLE dist_from_surface);
	void     SetTuxPosition (ETR_DOUBLE speed);
	ETR_DOUBLE   AdjustTimeStep (ETR_DOUBLE h, const TVector3d& vel);
	void     SolveOdeSystem ();
public:
	CControl ();

	// view:
	TViewMode viewmode;
	TVector3d viewpos;
	TVector3d plyr_pos;
	TVector3d viewdir;
	TVector3d viewup;
	TMatrix<4, 4> view_mat;
	bool view_init;
	// main:
	TVector3d cpos;
	TVector3d cvel;
	TVector3d last_pos;
	TVector3d cnet_force;
	TVector3d cdirection;
	TQuaternion corientation;
	ETR_DOUBLE way;

	bool orientation_initialized;
	TVector3d plane_nml;
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
	void UpdatePlayerPos (bool eps);
};

#endif
