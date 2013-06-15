/* --------------------------------------------------------------------
EXTREME TUXRACER

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

#ifndef GAME_CONFIG_SCREEN_H
#define GAME_CONFIG_SCREEN_H

#include "bh.h"
#include "states.h"

class CGameConfig : public State {
	void Enter();
	void Loop(ETR_DOUBLE time_step);
	void Keyb(unsigned int key, bool special, bool release, int x, int y);
	void Mouse(int button, int state, int x, int y);
	void Motion(int x, int y);
	void Exit();
public:
};

extern CGameConfig GameConfig;

#endif
