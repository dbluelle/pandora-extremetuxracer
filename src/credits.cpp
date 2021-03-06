/* --------------------------------------------------------------------
EXTREME TUXRACER

Copyright (C) 2010 Extreme Tuxracer Team

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
---------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <etr_config.h>
#endif

#include "credits.h"
#include "audio.h"
#include "ogl.h"
#include "particles.h"
#include "textures.h"
#include "font.h"
#include "gui.h"
#include "spx.h"
#include "game_type_select.h"
#include "winsys.h"

#define TOP_Y 160
#define BOTT_Y 64
#define OFFS_SCALE_FACTOR 1.2

CCredits Credits;


static ETR_DOUBLE y_offset = 0;
static bool moving = true;

void CCredits::LoadCreditList () {
	CSPList list(MAX_CREDITS);

	if (!list.Load (param.data_dir, "credits.lst")) {
		Message ("could not load credits list");
		return;
	}

	for (size_t i=0; i<list.Count(); i++) {
		const string& line = list.Line(i);
		TCredits credit;
		credit.text = SPStrN (line, "text", "");

		ETR_DOUBLE offset = SPFloatN (line, "offs", 0) * OFFS_SCALE_FACTOR * Winsys.scale;
		if (i>0) credit.offs = CreditList.back().offs + (int)offset;
		else credit.offs = offset;

		credit.col = SPIntN (line, "col", 0);
		credit.size = SPFloatN (line, "size", 1.0);
		CreditList.push_back(credit);
	}
}

void CCredits::DrawCreditsText () {
	int w = Winsys.resolution.width;
	int h = Winsys.resolution.height;
	ETR_DOUBLE offs = 0.0;
	if (moving) y_offset += g_game.time_step * 30;


	for (list<TCredits>::const_iterator i = CreditList.begin(); i != CreditList.end(); ++i) {
		offs = h - 100 - y_offset + i->offs;
		if (offs > h || offs < 0.0) // Draw only visible lines
			continue;

		if (i->col == 0)
			FT.SetColor (colWhite);
		else
			FT.SetColor (colDYell);
		FT.AutoSizeN (i->size);
		FT.DrawString (-1, (int)offs, i->text);
	}


	glDisable (GL_TEXTURE_2D);
	glColor(colBackgr);
	glRecti (0, 0, w, BOTT_Y);

	glBegin( GL_QUADS );
	glVertex2i(0, BOTT_Y);
	glVertex2i(w, BOTT_Y);
	glColor(colBackgr, 0);
	glVertex2i(w, BOTT_Y + 30);
	glVertex2i(0, BOTT_Y + 30);
	glEnd();

	glColor(colBackgr);
	glRecti (0, h - TOP_Y, w, h);

	glBegin( GL_QUADS );
	glVertex2i(w, h - TOP_Y);
	glVertex2i(0, h - TOP_Y);
	glColor(colBackgr, 0);
	glVertex2i(0, h - TOP_Y - 30);
	glVertex2i (w, h - TOP_Y - 30);
	glEnd();

	glEnable (GL_TEXTURE_2D);
	if (offs < TOP_Y) y_offset = 0;
}

void CCredits::Keyb (unsigned int key, bool special, bool release, int x, int y) {
	if (release) return;
	switch (key) {
		case SDLK_m:
			moving = !moving;
			break;
		case SDLK_u:
			param.ui_snow = !param.ui_snow;
			break;
		default:
			State::manager.RequestEnterState (GameTypeSelect);
	}
}

void CCredits::Mouse (int button, int state, int x, int y) {
	if (state == 1) State::manager.RequestEnterState (GameTypeSelect);
}

void CCredits::Motion(int x, int y) {
	if (param.ui_snow) push_ui_snow (cursor_pos);
}

void CCredits::Enter() {
	Music.Play (param.credits_music, -1);
	y_offset = 0;
	moving = true;
}

void CCredits::Loop() {
	int ww = Winsys.resolution.width;
	int hh = Winsys.resolution.height;

	Music.Update ();
	check_gl_error();
	ClearRenderContext ();
	ScopedRenderMode rm(GUI);
	SetupGuiDisplay ();

	DrawCreditsText ();
	if (param.ui_snow) {
		update_ui_snow ();
		draw_ui_snow();
	}
	Tex.Draw (BOTTOM_LEFT, 0, hh-256, 1);
	Tex.Draw (BOTTOM_RIGHT, ww-256, hh-256, 1);
	Tex.Draw (TOP_LEFT, 0, 0, 1);
	Tex.Draw (TOP_RIGHT, ww-256, 0, 1);
	Tex.Draw (T_TITLE_SMALL, CENTER, AutoYPosN (5), Winsys.scale);


	Reshape (ww, hh);
	Winsys.SwapBuffers();
}
