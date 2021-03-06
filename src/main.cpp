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

#ifdef HAVE_CONFIG_H
#include <etr_config.h>
#endif

#include "bh.h"
#include "textures.h"
#include "ogl.h"
#include "splash_screen.h"
#include "audio.h"
#include "font.h"
#include "tools.h"
#include "ogl_test.h"
#include "winsys.h"
#include <iostream>
#include <ctime>

TGameData g_game;

void InitGame (int argc, char **argv) {
	g_game.toolmode = NONE;
	g_game.argument = 0;
	if (argc == 4) {
		string group_arg = argv[1];
		if (group_arg == "--char") g_game.argument = 4;
		Tools.SetParameter(argv[2], argv[3]);
	} else if (argc == 2) {
		string group_arg = argv[1];
		if (group_arg == "9") g_game.argument = 9;
	}

	g_game.player = NULL;
	g_game.start_player = 0;
	g_game.course = NULL;
	g_game.mirrorred = false;
	g_game.character = NULL;
	g_game.location_id = 0;
	g_game.light_id = 0;
	g_game.snow_id = 0;
	g_game.cup = 0;
	g_game.theme_id = 0;
	g_game.force_treemap = 0;
	g_game.treesize = 3;
	g_game.treevar = 3;
}

// ====================================================================
// 					main
// ====================================================================

#if defined ( OS_WIN32_MINGW )
#undef main
#endif

#ifdef PANDORA
void enable_fastmath()
{
	static const unsigned int x = 0x04086060;
	static const unsigned int y = 0x03000000;
	int r;
	asm volatile (
		"fmrx	%0, fpscr			\n\t"	//r0 = FPSCR
		"and	%0, %0, %1			\n\t"	//r0 = r0 & 0x04086060
		"orr	%0, %0, %2			\n\t"	//r0 = r0 | 0x03000000
		"fmxr	fpscr, %0			\n\t"	//FPSCR = r0
		: "=r"(r)
		: "r"(x), "r"(y)
	);
}
#endif
int main( int argc, char **argv ) {
	// ****************************************************************
	cout << "\n----------- Extreme Tux Racer " ETR_VERSION_STRING " ----------------";
	cout << "\n----------- (C) 2010-2013 Extreme Tuxracer Team  --------\n\n";
#ifdef PANDORA
	enable_fastmath();
#endif
	srand (time (NULL));
	InitConfig (argv[0]);
	InitGame (argc, argv);
	Winsys.Init ();
	InitOpenglExtensions ();
	BuildGlobalVBO();

	// for checking the joystick and the OpgenGL version (the info is
	// written on the console):
	//	Winsys.PrintJoystickInfo ();
	//	PrintGLInfo ();

	// theses resources must or should be loaded before splashscreen starts
	Tex.LoadTextureList ();
	FT.LoadFontlist ();
	Winsys.SetFonttype ();
	Audio.Open ();
	Music.LoadMusicList ();
	Music.SetVolume (param.music_volume);

	switch (g_game.argument) {
		case 0:
			State::manager.Run(SplashScreen);
			break;
		case 4:
			g_game.toolmode = TUXSHAPE;
			State::manager.Run(Tools);
			break;
		case 9:
			State::manager.Run(OglTest);
			break;
	}

	Winsys.Quit();
	DeleteGlobalVBO();
#ifdef USE_GLES1
	closegles();
#endif
	return 0;
}
