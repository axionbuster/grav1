#include "Include.h"

int wWinMain(void* _0, void* _1, void* _2, int _3)
{
	// Raylib.
	InitWindow(600, 600, "a");
	SetTargetFPS(60);

	while (!WindowShouldClose())
	{
		BeginDrawing();
		{
			ClearBackground(WHITE);
		}
		EndDrawing();
	}

	CloseWindow();
	return 0;
}
