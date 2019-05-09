/*
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
 */

/** \file
 * \ingroup GHOST
 */

#include "GHOST_WindowSDL.h"
#include "SDL_mouse.h"
#include "glew-mx.h"

#include "GHOST_ContextSDL.h"

#include <assert.h>

GHOST_WindowSDL::GHOST_WindowSDL(GHOST_SystemSDL *system,
                                 const STR_String &title,
                                 GHOST_TInt32 left,
                                 GHOST_TInt32 top,
                                 GHOST_TUns32 width,
                                 GHOST_TUns32 height,
                                 GHOST_TWindowState state,
                                 const GHOST_TEmbedderWindowID parentWindow,
                                 GHOST_TDrawingContextType type,
                                 const bool stereoVisual,
                                 const bool exclusive)
    : GHOST_Window(width, height, state, stereoVisual, exclusive),
      m_system(system),
      m_valid_setup(false),
      m_invalid_window(false),
      m_sdl_custom_cursor(NULL)
{

  /* creating the window _must_ come after setting attributes */
  m_sdl_win = SDL_CreateWindow(title,
                               left,
                               top,
                               width,
                               height,
                               SDL_WINDOW_RESIZABLE | SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);

  /* now set up the rendering context. */
  if (setDrawingContextType(type) == GHOST_kSuccess) {
    m_valid_setup = true;
    GHOST_PRINT("Created window\n");
  }

  if (exclusive) {
    SDL_RaiseWindow(m_sdl_win);
  }

  setTitle(title);
}

GHOST_WindowSDL::~GHOST_WindowSDL()
{
  if (m_sdl_custom_cursor) {
    SDL_FreeCursor(m_sdl_custom_cursor);
  }

  releaseNativeHandles();

  SDL_DestroyWindow(m_sdl_win);
}

GHOST_Context *GHOST_WindowSDL::newDrawingContext(GHOST_TDrawingContextType type)
{
  if (type == GHOST_kDrawingContextTypeOpenGL) {
    GHOST_Context *context = new GHOST_ContextSDL(m_wantStereoVisual,
                                                  m_sdl_win,
                                                  0,  // profile bit
                                                  3,
                                                  3,
                                                  GHOST_OPENGL_SDL_CONTEXT_FLAGS,
                                                  GHOST_OPENGL_SDL_RESET_NOTIFICATION_STRATEGY);

    if (context->initializeDrawingContext())
      return context;
    else
      delete context;
  }

  return NULL;
}

GHOST_TSuccess GHOST_WindowSDL::invalidate(void)
{
  if (m_invalid_window == false) {
    m_system->addDirtyWindow(this);
    m_invalid_window = true;
  }

  return GHOST_kSuccess;
}

GHOST_TSuccess GHOST_WindowSDL::setState(GHOST_TWindowState state)
{
  switch (state) {
    case GHOST_kWindowStateNormal:
      SDL_SetWindowFullscreen(m_sdl_win, SDL_FALSE);
      SDL_RestoreWindow(m_sdl_win);
      break;
    case GHOST_kWindowStateMaximized:
      SDL_SetWindowFullscreen(m_sdl_win, SDL_FALSE);
      SDL_MaximizeWindow(m_sdl_win);
      break;
    case GHOST_kWindowStateMinimized:
      SDL_MinimizeWindow(m_sdl_win);
      break;
    case GHOST_kWindowStateFullScreen:
      SDL_SetWindowFullscreen(m_sdl_win, SDL_TRUE);
      break;
    default:
      break;
  }

  return GHOST_kSuccess;
}

GHOST_TWindowState GHOST_WindowSDL::getState() const
{
  Uint32 flags = SDL_GetWindowFlags(m_sdl_win);

  if (flags & SDL_WINDOW_FULLSCREEN)
    return GHOST_kWindowStateFullScreen;
  else if (flags & SDL_WINDOW_MAXIMIZED)
    return GHOST_kWindowStateMaximized;
  else if (flags & SDL_WINDOW_MINIMIZED)
    return GHOST_kWindowStateMinimized;
  return GHOST_kWindowStateNormal;
}

bool GHOST_WindowSDL::getValid() const
{
  return GHOST_Window::getValid() && m_valid_setup;
}

void GHOST_WindowSDL::setTitle(const STR_String &title)
{
  SDL_SetWindowTitle(m_sdl_win, title.ReadPtr());
}

void GHOST_WindowSDL::getTitle(STR_String &title) const
{
  title = SDL_GetWindowTitle(m_sdl_win);
}

void GHOST_WindowSDL::getWindowBounds(GHOST_Rect &bounds) const
{
  getClientBounds(bounds);
}

void GHOST_WindowSDL::getClientBounds(GHOST_Rect &bounds) const
{
  int x, y, w, h;
  SDL_GetWindowSize(m_sdl_win, &w, &h);
  SDL_GetWindowPosition(m_sdl_win, &x, &y);

  bounds.m_l = x;
  bounds.m_r = x + w;
  bounds.m_t = y;
  bounds.m_b = y + h;
}

GHOST_TSuccess GHOST_WindowSDL::setClientWidth(GHOST_TUns32 width)
{
  int height;
  SDL_GetWindowSize(m_sdl_win, NULL, &height);
  SDL_SetWindowSize(m_sdl_win, width, height);
  return GHOST_kSuccess;
}

GHOST_TSuccess GHOST_WindowSDL::setClientHeight(GHOST_TUns32 height)
{
  int width;
  SDL_GetWindowSize(m_sdl_win, &width, NULL);
  SDL_SetWindowSize(m_sdl_win, width, height);
  return GHOST_kSuccess;
}

GHOST_TSuccess GHOST_WindowSDL::setClientSize(GHOST_TUns32 width, GHOST_TUns32 height)
{
  SDL_SetWindowSize(m_sdl_win, width, height);
  return GHOST_kSuccess;
}

void GHOST_WindowSDL::screenToClient(GHOST_TInt32 inX,
                                     GHOST_TInt32 inY,
                                     GHOST_TInt32 &outX,
                                     GHOST_TInt32 &outY) const
{
  /* XXXSDL_WEAK_ABS_COORDS */
  int x_win, y_win;
  SDL_GetWindowPosition(m_sdl_win, &x_win, &y_win);

  outX = inX - x_win;
  outY = inY - y_win;
}
void GHOST_WindowSDL::clientToScreen(GHOST_TInt32 inX,
                                     GHOST_TInt32 inY,
                                     GHOST_TInt32 &outX,
                                     GHOST_TInt32 &outY) const
{
  /* XXXSDL_WEAK_ABS_COORDS */
  int x_win, y_win;
  SDL_GetWindowPosition(m_sdl_win, &x_win, &y_win);

  outX = inX + x_win;
  outY = inY + y_win;
}

/* mouse cursor */
static unsigned char sdl_std_cursor_mask_xterm[] = {
    0xef, 0x01, 0xff, 0x01, 0xff, 0x01, 0x7c, 0x00, 0x38, 0x00, 0x38, 0x00, 0x38, 0x00, 0x38, 0x00,
    0x38, 0x00, 0x38, 0x00, 0x38, 0x00, 0x38, 0x00, 0x7c, 0x00, 0xff, 0x01, 0xff, 0x01, 0xef, 0x01,
};
static unsigned char sdl_std_cursor_xterm[] = {
    0x00, 0x77, 0x00, 0x1c, 0x00, 0x08, 0x00, 0x08, 0x00, 0x08, 0x00, 0x08, 0x00, 0x08, 0x00, 0x08,
    0x00, 0x08, 0x00, 0x08, 0x00, 0x08, 0x00, 0x08, 0x00, 0x1c, 0x00, 0x77, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_xterm 9
#define sdl_std_cursor_HEIGHT_xterm 16
#define sdl_std_cursor_HOT_X_xterm -3
#define sdl_std_cursor_HOT_Y_xterm -7

static unsigned char sdl_std_cursor_mask_watch[] = {
    0xfc, 0x0f, 0xfc, 0x0f, 0xfc, 0x0f, 0xfe, 0x1f, 0xff, 0x3f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x3f, 0xfe, 0x1f, 0xfc, 0x0f, 0xfc, 0x0f, 0xfc, 0x0f,
};
static unsigned char sdl_std_cursor_watch[] = {
    0xf8, 0x07, 0xf8, 0x07, 0xf8, 0x07, 0xfc, 0x0f, 0x86, 0x18, 0x83, 0x30, 0x81, 0xe0, 0xc1, 0xe1,
    0xc1, 0xe1, 0x21, 0xe0, 0x13, 0x30, 0x06, 0x18, 0xfc, 0x0f, 0xf8, 0x07, 0xf8, 0x07, 0xf8, 0x07,
};
#define sdl_std_cursor_WIDTH_watch 16
#define sdl_std_cursor_HEIGHT_watch 16
#define sdl_std_cursor_HOT_X_watch -15
#define sdl_std_cursor_HOT_Y_watch -7

static unsigned char sdl_std_cursor_mask_umbrella[] = {
    0xe8, 0x76, 0xfb, 0xdf, 0xfd, 0x3f, 0xfe, 0xff, 0xff, 0x3f, 0xff, 0xff, 0xcf, 0x79, 0xc0, 0x01,
    0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x07, 0xc0, 0x07, 0xc0, 0x07, 0xc0, 0x07, 0x80, 0x03,
};
static unsigned char sdl_std_cursor_umbrella[] = {
    0x88, 0x04, 0x20, 0x0a, 0xc9, 0x32, 0xf2, 0x09, 0x4c, 0x06, 0x43, 0x18, 0x40, 0x00, 0x40, 0x00,
    0x40, 0x00, 0x40, 0x00, 0x40, 0x00, 0x40, 0x01, 0x40, 0x01, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_umbrella 16
#define sdl_std_cursor_HEIGHT_umbrella 16
#define sdl_std_cursor_HOT_X_umbrella -7
#define sdl_std_cursor_HOT_Y_umbrella -12

static unsigned char sdl_std_cursor_mask_top_side[] = {
    0xff, 0x7f, 0xff, 0x7f, 0xff, 0x7f, 0xff, 0x7f, 0xc0, 0x01, 0xe0, 0x03, 0xf0, 0x07, 0xf8, 0x0f,
    0xdc, 0x1d, 0xcc, 0x19, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01,
};
static unsigned char sdl_std_cursor_top_side[] = {
    0xff, 0x1f, 0xff, 0x1f, 0x00, 0x00, 0x40, 0x00, 0xe0, 0x00, 0x50, 0x01, 0x48, 0x02, 0x44, 0x04,
    0x40, 0x00, 0x40, 0x00, 0x40, 0x00, 0x40, 0x00, 0x40, 0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_top_side 15
#define sdl_std_cursor_HEIGHT_top_side 16
#define sdl_std_cursor_HOT_X_top_side -6
#define sdl_std_cursor_HOT_Y_top_side -14

static unsigned char sdl_std_cursor_mask_top_right_corner[] = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x00, 0xf0, 0xfc, 0xf7, 0xfc, 0xf7, 0xfc, 0xf7,
    0xc0, 0xf7, 0xe0, 0xf7, 0x70, 0xf7, 0x38, 0xf7, 0x1c, 0xf7, 0x0c, 0xf7, 0x00, 0xf0, 0x00, 0xf0,
};
static unsigned char sdl_std_cursor_top_right_corner[] = {
    0xff, 0x3f, 0xff, 0x3f, 0x00, 0x30, 0x00, 0x30, 0x00, 0x30, 0xfc, 0x31, 0x80, 0x31, 0x40, 0x31,
    0x20, 0x31, 0x10, 0x31, 0x08, 0x31, 0x04, 0x31, 0x00, 0x30, 0x00, 0x30, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_top_right_corner 16
#define sdl_std_cursor_HEIGHT_top_right_corner 16
#define sdl_std_cursor_HOT_X_top_right_corner -13
#define sdl_std_cursor_HOT_Y_top_right_corner -14

static unsigned char sdl_std_cursor_mask_top_left_corner[] = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x0f, 0x00, 0xef, 0x3f, 0xef, 0x3f, 0xef, 0x3f,
    0xef, 0x03, 0xef, 0x07, 0xef, 0x0e, 0xef, 0x1c, 0xef, 0x38, 0xef, 0x30, 0x0f, 0x00, 0x0f, 0x00,
};
static unsigned char sdl_std_cursor_top_left_corner[] = {
    0xff, 0x3f, 0xff, 0x3f, 0x03, 0x00, 0x03, 0x00, 0x03, 0x00, 0xe3, 0x0f, 0x63, 0x00, 0xa3, 0x00,
    0x23, 0x01, 0x23, 0x02, 0x23, 0x04, 0x23, 0x08, 0x03, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_top_left_corner 16
#define sdl_std_cursor_HEIGHT_top_left_corner 16
#define sdl_std_cursor_HOT_X_top_left_corner 0
#define sdl_std_cursor_HOT_Y_top_left_corner -14

static unsigned char sdl_std_cursor_mask_spraycan[] = {
    0x00, 0x0c, 0x18, 0x0d, 0x7c, 0x0d, 0x7c, 0x0d, 0x7e, 0x0d, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00,
    0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00, 0xff, 0x00,
};
static unsigned char sdl_std_cursor_spraycan[] = {
    0x00, 0x06, 0x80, 0x00, 0x2c, 0x06, 0x9e, 0x00, 0x16, 0x06, 0x3f, 0x00, 0x21, 0x00, 0x27, 0x00,
    0x25, 0x00, 0x27, 0x00, 0x25, 0x00, 0x27, 0x00, 0x27, 0x00, 0x21, 0x00, 0x21, 0x00, 0x3f, 0x00,
};
#define sdl_std_cursor_WIDTH_spraycan 12
#define sdl_std_cursor_HEIGHT_spraycan 16
#define sdl_std_cursor_HOT_X_spraycan -9
#define sdl_std_cursor_HOT_Y_spraycan -14

static unsigned char sdl_std_cursor_mask_sb_v_double_arrow[] = {
    0x38, 0x00, 0x7c, 0x00, 0xfe, 0x00, 0xff, 0x01, 0xff, 0x01, 0x7c, 0x00, 0x7c, 0x00, 0x7c,
    0x00, 0x7c, 0x00, 0x7c, 0x00, 0xff, 0x01, 0xff, 0x01, 0xfe, 0x00, 0x7c, 0x00, 0x38, 0x00,
};
static unsigned char sdl_std_cursor_sb_v_double_arrow[] = {
    0x10, 0x00, 0x38, 0x00, 0x7c, 0x00, 0xfe, 0x00, 0x28, 0x00, 0x28, 0x00, 0x28, 0x00, 0x28,
    0x00, 0x28, 0x00, 0x28, 0x00, 0x28, 0x00, 0xfe, 0x00, 0x7c, 0x00, 0x38, 0x00, 0x10, 0x00,
};
#define sdl_std_cursor_WIDTH_sb_v_double_arrow 9
#define sdl_std_cursor_HEIGHT_sb_v_double_arrow 15
#define sdl_std_cursor_HOT_X_sb_v_double_arrow -3
#define sdl_std_cursor_HOT_Y_sb_v_double_arrow -8

static unsigned char sdl_std_cursor_mask_sb_h_double_arrow[] = {
    0x18,
    0x0c,
    0x1c,
    0x1c,
    0xfe,
    0x3f,
    0xff,
    0x7f,
    0xff,
    0x7f,
    0xff,
    0x7f,
    0xfe,
    0x3f,
    0x1c,
    0x1c,
    0x18,
    0x0c,
};
static unsigned char sdl_std_cursor_sb_h_double_arrow[] = {
    0x00,
    0x00,
    0x08,
    0x08,
    0x0c,
    0x18,
    0xfe,
    0x3f,
    0x0f,
    0x78,
    0xfe,
    0x3f,
    0x0c,
    0x18,
    0x08,
    0x08,
    0x00,
    0x00,
};
#define sdl_std_cursor_WIDTH_sb_h_double_arrow 15
#define sdl_std_cursor_HEIGHT_sb_h_double_arrow 9
#define sdl_std_cursor_HOT_X_sb_h_double_arrow -7
#define sdl_std_cursor_HOT_Y_sb_h_double_arrow -4

static unsigned char sdl_std_cursor_mask_right_side[] = {
    0x00, 0xf0, 0x00, 0xf0, 0xc0, 0xf0, 0xc0, 0xf1, 0x80, 0xf3, 0x00, 0xf7, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0x00, 0xf7, 0x80, 0xf3, 0xc0, 0xf1, 0xc0, 0xf0, 0x00, 0xf0, 0x00, 0xf0,
};
static unsigned char sdl_std_cursor_right_side[] = {
    0x00, 0x30, 0x00, 0x30, 0x40, 0x30, 0x80, 0x30, 0x00, 0x31, 0x00, 0x32, 0xff, 0x37, 0x00,
    0x32, 0x00, 0x31, 0x80, 0x30, 0x40, 0x30, 0x00, 0x30, 0x00, 0x30, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_right_side 16
#define sdl_std_cursor_HEIGHT_right_side 15
#define sdl_std_cursor_HOT_X_right_side -13
#define sdl_std_cursor_HOT_Y_right_side -7

static unsigned char sdl_std_cursor_mask_right_ptr[] = {
    0x00, 0x03, 0x80, 0x03, 0xc0, 0x03, 0xe0, 0x03, 0xf0, 0x03, 0xf8, 0x03, 0xfc, 0x03, 0xfe, 0x03,
    0xff, 0x03, 0xff, 0x03, 0xf8, 0x03, 0xbc, 0x03, 0x3c, 0x03, 0x1e, 0x00, 0x1e, 0x00, 0x0c, 0x00,
};
static unsigned char sdl_std_cursor_right_ptr[] = {
    0x00, 0x80, 0x00, 0xc0, 0x00, 0xe0, 0x00, 0xf0, 0x00, 0xf8, 0x00, 0xfc, 0x00, 0xfe, 0x00, 0xff,
    0x00, 0xf8, 0x00, 0xd8, 0x00, 0x8c, 0x00, 0x0c, 0x00, 0x06, 0x00, 0x06, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_right_ptr 10
#define sdl_std_cursor_HEIGHT_right_ptr 16
#define sdl_std_cursor_HOT_X_right_ptr -7
#define sdl_std_cursor_HOT_Y_right_ptr -14

static unsigned char sdl_std_cursor_mask_question_arrow[] = {
    0xf8, 0x00, 0xfc, 0x01, 0xfe, 0x03, 0xff, 0x07, 0x8f, 0x07, 0x9f, 0x07, 0xde, 0x07, 0xfc, 0x03,
    0xf8, 0x01, 0xf8, 0x00, 0xf8, 0x00, 0xfc, 0x01, 0xfe, 0x03, 0xfc, 0x01, 0xf8, 0x00, 0x70, 0x00,
};
static unsigned char sdl_std_cursor_question_arrow[] = {
    0x7c, 0x00, 0xfe, 0x00, 0xc7, 0x01, 0x83, 0x01, 0x87, 0x01, 0xc6, 0x01, 0xe0, 0x00, 0x78, 0x00,
    0x38, 0x00, 0x28, 0x00, 0x28, 0x00, 0xee, 0x00, 0x6c, 0x00, 0x38, 0x00, 0x10, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_question_arrow 11
#define sdl_std_cursor_HEIGHT_question_arrow 16
#define sdl_std_cursor_HOT_X_question_arrow -4
#define sdl_std_cursor_HOT_Y_question_arrow -8

static unsigned char sdl_std_cursor_mask_pirate[] = {
    0xf0, 0x03, 0xf8, 0x07, 0xfc, 0x0f, 0xfe, 0x1f, 0xfe, 0x1f, 0xfc, 0x0f, 0xf8, 0x07, 0xf1, 0x83,
    0xf1, 0xe3, 0xf3, 0xf3, 0xef, 0x39, 0x1e, 0x1e, 0xe0, 0x01, 0xfe, 0xc7, 0xff, 0xff, 0x0f, 0x7c,
};
static unsigned char sdl_std_cursor_pirate[] = {
    0xe0, 0x01, 0xf0, 0x03, 0xf8, 0x07, 0xcc, 0x0c, 0xcc, 0x0c, 0xf8, 0x07, 0xf0, 0x03, 0xe0, 0x01,
    0xe1, 0x21, 0xe1, 0x61, 0xc2, 0x10, 0x1c, 0x0e, 0xe0, 0x01, 0xf8, 0x47, 0x0f, 0x7c, 0x01, 0x20,
};
#define sdl_std_cursor_WIDTH_pirate 16
#define sdl_std_cursor_HEIGHT_pirate 16
#define sdl_std_cursor_HOT_X_pirate -7
#define sdl_std_cursor_HOT_Y_pirate -4

static unsigned char sdl_std_cursor_mask_left_side[] = {
    0x0f, 0x00, 0x0f, 0x00, 0x0f, 0x03, 0x8f, 0x03, 0xcf, 0x01, 0xef, 0x00, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xef, 0x00, 0xcf, 0x01, 0x8f, 0x03, 0x0f, 0x03, 0x0f, 0x00, 0x0f, 0x00,
};
static unsigned char sdl_std_cursor_left_side[] = {
    0x03, 0x00, 0x03, 0x00, 0x83, 0x00, 0x43, 0x00, 0x23, 0x00, 0x13, 0x00, 0xfb, 0x3f, 0x13,
    0x00, 0x23, 0x00, 0x43, 0x00, 0x83, 0x00, 0x03, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_left_side 16
#define sdl_std_cursor_HEIGHT_left_side 15
#define sdl_std_cursor_HOT_X_left_side 0
#define sdl_std_cursor_HOT_Y_left_side -7

static unsigned char sdl_std_cursor_mask_left_ptr[] = {
    0x03, 0x00, 0x07, 0x00, 0x0f, 0x00, 0x1f, 0x00, 0x3f, 0x00, 0x7f, 0x00, 0xff, 0x00, 0xff, 0x01,
    0xff, 0x03, 0xff, 0x03, 0x7f, 0x00, 0xf7, 0x00, 0xf3, 0x00, 0xe0, 0x01, 0xe0, 0x01, 0xc0, 0x00,
};
static unsigned char sdl_std_cursor_left_ptr[] = {
    0x00, 0x00, 0x02, 0x00, 0x06, 0x00, 0x0e, 0x00, 0x1e, 0x00, 0x3e, 0x00, 0x7e, 0x00, 0xfe, 0x00,
    0xfe, 0x00, 0x3e, 0x00, 0x36, 0x00, 0x62, 0x00, 0x60, 0x00, 0xc0, 0x00, 0xc0, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_left_ptr 10
#define sdl_std_cursor_HEIGHT_left_ptr 16
#define sdl_std_cursor_HOT_X_left_ptr -8
#define sdl_std_cursor_HOT_Y_left_ptr -14

static unsigned char sdl_std_cursor_mask_exchange[] = {
    0xe3, 0x07, 0xf7, 0x0f, 0xff, 0x1f, 0xff, 0x3f, 0x3f, 0x38, 0xff, 0x30, 0xff, 0x00, 0xff, 0x00,
    0x00, 0xff, 0x00, 0xff, 0x0c, 0xfe, 0x1c, 0xfc, 0xfc, 0xff, 0xf8, 0xff, 0xf0, 0xef, 0xe0, 0xc7,
};
static unsigned char sdl_std_cursor_exchange[] = {
    0xf1, 0x03, 0xfb, 0x07, 0x1f, 0x0c, 0x09, 0x08, 0x19, 0x00, 0x3f, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x3f, 0x00, 0x26, 0x04, 0x24, 0x0c, 0x3e, 0xf8, 0x37, 0xf0, 0x23, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_exchange 16
#define sdl_std_cursor_HEIGHT_exchange 16
#define sdl_std_cursor_HOT_X_exchange -6
#define sdl_std_cursor_HOT_Y_exchange -8

static unsigned char sdl_std_cursor_mask_crosshair[] = {
    0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01,
};
static unsigned char sdl_std_cursor_crosshair[] = {
    0x80, 0x00, 0x80, 0x00, 0x80, 0x00, 0x80, 0x00, 0x80, 0x00, 0x80, 0x00, 0x80, 0x00, 0x7f, 0xff,
    0x80, 0x00, 0x80, 0x00, 0x80, 0x00, 0x80, 0x00, 0x80, 0x00, 0x80, 0x00, 0x80, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_crosshair 16
#define sdl_std_cursor_HEIGHT_crosshair 16
#define sdl_std_cursor_HOT_X_crosshair -7
#define sdl_std_cursor_HOT_Y_crosshair -8

static unsigned char sdl_std_cursor_mask_bottom_side[] = {
    0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xc0, 0x01, 0xcc, 0x19, 0xdc, 0x1d,
    0xf8, 0x0f, 0xf0, 0x07, 0xe0, 0x03, 0xc0, 0x01, 0xff, 0x7f, 0xff, 0x7f, 0xff, 0x7f, 0xff, 0x7f,
};
static unsigned char sdl_std_cursor_bottom_side[] = {
    0x40, 0x00, 0x40, 0x00, 0x40, 0x00, 0x40, 0x00, 0x40, 0x00, 0x40, 0x00, 0x44, 0x04, 0x48, 0x02,
    0x50, 0x01, 0xe0, 0x00, 0x40, 0x00, 0x00, 0x00, 0xff, 0x1f, 0xff, 0x1f, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_bottom_side 15
#define sdl_std_cursor_HEIGHT_bottom_side 16
#define sdl_std_cursor_HOT_X_bottom_side -6
#define sdl_std_cursor_HOT_Y_bottom_side -1

static unsigned char sdl_std_cursor_mask_bottom_right_corner[] = {
    0x00, 0xf0, 0x00, 0xf0, 0x0c, 0xf7, 0x1c, 0xf7, 0x38, 0xf7, 0x70, 0xf7, 0xe0, 0xf7, 0xc0, 0xf7,
    0xfc, 0xf7, 0xfc, 0xf7, 0xfc, 0xf7, 0x00, 0xf0, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
};
static unsigned char sdl_std_cursor_bottom_right_corner[] = {
    0x00, 0x30, 0x00, 0x30, 0x04, 0x31, 0x08, 0x31, 0x10, 0x31, 0x20, 0x31, 0x40, 0x31, 0x80, 0x31,
    0xfc, 0x31, 0x00, 0x30, 0x00, 0x30, 0x00, 0x30, 0xff, 0x3f, 0xff, 0x3f, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_bottom_right_corner 16
#define sdl_std_cursor_HEIGHT_bottom_right_corner 16
#define sdl_std_cursor_HOT_X_bottom_right_corner -13
#define sdl_std_cursor_HOT_Y_bottom_right_corner -1

static unsigned char sdl_std_cursor_mask_bottom_left_corner[] = {
    0x0f, 0x00, 0x0f, 0x00, 0xef, 0x30, 0xef, 0x38, 0xef, 0x1c, 0xef, 0x0e, 0xef, 0x07, 0xef, 0x03,
    0xef, 0x3f, 0xef, 0x3f, 0xef, 0x3f, 0x0f, 0x00, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
};
static unsigned char sdl_std_cursor_bottom_left_corner[] = {
    0x03, 0x00, 0x03, 0x00, 0x23, 0x08, 0x23, 0x04, 0x23, 0x02, 0x23, 0x01, 0xa3, 0x00, 0x63, 0x00,
    0xe3, 0x0f, 0x03, 0x00, 0x03, 0x00, 0x03, 0x00, 0xff, 0x3f, 0xff, 0x3f, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_bottom_left_corner 16
#define sdl_std_cursor_HEIGHT_bottom_left_corner 16
#define sdl_std_cursor_HOT_X_bottom_left_corner 0
#define sdl_std_cursor_HOT_Y_bottom_left_corner -1

static unsigned char sdl_std_cursor_mask_arrow[] = {
    0x00, 0xe0, 0x00, 0xf8, 0x00, 0xfe, 0x80, 0x7f, 0xe0, 0x7f, 0xf8, 0x3f, 0xfc, 0x3f, 0xfc, 0x1f,
    0xe0, 0x1f, 0xf0, 0x0f, 0xf8, 0x0f, 0x7c, 0x07, 0x3e, 0x07, 0x1f, 0x02, 0x0e, 0x00, 0x04, 0x00,
};
static unsigned char sdl_std_cursor_arrow[] = {
    0x00, 0x30, 0x00, 0x3c, 0x00, 0x1f, 0xc0, 0x1f, 0xf0, 0x0f, 0xfc, 0x0f, 0xc0, 0x07, 0xe0, 0x07,
    0x70, 0x03, 0x38, 0x03, 0x1c, 0x01, 0x0e, 0x01, 0x07, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00,
};
#define sdl_std_cursor_WIDTH_arrow 16
#define sdl_std_cursor_HEIGHT_arrow 16
#define sdl_std_cursor_HOT_X_arrow -13
#define sdl_std_cursor_HOT_Y_arrow -14
/* end cursor data */

static SDL_Cursor *sdl_std_cursor_array[(int)GHOST_kStandardCursorNumCursors] = {0};

/* utility function mostly a copy of SDL_CreateCursor but allows us to change
 * color and supports blenders flipped bits */
static SDL_Cursor *sdl_ghost_CreateCursor(
    const Uint8 *data, const Uint8 *mask, int w, int h, int hot_x, int hot_y)
{
  SDL_Surface *surface;
  SDL_Cursor *cursor;
  int x, y;
  Uint32 *pixel;
  Uint8 datab = 0, maskb = 0;
  const Uint32 black = 0xFF000000;
  const Uint32 white = 0xFFFFFFFF;
  const Uint32 transparent = 0x00000000;

  /* Make sure the width is a multiple of 8 */
  w = ((w + 7) & ~7);

  /* Create the surface from a bitmap */
  surface = SDL_CreateRGBSurface(0, w, h, 32, 0x00FF0000, 0x0000FF00, 0x000000FF, 0xFF000000);
  if (!surface) {
    return NULL;
  }
  for (y = 0; y < h; ++y) {
    pixel = (Uint32 *)((Uint8 *)surface->pixels + y * surface->pitch);
    for (x = 0; x < w; ++x) {
      if ((x % 8) == 0) {
        datab = *data++;
        maskb = *mask++;

        /* reverse bit order */
        datab = (datab * 0x0202020202ULL & 0x010884422010ULL) % 1023;
        maskb = (maskb * 0x0202020202ULL & 0x010884422010ULL) % 1023;
      }
      if (maskb & 0x80) {
        *pixel++ = (datab & 0x80) ? white : black;
      }
      else {
        *pixel++ = (datab & 0x80) ? white : transparent;
      }
      datab <<= 1;
      maskb <<= 1;
    }
  }

  cursor = SDL_CreateColorCursor(surface, hot_x, hot_y);

  SDL_FreeSurface(surface);

  return cursor;
}

/* TODO, this is currently never freed but it wont leak either. */
static void sdl_cursor_init(void)
{

#define DEF_CURSOR(name, ind) \
  { \
    sdl_std_cursor_array[(int)ind] = sdl_ghost_CreateCursor( \
        sdl_std_cursor_##name, \
        sdl_std_cursor_mask_##name, \
        sdl_std_cursor_WIDTH_##name, \
        sdl_std_cursor_HEIGHT_##name, \
        (sdl_std_cursor_WIDTH_##name + (sdl_std_cursor_HOT_X_##name)) - 1, \
        (sdl_std_cursor_HEIGHT_##name + (sdl_std_cursor_HOT_Y_##name)) - 1); \
    assert(sdl_std_cursor_array[(int)ind] != NULL); \
  } \
  (void)0

  DEF_CURSOR(left_ptr, GHOST_kStandardCursorDefault);
  DEF_CURSOR(right_ptr, GHOST_kStandardCursorRightArrow);
  DEF_CURSOR(left_ptr, GHOST_kStandardCursorLeftArrow);
  DEF_CURSOR(umbrella, GHOST_kStandardCursorInfo);  // TODO, replace this one.
  DEF_CURSOR(pirate, GHOST_kStandardCursorDestroy);
  DEF_CURSOR(question_arrow, GHOST_kStandardCursorHelp);
  DEF_CURSOR(exchange, GHOST_kStandardCursorCycle);
  DEF_CURSOR(spraycan, GHOST_kStandardCursorSpray);
  DEF_CURSOR(watch, GHOST_kStandardCursorWait);
  DEF_CURSOR(xterm, GHOST_kStandardCursorText);
  DEF_CURSOR(crosshair, GHOST_kStandardCursorCrosshair);
  DEF_CURSOR(sb_v_double_arrow, GHOST_kStandardCursorUpDown);
  DEF_CURSOR(sb_h_double_arrow, GHOST_kStandardCursorLeftRight);
  DEF_CURSOR(top_side, GHOST_kStandardCursorTopSide);
  DEF_CURSOR(bottom_side, GHOST_kStandardCursorBottomSide);
  DEF_CURSOR(left_side, GHOST_kStandardCursorLeftSide);
  DEF_CURSOR(right_side, GHOST_kStandardCursorRightSide);
  DEF_CURSOR(top_left_corner, GHOST_kStandardCursorTopLeftCorner);
  DEF_CURSOR(top_right_corner, GHOST_kStandardCursorTopRightCorner);
  DEF_CURSOR(bottom_right_corner, GHOST_kStandardCursorBottomRightCorner);
  DEF_CURSOR(bottom_left_corner, GHOST_kStandardCursorBottomLeftCorner);
  DEF_CURSOR(arrow, GHOST_kStandardCursorCopy);
  // DEF_CURSOR(arrow, GHOST_kStandardCursorCustom);
  DEF_CURSOR(arrow, GHOST_kStandardCursorPencil);

#undef DEF_CURSOR
}

GHOST_TSuccess GHOST_WindowSDL::setWindowCursorGrab(GHOST_TGrabCursorMode mode)
{
  return GHOST_kSuccess;
}

GHOST_TSuccess GHOST_WindowSDL::setWindowCursorShape(GHOST_TStandardCursor shape)
{
  if (sdl_std_cursor_array[0] == NULL) {
    sdl_cursor_init();
  }

  SDL_SetCursor(sdl_std_cursor_array[(int)shape]);
  return GHOST_kSuccess;
}

GHOST_TSuccess GHOST_WindowSDL::setWindowCustomCursorShape(GHOST_TUns8 bitmap[16][2],
                                                           GHOST_TUns8 mask[16][2],
                                                           int hotX,
                                                           int hotY)
{
  return setWindowCustomCursorShape(
      (GHOST_TUns8 *)bitmap, (GHOST_TUns8 *)mask, 16, 16, hotX, hotY, 0, 1);
}

GHOST_TSuccess GHOST_WindowSDL::setWindowCustomCursorShape(GHOST_TUns8 *bitmap,
                                                           GHOST_TUns8 *mask,
                                                           int sizex,
                                                           int sizey,
                                                           int hotX,
                                                           int hotY,
                                                           int fg_color,
                                                           int bg_color)
{
  if (m_sdl_custom_cursor) {
    SDL_FreeCursor(m_sdl_custom_cursor);
  }

  m_sdl_custom_cursor = sdl_ghost_CreateCursor(
      (const Uint8 *)bitmap, (const Uint8 *)mask, sizex, sizex, hotX, hotY);

  SDL_SetCursor(m_sdl_custom_cursor);
  return GHOST_kSuccess;
}

GHOST_TSuccess GHOST_WindowSDL::setWindowCursorVisibility(bool visible)
{
  SDL_ShowCursor(visible);
  return GHOST_kSuccess;
}

GHOST_TUns16 GHOST_WindowSDL::getDPIHint()
{
  int displayIndex = SDL_GetWindowDisplayIndex(m_sdl_win);
  if (displayIndex < 0) {
    return 96;
  }

  float ddpi;
  if (SDL_GetDisplayDPI(displayIndex, &ddpi, NULL, NULL) != 0) {
    return 96;
  }

  return (int)ddpi;
}
