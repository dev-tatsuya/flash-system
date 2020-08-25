# coding: UTF-8

from pygame.locals import *

_key_name = {"backspace": K_BACKSPACE, "tab": K_TAB, "clear": K_CLEAR, "return": K_RETURN, "pause": K_PAUSE, "escape": K_ESCAPE,
     "space": K_SPACE, "exclaim": K_EXCLAIM, "quotedbl": K_QUOTEDBL, "hash": K_HASH, "dollar": K_DOLLAR,
     "ampersand": K_AMPERSAND, "quote": K_QUOTE, "left parenthesis": K_LEFTPAREN, "right parenthesis": K_RIGHTPAREN,
     "asterisk": K_ASTERISK, "plus sign": K_PLUS, "comma": K_COMMA, "minus sign": K_MINUS, "period": K_PERIOD,
     "forward slash": K_SLASH, "0": K_0, "1": K_1, "2": K_2, "3": K_3, "4": K_4, "5": K_5, "6": K_6, "7": K_7, "8": K_8,
     "9": K_9, "colon": K_COLON, "semicolon": K_SEMICOLON, "less-than sign": K_LESS, "equals sign": K_EQUALS,
     "greater-than sign": K_GREATER, "question mark": K_QUESTION, "at": K_AT, "left bracket": K_LEFTBRACKET,
     "backslash": K_BACKSLASH, "right bracket": K_RIGHTBRACKET, "caret": K_CARET, "underscore": K_UNDERSCORE,
     "grave": K_BACKQUOTE, "a": K_a, "b": K_b, "c": K_c, "d": K_d, "e": K_e, "f": K_f, "g": K_g, "h": K_h, "i": K_i,
     "j": K_j, "k": K_k, "l": K_l, "m": K_m, "n": K_n, "o": K_o, "p": K_p, "q": K_q, "r": K_r, "s": K_s, "t": K_t,
     "u": K_u, "v": K_v, "w": K_w, "x": K_x, "y": K_y, "z": K_z, "delete": K_DELETE, "keypad 0": K_KP0,
     "keypad 1": K_KP1, "keypad 2": K_KP2, "keypad 3": K_KP3, "keypad 4": K_KP4, "keypad 5": K_KP5, "keypad 6": K_KP6,
     "keypad 7": K_KP7, "keypad 8": K_KP8, "keypad 9": K_KP9, "keypad period": K_KP_PERIOD,
     "keypad divide": K_KP_DIVIDE, "keypad multiply": K_KP_MULTIPLY, "keypad minus": K_KP_MINUS,
     "keypad plus": K_KP_PLUS, "keypad enter": K_KP_ENTER, "keypad equals": K_KP_EQUALS, "up arrow": K_UP,
     "down arrow": K_DOWN, "right arrow": K_RIGHT, "left arrow": K_LEFT, "insert": K_INSERT, "home": K_HOME,
     "end": K_END, "page up": K_PAGEUP, "page down": K_PAGEDOWN, "F1": K_F1, "F2": K_F2, "F3": K_F3, "F4": K_F4,
     "F5": K_F5, "F6": K_F6, "F7": K_F7, "F8": K_F8, "F9": K_F9, "F10": K_F10, "F11": K_F11, "F12": K_F12, "F13": K_F13,
     "F14": K_F14, "F15": K_F15, "numlock": K_NUMLOCK, "capslock": K_CAPSLOCK, "scrollock": K_SCROLLOCK,
     "right shift": K_RSHIFT, "left shift": K_LSHIFT, "right control": K_RCTRL, "left control": K_LCTRL,
     "right alt": K_RALT, "left alt": K_LALT, "right meta": K_RMETA, "left meta": K_LMETA, "left Windows key": K_LSUPER,
     "right Windows key": K_RSUPER, "mode shift": K_MODE, "help": K_HELP, "print screen": K_PRINT, "sysrq": K_SYSREQ,
     "break": K_BREAK, "menu": K_MENU, "power": K_POWER, "Euro": K_EURO}

def get_key(key):
     if type(key) == int:
          return key
     elif type(key) == str:
          return _key_name[key]
     else:
          raise Exception("key must be int or string.")