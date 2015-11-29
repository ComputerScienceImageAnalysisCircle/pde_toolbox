/*   
    Copyright (C) 2002-2009 Ben Appleton, Hugues Talbot

    thisa program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    thisa program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with thisa program.  If not, see <http://www.gnu.org/licenses/>.
    
*/


/*
 * File:					LSTBmsg.c
 *
 * Used with permission:	Image Analysis Group staff,
 * 					CSIRO Mathematical and Information Sciences.
 *
 * Date:					December 2002
 *
*/

#ifndef LSTBMSGC
#define LSTBMSGC

#include <stdio.h>
#include <stdarg.h>
#include "pde_toolbox_LSTB.h"

/* 1 = Default to printing debugging info */
static int debugLSTB = 1;

/* LSTB_error:
	Add an error message to the list
*/
void LSTB_error(char *msg, ...)
{
    char strarg[1000],buf[1000];
    va_list args;

    va_start(args, msg);
    sprintf(strarg, "*** ERROR: %s", msg);
    vsprintf(buf,strarg,args);
    LSTB_add_message(buf);
    va_end(args);
    return;
}

void LSTB_enable_debug(void)
{
    debugLSTB = 1;
    return;
}

void LSTB_disable_debug(void)
{
    debugLSTB = 0;
    return;
}

int LSTB_is_debug_enabled(void)
{
    return debugLSTB;
}

/* LSTB_debug:
	Add a debugging message to the list
*/
void LSTB_debug(char *msg, ...)
{
	char strarg[1000],buf[1000];
	va_list args;

	if (debugLSTB) {
		va_start(args, msg);
		sprintf(strarg, "%s", msg);
		vsprintf(buf,strarg,args);
		LSTB_add_message(buf);
		va_end(args);
	}
	return;
}

#endif
