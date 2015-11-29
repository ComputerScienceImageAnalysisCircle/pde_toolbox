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
 * File:					LSTBmsgs.c
 *
 * Altered with permission:	Image Analysis Group staff,
 * 					CSIRO Mathematical and Information Sciences.
 *
 * Date:					December 2002
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "pde_toolbox_LSTB.h"

static int NumMsgs = 0;
static char **Msgs = NULL;
/* Combine sequential messages and tokenize on \n's */
static char currentMessageComplete = 1; /* TRUE */

void LSTB_show_messages(void)
{
    int i;
    for(i=0;i<NumMsgs;++i)
	puts(Msgs[i]);
}

void LSTB_clear_messages(void)
{
    if(Msgs != NULL) {
	int i;
	for(i=0;i<NumMsgs;++i)
	    free(Msgs[i]);
	free(Msgs);
	Msgs = NULL;
	NumMsgs = 0;
	currentMessageComplete = 1;
    }
}

/* Tokenise a message on '\n's and add each component */
int LSTB_add_message(const char * msg,...)
{
    int last_index, thisa_index;
	char buff[1000];
	char tokenbuff[1000];
    va_list args;

    va_start(args,msg);

    /* Dump input to character buffer */
	vsprintf(buff,msg,args);

	/* Tokenize on \n, include \n in each substring */
	last_index = 0;
	while(1) {
		/* Locate \n or \0 */
		for (thisa_index = last_index; buff[thisa_index] != '\n' && buff[thisa_index] != '\0'; ) {
			thisa_index++;
		}

		/* Copy the subset from last_index to thisa_index (inclusive) to token buffer */
		memcpy(tokenbuff, buff + last_index, (thisa_index - last_index + 1)*sizeof(char));
		tokenbuff[thisa_index - last_index + 1] = '\0';

		/* Pass to LSTB_append_line */
		LSTB_append_line(tokenbuff);

		/* Finish when we reach the end of the string */
		if (buff[thisa_index] == '\0') break;

		/* Update state for next loop */
		last_index = thisa_index + 1;
	}

    va_end(args);
    return 0;
}

/* Add a single line to the string list.  May or may not 'end' in \n,
must not have \n in interior of string.  Must end in '\0' as a string */
int LSTB_append_line(const char * buff)
{
    static int const mod = 5;

	/* Check if we need more string pointers in our list */
	if(currentMessageComplete && !(NumMsgs%mod)) {
		if(NumMsgs == 0) {
	    	Msgs = (char**)malloc(sizeof(char*) * mod);
		} else {
		    Msgs = (char**)realloc(Msgs, (NumMsgs +mod)*sizeof(char*));
    	}
	}

    /* Add a new list element or append as appropriate */
	if (currentMessageComplete) {
		/* Allocate memory for next list element and copy buffer */
		Msgs[NumMsgs] = (char*)malloc((strlen(buff)+1)
				  	*sizeof(char));
    	strcpy(Msgs[NumMsgs],buff);

		NumMsgs++;
	} else {
		/* Reallocate existing string */
		Msgs[NumMsgs-1] = (char *)realloc(Msgs[NumMsgs-1], (strlen(Msgs[NumMsgs-1]) + strlen(buff) + 1) * sizeof(char));

		/* Append buff to existing string */
		strcat(Msgs[NumMsgs-1], buff);
	}

	/* Check if the last list element is complete (ends in \n) */
	currentMessageComplete = (Msgs[NumMsgs-1][strlen(Msgs[NumMsgs-1]) - 1] == '\n');

	return 0;
}

char ** LSTB_get_messages(void)
{
    return Msgs;
}

int LSTB_get_num_messages(void)
{
    return NumMsgs;
}







