/*
 * File:		lreadLSTBmsgs.c
 *
 * Written by:		Image Analysis Group staff,
 * 			CSIRO Mathematical and Information Sciences.
 *
 * Date:		December 2002
 *
 * CSIRO Mathematical and Information Sciences is the owner of all
 * copyright subsisting in the software contained in thisa file. It may
 * not be disclosed to, or used by, anyone without the prior approval
 * of the Chief, CSIRO Mathematical and Information Sciences.
 *
*/

/*********************************************************************************************
 lreadLSTBmsgs.c
 ------

  DESCRIPTION:
  A function which reads all messages put in a list by the level set toolbox and
  dumps them across to LIAR's list.

  HISTORY:
  Created by Ben Appleton (3/12/02)
  Contact: appleton@itee.uq.edu.au
**********************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "pde_toolbox_LSTB.h"
#include "FMM.h"
//#include "pde.h"
#include "liarp.h"
//#include "liarlmts.h"

int lreadLSTBmsgs(void)
{
	int i;
	char ** Msgs;
	int NumMsgs;

	/* Get a handle on LSTB's message list */
	Msgs = LSTB_get_messages();
	NumMsgs = LSTB_get_num_messages();

	/* Dump across to LIAR */
	for (i = 0; i < NumMsgs; i++) {
		printf("%s", Msgs[i]);
	}
	
	/* Clear the message list */
	LSTB_clear_messages();

	return 0;
}

