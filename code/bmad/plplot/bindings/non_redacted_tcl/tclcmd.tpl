/*--------------------------------------------------------------------------*\
 * %cmd%Cmd
 *
 * Processes %cmd% Tcl command.
\*--------------------------------------------------------------------------*/

static int
%cmd%Cmd( ClientData PL_UNUSED( clientData ), Tcl_Interp *interp, int argc, const char *argv[] )
{
<argdecls>

    pl_errcode = 0; errmsg[0] = '\0';

    if ( (argc == 2) && (strlen(argv[1])>0) && (strncmp(argv[1],"-help",strlen(argv[1])) == 0) ) {
	Tcl_AppendResult( interp, "command syntax: \"",
			  "%cmd% %args%", "\"",
			  (char *) NULL);
	return TCL_ERROR;
    }

    if ( (!%isref% && %ndefs% && (argc < (1 + %nargs% - %ndefs%))) ||
         (!%isref% && !%ndefs% && (argc != (%nargs% + 1))) ||
         ( %isref% && (argc != 1) && (argc != (%nargs% + 1))) ) {
	Tcl_AppendResult( interp, "wrong # args: should be \"",
			  "%cmd% %args%", "\"",
			  (char *) NULL);
	return TCL_ERROR;
    }

<getargs>

<plcmd>

<fetch_result>

    if (pl_errcode != 0) {
	Tcl_AppendResult(interp, errmsg, (char *) NULL);
	return TCL_ERROR;
    }

    plflush();
    return TCL_OK;
}

