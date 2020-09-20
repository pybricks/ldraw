//---------------------------------------------------------------------------

#pragma hdrstop

/*
 * Program Name: LSynth
 * Program Description:
 *   LDraw CAD system compatible flexible part synthesizer.  This program
 *   reads your LDraw file, searching for unofficial META commands (structured
 *   comments), that specify things you want synthesized.
 *
 *   LSynth has two primary forms of synthesis.  The first kind of synthesis
 *   creates hose like things including:
 *     rubbed, rubber, pneumatic, and flex system hoses, as well as electric
 *     and fiber optic cables, flex system cables, flexible axles, string, and
 *     minifig chain.
 *
 *   It also creates things that travel around circular lego parts.  Things
 *   like:
 *     rubber band, rubber belt, technic chain, technic plastic tread, technic
 *     rubber tread.
 *
 *   The files tube.c, tube.h, curve.c and curve.h perform hose synthesis.
 *   The files band.c and band.h perform band synthesis.
 *
 *   This file (main.c) contains the main entry/exit points for the program.
 *   It opens and scans the LDraw file provided, identifies synthesis
 *   synthesis specifications and hands them off to the appropriate synthesis
 *   methodology.
 */

#pragma hdrstop

#ifdef WIN32
// Added a messagebox for lsynth.mpd version mismatch.
#include <windows.h>
#endif

#include <ctype.h>

#include "lsynthcp.h"
#include "hose.h"
#include "band.h"

typedef struct {
  char name[126];
  char nickname[128];
  char method[128];
} product_t;

product_t products[256];
int       n_products = 0;

char version[] = "3.1";
char beta[] = ""; // " Beta I";

char mpdversion[32] = "UNKNOWN";

int ldraw_part = 0;
int group_size;

//---------------------------------------------------------------------------
void messagebox( const char* title, const char* message )
{
  char cmd[1024];
  int i;

#if defined( __APPLE__ ) && defined( __MACH__ )
  \\ Skip activate if we dont want to focus on dialog.  Finder will bounce instead.
  sprintf(cmd, "osascript -e \'tell app \"Finder\" to activate display dialog \"%s\"\'", message);
  i = system(cmd);

#elif WIN32
  MessageBox( NULL, message, title, MB_OK );

#else
  // Unknown OS.  We'll probably have to settle for commandline warnings. 
  // But first attempt to launch an oldsKool stylee xmessage.  
  // Gnome systems often come with the much prettier zenity.  Try it first.
  sprintf(cmd, "zenity --warning --title=\"%s\" --text=\"%s\" 2>1 >/dev/null", title, message);
  i = system(cmd);
  i = i >> 8; // Get exit value of zenity.
  if ((i == 0) || (i == 1)) // (0=OK, 1=X)  
    return;

  // No gnome?  Try kdialog from KDE.
  sprintf(cmd, "kdialog --title=\"%s\" --sorry \"%s\" 2>1 >/dev/null", title, message);
  i = system(cmd);
  i = i >> 8; // Get exit value of kdialog.
  if ((i == 0) || (i == 1)) // (0=OK, 1=X)  
    return;

  // The unpatched Athena widgets are hideous, but there's a pretty good chance it'll work.
  sprintf(cmd, "xmessage -bg lightgrey -fn 9x15bold -buttons OK -center -title \"%s\" \"%s\"", title, message);
  // I should probably skip the fork and just use system to make it a modal messagebox.
  if(fork()==0){
    close(1); close(2);
    system(cmd);
    exit(0);
  }
  
  // Still here?  Sorry, yer stuck with command line output.
#endif
}

//---------------------------------------------------------------------------
   /* If this code works, it was written by Lars C. Hassing. */
   /* If not, I don't know who wrote it.                     */

   /* Like fgets, except that 1) any line ending is accepted (\n (unix),
   \r\n (DOS/Windows), \r (Mac (OS9)) and 2) Str is ALWAYS zero terminated 
   (even if no line ending was found) */
   //--------------------------------------------------------------------------- 
   char *L3fgets(char *Str, int n, FILE *fp) 
   { 
      register int   c; 
      int            nextc; 
      register char *s = Str; 
    
      while (--n > 0) 
      {
         if ((c = getc(fp)) == EOF)
            break;
         if (c == '\032')
            continue;              /* Skip CTRL+Z                               */
         if (c == '\r' || c == '\n')
         {
            *s++ = '\n';
            /* We got CR or LF, eat next character if LF or CR respectively */
            if ((nextc = getc(fp)) == EOF)
               break;
            if (nextc == c || (nextc != '\r' && nextc != '\n'))
               ungetc(nextc, fp);  /* CR-CR or LF-LF or ordinary character      */
            break;
         }
         *s++ = c;
      }
      *s = 0;

      /* if (ferror(fp)) return NULL; if (s == Str) return NULL; */
      if (s == Str)
         return NULL;

      return Str;
   }

   //---------------------------------------------------------------------------

char *
fgetline(
  char *line,
  int   len,
  FILE *file)
{
  char *rc;
  while ((rc = L3fgets(line,len,file))) {
    char *nonwhite;

    nonwhite = line + strspn(line," \t");

    if (strncasecmp(nonwhite,"0 ROTATION C",strlen("0 ROTATION C")) == 0 ||
        strncasecmp(nonwhite,"0 COLOR",strlen("0 COLOR")) == 0) {
      continue;
    }

    nonwhite = line + strspn(line," \t");
    if (strncasecmp(nonwhite,"0 WRITE ",strlen("0 WRITE ")) == 0) {
      strcpy(nonwhite + 2, nonwhite + strlen("0 WRITE "));
    }
    break;
  }
  return rc;
}

void
strclean(char *str)
{
  if (strncasecmp(str,"0 WRITE ",strlen("0 WRITE ")) == 0) {
    strcpy(str + 2, str + strlen("0 WRITE "));
  }
}

//---------------------------------------------------------------------------

/*
 * Skip over MLCad ROTATION and COLOR statements
 */

static int skip_rot(char *line, int n_line, FILE *dat, FILE *temp)
{
  char *nonwhite;

  nonwhite = line + strspn(line,"\t");

  while (strncasecmp(nonwhite,"0 ROTATION C",strlen("0 ROTATION C")) == 0 ||
         strncasecmp(nonwhite,"0 COLOR",strlen("0 COLOR")) == 0) {
    fputs(line,temp);
    L3fgets(line,n_line,dat);  /* FIXME: check fgets rc */
    nonwhite = line + strspn(line," \t");
  }

  return 0;
}


/*****************************************************************************
 *
 * Read in and parse up synthesis descriptions and constraints.
 *
 ****************************************************************************/

int
parse_descr(char *fullpath_progname)
{
  char filename[256];
  FILE *mpd;
  char line[256];

  strcpy(filename,fullpath_progname);
  {
    char *l, *p;

    for (l = p = filename; *p; *p++) {
      if (*p == '\\' || *p == '/') {
        l = p+1;
      }
    }
    *l = '\0';
  }
  strcat(filename,"lsynth.mpd");

  mpd = fopen(filename,"r");

  if (mpd == NULL) {
    printf("Failed to open lsynth.mpd for reading.\n");
    messagebox("LSynth", "Failed to open lsynth.mpd for reading.");
    return -1;
  }

  while(fgetline(line,sizeof(line),mpd)) {
    char stretch[64];
    char type[64];
    char product[126], nickname[128], method[128];
    int  d,st,i;
    PRECISION s,t;
    int got_end = 0;

    strclean(line);

    if (sscanf(line,"0 !VERSION %d.%d\n", &i, &d) == 2) {
      sprintf(mpdversion, "%d.%d", i, d);
    }
    else if (!strncmp(line,"0 SYNTH ", 8))
    {
      if (strcmp(mpdversion, version))
      {
        char s[256];
        sprintf(s, "\nWarning: lsynth.mpd version %s does not match executable version %s!", 
                mpdversion, version);
        printf("%s\n\n", s);
        messagebox("LSynth", s);
        strcpy(mpdversion, version); // Only warn once.
      }
    }

    if (sscanf(line,"0 SYNTH PART %s %s %s\n",product, nickname, method) == 3) {
      strcpy(products[n_products].name,product);
      strcpy(products[n_products].nickname,nickname);
      strcpy(products[n_products].method,method);
      n_products++;
    } else if (sscanf(line,"0 SYNTH BEGIN DEFINE %s HOSE %s %d %d %f\n",
        type,stretch,&d,&st,&t) == 5) {
      if (strcasecmp(stretch,"STRETCH") == 0) {
        hose_types[n_hose_types].fill = STRETCH;
      } else if (strcasecmp(stretch,"FIXED") == 0) {
        hose_types[n_hose_types].fill = FIXED;
      } else if ((strncasecmp(stretch,"FIXED",strlen("FIXED")) == 0) &&
                 (sscanf(stretch, "FIXED%d", &i) == 1) && (i >1)) {
        hose_types[n_hose_types].fill = i;
      } else {
        printf("Error: Unrecognized fill type %s for hose type %s.  Aborting\n",
          stretch,type);
        fclose(mpd);
        return -1;
      }

      strcpy(hose_types[n_hose_types].type,type);
      hose_types[n_hose_types].diameter = d;
      hose_types[n_hose_types].stiffness = st;
      hose_types[n_hose_types].twist = t;

      for (i = 0; i < 3; i++) {
        part_t *part;
        int     got_part = 0;

        if (i == 0) {
          part = &hose_types[n_hose_types].start;
        } else if (i == 1) {
          part = &hose_types[n_hose_types].mid;
        } else {
          part = &hose_types[n_hose_types].end;
        }

        while (fgetline(line,sizeof(line),mpd)) {

          int n;

          n = sscanf(line,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
                &part->attrib,
                &part->offset[0],    &part->offset[1],    &part->offset[2],
                &part->orient[0][0], &part->orient[0][1], &part->orient[0][2],
                &part->orient[1][0], &part->orient[1][1], &part->orient[1][2],
                &part->orient[2][0], &part->orient[2][1], &part->orient[2][2],
                 part->type);

          if (n == 14) {
            got_part = 1;
            break;
          }
        }
        if ( ! got_part) {
          printf("Error: Unexpected end of file\n");
          fclose(mpd);
          return -1;
        }
      }

      // Assume no alternate mid part.
      strcpy(hose_types[n_hose_types].alt.type, "");

      got_end = 0;
      while (fgetline(line,sizeof(line),mpd)) {
	part_t *part;
	int n;

        if (strcasecmp(line,"0 SYNTH END\n") == 0) {
	  got_end = 1;
	  break;
	}

	// Look for an alternate mid part
	part = &hose_types[n_hose_types].alt;
	
	n = sscanf(line,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
                &part->attrib,
                &part->offset[0],    &part->offset[1],    &part->offset[2],
                &part->orient[0][0], &part->orient[0][1], &part->orient[0][2],
                &part->orient[1][0], &part->orient[1][1], &part->orient[1][2],
                &part->orient[2][0], &part->orient[2][1], &part->orient[2][2],
                 part->type);

	if (n != 14)
	  strcpy(hose_types[n_hose_types].alt.type, ""); // Skip comments
	else
        {
#ifdef DEBUGGING_HOSES
	  printf("Found HOSE alt segment %s\n", hose_types[n_hose_types].alt.type);
#endif
        }
      }
      if ( ! got_end) {
	printf("Error: Unexepcted end of file\n");
	fclose(mpd);
	return -1;
      }
      n_hose_types++;
    } else if (strcasecmp(line,"0 SYNTH BEGIN DEFINE HOSE CONSTRAINTS\n") == 0) {
      while(fgetline(line,sizeof(line),mpd)) {
        part_t *part = &hose_constraints[n_hose_constraints];
        if (sscanf(line,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          &part->attrib,
          &part->offset[0],    &part->offset[1],    &part->offset[2],
          &part->orient[0][0], &part->orient[0][1], &part->orient[0][2],
          &part->orient[1][0], &part->orient[1][1], &part->orient[1][2],
          &part->orient[2][0], &part->orient[2][1], &part->orient[2][2],
           part->type) == 14) {
          n_hose_constraints++;
        } else if (strcasecmp(line,"0 SYNTH END\n") == 0) {
          break;
        }
      }
    } else if (sscanf(line,"0 SYNTH BEGIN DEFINE %s BAND %s %f %f\n",
                 type,stretch,&s,&t) == 4) {
      int n;

      if (strcasecmp(stretch,"STRETCH") == 0) {
        band_types[n_band_types].fill = STRETCH;
        n = 2;
      } else if (strcasecmp(stretch,"FIXED") == 0) {
        band_types[n_band_types].fill = FIXED;
        n = 2;
      } else if (strcasecmp(stretch,"FIXED3") == 0) {
        band_types[n_band_types].fill = FIXED3;
        n = 4;
      } else {
        printf("Error: Unrecognized fill type %s for hose type %s.  Aborting\n",
          stretch,type);
        fclose(mpd);
        return -1;
      }

      strcpy(band_types[n_band_types].type,type);
      band_types[n_band_types].scale = s;
      band_types[n_band_types].thresh = t;

      for (i = 0; i < n; i++) {
        part_t *part;
        int     got_part = 0;

        if (i == 0) {
          part = &band_types[n_band_types].tangent;
        } else if (i == 1) {
          part = &band_types[n_band_types].arc;
        } else if (i == 2) {
          part = &band_types[n_band_types].start_trans;
        } else {
          part = &band_types[n_band_types].end_trans;
        }

        while (fgetline(line,sizeof(line),mpd)) {
          if (sscanf(line,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
            &part->attrib,
            &part->offset[0],    &part->offset[1],    &part->offset[2],
            &part->orient[0][0], &part->orient[0][1], &part->orient[0][2],
            &part->orient[1][0], &part->orient[1][1], &part->orient[1][2],
            &part->orient[2][0], &part->orient[2][1], &part->orient[2][2],
             part->type) == 14) {
            got_part = 1;
            break;
          }
        }
        if ( ! got_part) {
          printf("Error: Unexpected end of file\n");
          fclose(mpd);
          return -1;
        }
      }

      if (L3fgets(line,sizeof(line),mpd)) {
        if (strcasecmp(line,"0 SYNTH END\n") != 0) {
          printf("Error: Expected SYNTH END, got this instead\n");
          printf(line);
          fclose(mpd);
          return -1;
        }
      } else {
        printf("Error: Unexepcted end of file\n");
        fclose(mpd);
        return -1;
      }
      n_band_types++;
    } else if (sscanf(line,"0 SYNTH BEGIN DEFINE %s BAND %s %f %f\n",
                 type,stretch,&s,&t) == 4) {
      int n;

      if (strcasecmp(stretch,"STRETCH") == 0) {
        band_types[n_band_types].fill = STRETCH;
        n = 2;
      } else if (strcasecmp(stretch,"FIXED") == 0) {
        band_types[n_band_types].fill = FIXED;
        n = 2;
      } else if (strcasecmp(stretch,"FIXED3") == 0) {
        band_types[n_band_types].fill = FIXED3;
        n = 4;
      } else {
        printf("Error: Unrecognized fill type %s for hose type %s.  Aborting\n",
          stretch,type);
        fclose(mpd);
        return -1;
      }

      strcpy(band_types[n_band_types].type,type);
      band_types[n_band_types].scale = s;
      band_types[n_band_types].thresh = t;
      band_types[n_band_types].pulley = 0;

      for (i = 0; i < n; i++) {
        part_t *part;
        int     got_part = 0;

        if (i == 0) {
          part = &band_types[n_band_types].tangent;
        } else if (i == 1) {
          part = &band_types[n_band_types].arc;
        } else if (i == 2) {
          part = &band_types[n_band_types].start_trans;
        } else {
          part = &band_types[n_band_types].end_trans;
        }

        while(fgetline(line,sizeof(line),mpd)) {
          if (sscanf(line,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
            &part->attrib,
            &part->offset[0],    &part->offset[1],    &part->offset[2],
            &part->orient[0][0], &part->orient[0][1], &part->orient[0][2],
            &part->orient[1][0], &part->orient[1][1], &part->orient[1][2],
            &part->orient[2][0], &part->orient[2][1], &part->orient[2][2],
             part->type) == 14) {
            got_part = 1;
            break;
          }
        }
        if ( ! got_part) {
          printf("Error: Unexpected end of file\n");
          fclose(mpd);
          return -1;
        }
      }

      if (L3fgets(line,sizeof(line),mpd)) {
        if (strcasecmp(line,"0 SYNTH END\n") != 0) {
          printf("Error: Expected SYNTH END, got this instead\n");
          printf(line);
          fclose(mpd);
          return -1;
        }
      } else {
        printf("Error: Unexepcted end of file\n");
        fclose(mpd);
        return -1;
      }
      n_band_types++;
    } else if (sscanf(line,"0 SYNTH BEGIN DEFINE %s PULLEY %s %f %f\n",
                 type,stretch,&s,&t) == 4) {
      int n;

      if (strcasecmp(stretch,"STRETCH") == 0) {
        band_types[n_band_types].fill = STRETCH;
        n = 2;
      } else if (strcasecmp(stretch,"FIXED") == 0) {
        band_types[n_band_types].fill = FIXED;
        n = 2;
      } else if (strcasecmp(stretch,"FIXED3") == 0) {
        band_types[n_band_types].fill = FIXED3;
        n = 4;
      } else {
        printf("Error: Unrecognized fill type %s for hose type %s.  Aborting\n",
          stretch,type);
        fclose(mpd);
        return -1;
      }

      strcpy(band_types[n_band_types].type,type);
      band_types[n_band_types].scale = s;
      band_types[n_band_types].thresh = t;

      band_types[n_band_types].pulley = 1;

      for (i = 0; i < n; i++) {
        part_t *part;
        int     got_part = 0;

        if (i == 0) {
          part = &band_types[n_band_types].tangent;
        } else if (i == 1) {
          part = &band_types[n_band_types].arc;
        } else if (i == 2) {
          part = &band_types[n_band_types].start_trans;
        } else {
          part = &band_types[n_band_types].end_trans;
        }

        while (fgetline(line,sizeof(line),mpd)) {
          if (sscanf(line,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
            &part->attrib,
            &part->offset[0],    &part->offset[1],    &part->offset[2],
            &part->orient[0][0], &part->orient[0][1], &part->orient[0][2],
            &part->orient[1][0], &part->orient[1][1], &part->orient[1][2],
            &part->orient[2][0], &part->orient[2][1], &part->orient[2][2],
             part->type) != 14) {
            got_part = 1;
            break;
          }
        }
        if ( ! got_part) {
          printf("Error: Unexpected end of file\n");
          fclose(mpd);
          return -1;
        }
      }

      if (L3fgets(line,sizeof(line),mpd)) {
        if (strcasecmp(line,"0 SYNTH END\n") != 0) {
          printf("Error: Expected SYNTH END, got this instead\n");
          printf(line);
          fclose(mpd);
          return -1;
        }
      } else {
        printf("Error: Unexepcted end of file\n");
        fclose(mpd);
        return -1;
      }
      n_band_types++;

    } else if (strcasecmp(line,"0 SYNTH BEGIN DEFINE BAND CONSTRAINTS\n") == 0) {
      while(fgetline(line,sizeof(line),mpd)) {
        part_t *part = &band_constraints[n_band_constraints];
        if (sscanf(line,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s\n",
          &part->attrib,
          &part->offset[0],    &part->offset[1],    &part->offset[2],
          &part->orient[0][0], &part->orient[0][1], &part->orient[0][2],
          &part->orient[1][0], &part->orient[1][1], &part->orient[1][2],
          &part->orient[2][0], &part->orient[2][1], &part->orient[2][2],
           part->type) == 14) {
          n_band_constraints++;
        } else if (strcasecmp(line,"0 SYNTH END\n") == 0) {
          break;
        }
      }
    }
  }

  fclose(mpd);
  return 0;
}

void list_products()
{
  int i;

  printf("\n\nComplete parts LSynth can create\n");
  for (i = 0; i < n_products; i++) {
    printf("  %s %s (%s)\n",
      products[i].name,
      products[i].nickname,
      products[i].method);
  }
}

void product_ini(void)
{
  int i;

  for (i = 0; i < n_products; i++) {
    printf("%-20s = SYNTH BEGIN %s 16\n",
      products[i].nickname,
      products[i].nickname);
  }

  for (i = 0; i < n_products; i++) {
    printf("%-20s = SYNTH BEGIN %s 16\n",
      products[i].name,
      products[i].name);
  }
}

char *
isproduct(char *type)
{
  int i;

  for (i = 0; i < n_products; i++) {
    if (strncasecmp(products[i].name,type,strlen(products[i].name)) == 0) {
      return products[i].name;
    }
    if (strncasecmp(products[i].nickname,type,strlen(products[i].nickname)) == 0) {
      return products[i].name;
    }
  }
  return NULL;
}

char *
product_method(char *type)
{
  int i;

  for (i = 0; i < n_products; i++) {
    if (strncasecmp(products[i].name,type,strlen(products[i].name)) == 0) {
      return products[i].method;
    }
    if (strncasecmp(products[i].nickname,type,strlen(products[i].nickname)) == 0) {
      return products[i].method;
    }
  }
  return NULL;
}

char *
product_nickname(char *type)
{
  int i;

  for (i = 0; i < n_products; i++) {
    if (strncasecmp(products[i].name,type,strlen(products[i].name)) == 0) {
      return products[i].nickname;
    }
    if (strncasecmp(products[i].nickname,type,strlen(products[i].nickname)) == 0) {
      return products[i].nickname;
    }
  }
  return NULL;
}


/*
 * Skip over results rom previous syntesis efforts.
 */

int skip_synthesized(FILE *dat, char *line, int sizeof_line)
{
  int rc;
  char *nonwhite;

  while (L3fgets(line,sizeof(line),dat)) {
    nonwhite = line + strspn(line," \t");
    if (strcasecmp(nonwhite,"0 SYNTHESIZED END\n") == 0) {
      return 0;
    }
  }
  return -1;
}

/*
 * Gather constraints from hose synthesis
 */

static
int synth_hose_class(
  char *method,
  int   hose_color,
  FILE *dat,
  FILE *temp,
  char *group)
{
  char   line[512];
  char  *nonwhite;
  part_t constraints[128];
  int    constraint_n = 0;
  int    color;
  char   start_type[64];
  float  x,y,z, a,b,c, d,e,f, g,h,i;
  int    rc = 0;
  int    hide = 1;
  int    ghost = 0;
  int    group_count = 0;
  char   group_name[512];

  if (group) {
    strcpy(group_name,group);
    group = group_name;
  } else {
    group_name[0] = '\0';
  }

  memset(constraints, 0, 128*sizeof(part_t));

  /* gather up the constraints */

  while (L3fgets(line,sizeof(line), dat)) {
    nonwhite = line + strspn(line," \t");

    skip_rot(line,sizeof(line),dat,temp);

    nonwhite = line + strspn(line," \t");

    if (strncasecmp(nonwhite,"0 GHOST ",strlen("0 GHOST ")) == 0) {
      ghost = 1;
      nonwhite += strlen("0 GHOST ");

      nonwhite += strspn(nonwhite," \t");
    }

    strclean(nonwhite);

    if (sscanf(nonwhite,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s",
        &color, &x,&y,&z, &a,&b,&c, &d,&e,&f, &g,&h,&i, start_type) == 14) {

      if (ishoseconstraint(start_type)) {
        part_t *constr = &constraints[constraint_n];

        if ( ! ldraw_part) {
          if (hide) {
            fputs("0 ",temp);
          }
          fputs(line,temp);
        }
        constr->offset[0] = x;   constr->offset[1] = y;   constr->offset[2] = z;

        constr->orient[0][0] = a;constr->orient[0][1] = b;constr->orient[0][2] = c;
        constr->orient[1][0] = d;constr->orient[1][1] = e;constr->orient[1][2] = f;
        constr->orient[2][0] = g;constr->orient[2][1] = h;constr->orient[2][2] = i;

        strcpy(constraints[constraint_n].type,start_type);
        constraint_n++;
      } else {
        if (group) {
          fprintf(temp,"0 MLCAD BTG %s\n",group);
          group_count++;
        }
        fputs(line,temp);
      }
    } else if (strcasecmp(nonwhite,"0 SYNTH HIDE\n") == 0) {
      fputs(line,temp);
      hide = 1;
    } else if (strcasecmp(nonwhite,"0 SYNTH SHOW\n") == 0) {
      fputs(line,temp);
      hide = 0;
    } else if (strncasecmp(nonwhite,"0 MLCAD BTG ",strlen("0 MLCAD BTG ")) == 0) {
      if (! group) {
        nonwhite += strlen("0 MLCAD BTG ");
        strcpy(group_name,nonwhite);
        group = group_name;
        group[strlen(group)-1] = '\0';
      }
    } else if (strncasecmp(nonwhite,"0 GROUP ",strlen("0 GROUP ")) == 0) {
      fputs(line,temp);
    } else if (strcasecmp(line,"0 SYNTHESIZED BEGIN\n") == 0) {
      rc = skip_synthesized(dat, line, sizeof(line));
      if (rc < 0) {
        break;
      }
    } else if (strncasecmp(nonwhite,"0 SYNTH END\n",strlen("0 SYNTH END\n")) == 0 ||
               strncasecmp(nonwhite,"0 WRITE SYNTH END\n",strlen("0 WRITE SYNTH END\n")) == 0) {
      rc = synth_hose(method,constraint_n,constraints,ghost,group,group_count,hose_color,temp);
      if ( ! ldraw_part ) {
        fputs(line,temp);
      }
      break;
    } else {
      fputs(line,temp);
    }
  }
  return rc;
}

/*
 * Gather up rubber band constraints
 */

static
int synth_band_class(
  char *method,
  int   color,
  FILE *dat,
  FILE *temp,
  char *group)
{
  char line[512];
  char *nonwhite;
  LSL_band_constraint constraints[128];
  int  constraint_n = 0;
  int  rc = 0;
  int  hide = 0;
  int  ghost = 0;

  memset(constraints, 0, sizeof(LSL_band_constraint)*128);

  /* gather up the constraints */

  while (L3fgets(line,sizeof(line), dat)) {
    float x,y,z, a,b,c, d,e,f, g,h,i;
    char start_type[64];
    int t;

    nonwhite = line + strspn(line, " \t");

    skip_rot(line,sizeof(line),dat,temp);

    nonwhite = line + strspn(line, " \t");

    if (strncasecmp(nonwhite,"0 GHOST ",strlen("0 GHOST ")) == 0) {
      ghost = 1;
      nonwhite += strlen("0 GHOST ");

      nonwhite += strspn(nonwhite," \t");
    }

    if (sscanf(nonwhite,"1 %d %f %f %f %f %f %f %f %f %f %f %f %f %s",
        &t, &x,&y,&z, &a,&b,&c, &d,&e,&f, &g,&h,&i, start_type) == 14) {

      if (isbandconstraint(start_type)) {
        part_t *cp = &constraints[constraint_n].part;

        if (hide) {
          fputs("0 ",temp);
        }
        fputs(line,temp);

        constraints[constraint_n].part.attrib = color;
        strcpy(constraints[constraint_n].part.type,start_type);

        cp->offset[0] = x;   cp->offset[1] = y;   cp->offset[2] = z;
        cp->orient[0][0] = a;cp->orient[0][1] = b;cp->orient[0][2] = c;
        cp->orient[1][0] = d;cp->orient[1][1] = e;cp->orient[1][2] = f;
        cp->orient[2][0] = g;cp->orient[2][1] = h;cp->orient[2][2] = i;

        strcpy(constraints[constraint_n].part.type,start_type);
        constraint_n++;
      }
    } else {

      strclean(nonwhite);

      if (strcasecmp(nonwhite,"0 SYNTH INSIDE\n") == 0) {
        fputs(line,temp);
        strcpy(constraints[constraint_n].part.type,"INSIDE");
        constraint_n++;
      } else if (strcasecmp(nonwhite,"0 SYNTH OUTSIDE\n") == 0) {
        fputs(line,temp);
        strcpy(constraints[constraint_n].part.type,"OUTSIDE");
        constraint_n++;
      } else if (strcasecmp(nonwhite,"0 SYNTH CROSS\n") == 0) {
        fputs(line,temp);
        strcpy(constraints[constraint_n].part.type,"CROSS");
        constraint_n++;
      } else if (strcasecmp(nonwhite,"0 SYNTH HIDE\n") == 0) {
        fputs(line,temp);
        hide = 1;
      } else if (strcasecmp(nonwhite,"0 SYNTH SHOW\n") == 0) {
        fputs(line,temp);
        hide = 0;
      } else if (strncasecmp(nonwhite,"0 MLCAD BTG ",strlen("0 MLCAD BTG ")) == 0) {
        if (! group) {
          group = nonwhite + strlen("0 MLCAD BTG ");
          group[strlen(group)-1] = '\0';
        }
        fputs(line,temp);
      } else if (strncasecmp(nonwhite,"0 GROUP ",strlen("0 GROUP ")) == 0) {
      } else if (strcasecmp(line,"0 SYNTHESIZED BEGIN\n") == 0) {
        rc = skip_synthesized(dat, line, sizeof(line));
        if (rc < 0) {
          break;
        }
      } else if (strncasecmp(nonwhite,"0 SYNTH END\n",strlen("0 SYNTH END\n")) == 0 ||
                 strncasecmp(nonwhite,"0 WRITE SYNTH END\n",strlen("0 WRITE SYNTH END\n")) == 0) {

        rc = synth_band(method,constraint_n,constraints,color,temp,ghost,group);
        if ( ! ldraw_part ) {
          fputs(line,temp);
        }
        break;
      } else {
        fputs(line,temp);
      }
    }
  }
  return 0;
}

PRECISION max_bend = 0.05;
PRECISION max_twist = 0.0174;
PRECISION band_res = 1;

//---------------------------------------------------------------------------

char * stripquotes(char *s)
{
  char *p;
  int i;

  // Strip away leading whitespace (spaces and tabs).
  s += strspn(s, " \t");

  // Remove leading quotes
  if (*s == '\"')
    s++;

  // Allocate memory so we can modify the end of the string.
  s = strdup(s);

  // Eliminate trailing whitespace

  for (p = s + (strlen(s)-1); p >= s; p--) {
    if ((*p == ' ') || (*p == '\t')) {
      *p = 0;
    } else {
      break;
    }
  }

  // Remove trailing quotes.
  if ((p = strrchr(s, '\"')) != NULL) {
    *p = 0;
  }

  return(s);
}

#pragma argsused
int main(int argc, char* argv[])
{
  int   dat_argc = 1;
  char *dat_name;
  char *dst_name;
  char *synth_name = NULL;
  char  filename[512];
  FILE *outfile;
  FILE *synthfile;
  FILE *dat;
  char line[512];
  char *nonwhite;
  int   synthcount = 0;
  int   subfiles = 0;
  char *product = NULL;
  char *method = NULL;
  /*
   * Read in the descriptions and constraints for synthesis
   */

  if (parse_descr(argv[0])) {
    return 1;
  }

  /*
   * Output MLCad.ini for LSynth
   */

  if (argc == 2 && strcasecmp(argv[1],"-m") == 0) {
    char path[256];
    char *l,*p;
    int i;

    strcpy(path,argv[0]);

    for (i = 0; i < 2; i++) {
      for (p = path; *p; p++) {
        if (*p == '\\') {
          l = p;
        }
      }
      *l = '\0';
    }

    printf("[LSYNTH]\n");
    printf("%%PATH = \"%s\"\n",path);
    product_ini();
    hose_ini();
    band_ini();
    printf("Tangent Statement: INSIDE = SYNTH INSIDE\n");
    printf("Tangent Statement: OUTSIDE = SYNTH OUTSIDE\n");
    printf("Tangent Statement: CROSS = SYNTH CROSS\n");
    printf("Visibility Statement: SHOW = SYNTH SHOW\n");
    printf("Visibility Statement: HIDE = SYNTH HIDE\n");
    return 1;
  }

  if (argc == 2 && strcasecmp(argv[1],"-p") == 0) {
    printf(argv[0]);
    exit(0);
  }

  printf("LSynth version %s%s by Kevin Clague, kevin_clague@yahoo.com\n",version,beta);
  printf("\n");  // ("                   and Don Heyse\n"); // Much cleaner.

  if (argc == 2 && strcasecmp(argv[1],"-v") == 0) {
    return 1;
  }

  /*
   * Print out help display
   */

  if (argc == 2 && strcasecmp(argv[1],"-h") == 0 || argc == 1) {
    extern void list_hose_types(void);
    extern void list_band_types(void);
    extern void list_band_constraints(void);

    printf("LSynth is an LDraw compatible flexible part synthesizer\n");
    printf("  usage: lsynthcp [-v] [-h] [-m] [-l] [-p] <src> <dst>\n");
    printf("    -v - prints lsynthcp version\n");
    printf("    -h - prints this help message\n");
    printf("    -m - prints out the LSynth portion of the MLcad.ini for using\n");
    printf("         this program\n");
    printf("    -l - format the output as an official ldraw part\n");
    printf("    -p - prints out the full path name of the this executable\n");
    printf("The easiest way to use LSynth is from within MLcad.  You need to\n");
    printf("make additions to MLCad.ini.  Please see Willy Tscager's tutorial\n");
    printf("page http://www.holly-wood.it/mlcad/lsynth-en.html.\n");
    printf("\n");
    printf("To create a flexible part, you put specifications for the part\n");
    printf("directly into your LDraw file, where the part is needed.\n");

    list_products();
    list_hose_types();
    list_hose_constraints();
    list_band_types();
    list_band_constraints();
    return 1;
  }

  if (strcasecmp(argv[1],"-l") == 0) {
    ldraw_part = 1;
    dat_argc++;
  }
  dat_name = argv[dat_argc];
  dst_name = argv[dat_argc+1];

  if (argc < 3) {
    printf("usage: lsynth <input_file> <output_file>\n");
    return 1;
  }

  dat = fopen(dat_name,"r");

  if (dat == NULL) {
    printf("%s: Failed to open file %s for reading\n",argv[0],dat_name);
    return -1;
  }

  outfile = fopen(dst_name,"w");

  if (outfile == NULL) {
    printf("%s: Failed to open file %s for writing\n",argv[0],dst_name);
    return -1;
  }

  /*
   * Scan the input file looking for synthesis specifications
   */

  while (L3fgets(line,sizeof(line), dat)) {

    int t1;

    nonwhite = line + strspn(line, " \t");
    strclean(nonwhite);

    if (strncasecmp(nonwhite,"0 SYNTH BEGIN ",
                    strlen("0 SYNTH BEGIN ")) == 0) {
      char *option;
      char *group;
      char  tmp[256];
      int   color;

      nonwhite += strlen("0 SYNTH BEGIN ");

      synthfile = outfile; // By default synth data goes to main outfile.
      synthcount++; // Count the synth parts.

      // We may be writing synth data to separate subfiles.
      if (subfiles) {
        // Check if the subfile is named in the SYNTH BEGIN line.
        option = strstr(nonwhite, "NAME=");
        if (option) {
          option += strlen("NAME=");
          // Allocate memory for name, strip quotes and trailing spaces.
          synth_name = stripquotes(option);
        } else {
         // Otherwise just number the subfile.  Use lsynthN.ldr for stdout?
          char *p;

          // Remove the extension from dat_name if not already gone.
          if ((p = strrchr(dst_name, '.')) != NULL)
            *p = 0;
          // Build the new subfilename.
          sprintf(filename, "%s%d.ldr", dst_name, synthcount);
          synth_name = strdup(filename);
        }

        synthfile = fopen(synth_name,"w");
        if (synthfile == NULL) {
          printf("%s: Failed to open file %s for writing\n",argv[0],synth_name);
          return -1;
        }

        fputs(line, synthfile);
      }
      option = strstr(nonwhite, "LENGTH=");
      if (option) {
        option += strlen("LENGTH=");
      }
      option = strstr(nonwhite, "UNITS=");
      if (option) {
        option += strlen("UNITS=");
      }
      group = strstr(line, "GROUP=");
      if (group) {
        char *s;
        group += strlen("GROUP=");
        s = group;
        while (*s && *s != ' ' && *s != '\n') {
          s++;
        }
        *s = '\0';
      }

      /* check to see if it is a known synth command */

      if (sscanf(nonwhite,"%s %d",tmp,&color) != 2) {
        return -1;
      }

      product = isproduct(tmp);
      if (product) {
        fprintf(outfile,"0 LPUB PLI BEGIN SUB %s %d\n",product,color);
        method = product_method(product);
      } else {
        method = tmp;
      }

      if ( ! ldraw_part ) {
        fputs(line,outfile);
      }

      if (ishosetype(method)) {
        synth_hose_class(method,color,dat,synthfile,group);

      } else if (isbandtype(method)) {
        synth_band_class(method,color,dat,synthfile,group);

      } else {
        printf("Unknown synthesis type %s\n",nonwhite);
      }

      if (product) {
        fprintf(outfile,"0 LPUB PLI END\n");
        printf("Synthesized %s (%s)\n",product,product_nickname(product));
      } else {
        printf("Synthesized %s\n",method);
      }

      // Close subfile and cleanup

      if (synth_name) {
        fclose(synthfile);
        free(synth_name);
        synth_name = NULL;
      }
    } else {
      float foo,bar;

      if (sscanf(nonwhite,"0 SYNTH HOSE_RES %f %f",&foo,&bar) == 1) {
        max_bend = foo;
        max_twist = bar;
        if ( ! ldraw_part ) {
          fputs(line,outfile);
        }
      } else if (sscanf(nonwhite,"0 SYNTH BAND_RES %f",&foo) == 1) {
        band_res = foo;
        if ( ! ldraw_part ) {
          fputs(line,outfile);
        }
      } else {
        fputs(line,outfile);
      }
    }
  }
  fclose(dat);
  fclose(outfile);

  printf("lynthcp complete\n");
  return 0;
}







































































































































































































































































































































































































































































































































































































































































































































































































































































