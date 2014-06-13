/* Program to create a preprocessor include file that identifies
   the OpenGL and GLU version number and supported extensions.
   Based on an idea and code donated by David Scriven, U. of British Columbia.
*/

#if defined(MACOSX)
#include <OpenGL/gl.h>
#include <OpenGL/glext.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include <stdio.h>
#include <string.h>

char* f90glext[] =
{"GL_EXT_abgr",
 "GL_EXT_blend_color",
 "GL_EXT_blend_logic_op",
 "GL_EXT_blend_minmax",
 "GL_EXT_blend_subtract",
 "GL_EXT_color_table",
 "GL_EXT_paletted_texture",
 "GL_EXT_point_parameters",
 "GL_EXT_polygon_offset",
 "GL_EXT_rescale_normal",
 "GL_EXT_shared_texture_palette",
 "GL_EXT_subtexture",
 "GL_EXT_texture3D",
 "GL_EXT_texture_object",
 "GL_EXT_vertex_array",
 "GL_SGIS_multisample",
 "GL_MESA_resize_buffers",
 "GL_MESA_window_pos",
 "GL_ARB_point_parameters"
 };
 
 char* Level[] =
 { 
   "1.0",
   "1.1",
   "1.2",
   "1.3",
   "1.4",
   "1.5",
   "1.6",
   "1.7",
   "1.8",
   "1.9"
 };
 
int number_of_extensions = 19;

int
main(int argc, char *argv[])
{
   char *displayName = ":0";
   const char *glVersion;
   const char *gluVersion;
   const GLubyte *glExtensions;
   char VerGL[8];
   char VerGLU[8];
   int  GL_Level = 0;
   int  GLU_Level = 0;
   FILE *glopts;

   int i,j,k;

   int iwin;

   glutInit(&argc,argv);
   glutInitDisplayMode(GLUT_RGB);
   iwin = glutCreateWindow("a");

   glopts  = fopen("glopt.h","w");
   if (glopts == NULL)
   {
      fprintf(stderr,"error opening file glopt.h.\n"); 
      return -1;
   }
  
/* Get the versions and extensions with GetStrings rather than from
   which macros are defined in gl.h, because gl.h is sometimes wrong */

   glVersion = (const char *) glGetString(GL_VERSION);
   gluVersion = (const char *) gluGetString(GLU_VERSION);
   glExtensions = (const GLubyte *) glGetString(GL_EXTENSIONS);

   for (k = 0; k < number_of_extensions; k++)
   {
     const GLubyte *start = glExtensions;
     GLubyte *where, *terminator;
     const char *extension = f90glext[k];
     int present = 0;  
     
     for (;;) 
     {
       where = (GLubyte *) strstr((const char *) start, extension);
       if (!where)
         break;
       terminator = where + strlen(extension);
       if (where == start || *(where - 1) == ' ')
       if (*terminator == ' ' || *terminator == '\0')
       {
         present = 1;
         break;
       }
       start = terminator;
     }
     if(present)fprintf(glopts,"#define %s\n",extension);
     else fprintf(glopts,"#undef %s\n",extension);
   }
  
   strncpy(VerGL, glVersion,3);
   *(VerGL+3) = '\0';
   strncpy(VerGLU, gluVersion,3);
   *(VerGLU+3) = '\0';

   GL_Level = 2; /* default to the highest level supported by f90gl */
   for(i=0; i < 10; i++)
   {
      char *p = strstr(VerGL, Level[i]);
      if(p != NULL) {GL_Level = i;}
      p = strstr(VerGLU, Level[i]);   
      if(p != NULL) {GLU_Level = i;}   
   }
   
   for(j=0; j <= GL_Level; j++)
     fprintf(glopts,"#define OPENGL_1_%d\n",j);
   
   if(GLU_Level < 2)
     fprintf(glopts,"#define GLU_1_1_TESS\n");
     
   fclose(glopts);
   
   return 0;
}

