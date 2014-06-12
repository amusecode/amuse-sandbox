/* This program prints the sizes of C types interfaced by f90gl */
#include <stdio.h>
#include <GL/gl.h>
main()
{
   printf("Sizes of types in bytes:\n");
   printf("pointers    %d\n",sizeof(int *));
   printf("GLbyte      %d\n",sizeof(GLbyte));
   printf("GLubyte     %d\n",sizeof(GLubyte));
   printf("GLshort     %d\n",sizeof(GLshort));
   printf("GLushort    %d\n",sizeof(GLushort));
   printf("GLint       %d\n",sizeof(GLint));
   printf("GLuint      %d\n",sizeof(GLuint));
   printf("GLenum      %d\n",sizeof(GLenum));
   printf("GLbitfield  %d\n",sizeof(GLbitfield));
   printf("GLsizei     %d\n",sizeof(GLsizei));
   printf("GLfloat     %d\n",sizeof(GLfloat));
   printf("GLclampf    %d\n",sizeof(GLclampf));
   printf("GLdouble    %d\n",sizeof(GLdouble));
   printf("GLclampd    %d\n",sizeof(GLclampd));
   printf("C int       %d\n",sizeof(int));
   printf("C uns. int  %d\n",sizeof(unsigned int));
   printf("C float     %d\n",sizeof(float));
}
