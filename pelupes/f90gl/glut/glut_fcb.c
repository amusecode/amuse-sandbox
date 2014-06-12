#include "fppr.inc"

#ifdef FREEGLUT

/* #include "glutint.h" */
#include <malloc.h>
#include "GL/glut.h"
#include "stdio.h"

/* environments that don't want to use APIENTRY */

#ifndef WIN32
#define APIENTRY
#endif

#ifdef ABS32
#define APIENTRY
#endif

#ifdef PGI32
#define APIENTRY
#endif

#ifdef LF9X
#define APIENTRY
#endif

#ifdef IFORT
#undef APIENTRY
#define APIENTRY
#endif

#include "glutf90.h"


/* FCB stands for Fortran CallBack. */

/* There is only one idleFunc, menuStateFunc, and menuStatusFunc, so they
   can be saved in the wrappers for Fortran rather than the C structures. */

/* Set a Fortran callback function. */
typedef struct _GLUTwindow GLUTwindow;
struct _GLUTwindow{
  int type;
  int num;
  int value;
  GLUTwindow *next;

  union{
  GLUTdisplayFCB fdisplay;  /* Fortran display  */
  GLUTdisplayFCB fdisplay_ov;  /* Fortran display  */
  GLUTreshapeFCB freshape;  /* Fortran reshape  */
  GLUTmouseFCB fmouse;  /* Fortran mouse  */
  GLUTmotionFCB fmotion;  /* Fortran motion  */
  GLUTpassiveFCB fpassive;  /* Fortran passive  */
  GLUTentryFCB fentry;  /* Fortran entry  */
  GLUTkeyboardFCB fkeyboard;  /* Fortran keyboard  */
  GLUTkeyboardFCB fkeyboardUp;  /* Fortran keyboard up */
  GLUTwindowStatusFCB fwindowStatus;  /* Fortran window status */
  GLUTvisibilityFCB fvisibility;  /* Fortran visibility */
  GLUTspecialFCB fspecial;  /* special key */
  GLUTspecialFCB fspecialUp;  /* special key up */
  GLUTbuttonBoxFCB fbuttonBox;  /* button box */
  GLUTdialsFCB fdials;  /* dials */
  GLUTspaceMotionFCB fspaceMotion;  /* Spaceball motion */
  GLUTspaceRotateFCB fspaceRotate;  /* Spaceball rotate */
  GLUTspaceButtonFCB fspaceButton;  /* Spaceball button */
  GLUTtabletMotionFCB ftabletMotion;  /* tablet motion */
  GLUTtabletButtonFCB ftabletButton;  /* tablet button */
#ifdef _WIN32
  GLUTjoystickFCB fjoystick;  /* joystick */
#endif
  GLUTselectFCB fselect; 
  GLUTtimerFCB ffunc;
  };
}; 

int __timernum=0;

GLUTwindow *head= (GLUTwindow *) NULL;
GLUTwindow *end= (GLUTwindow *) NULL;

/* this initialises a node, allocates memory for the node, and returns   */
/* a pointer to the new node. Must pass it the node details, name and id */
GLUTwindow *initnode( int type, int num )
{
   GLUTwindow *ptr;
   ptr = (GLUTwindow *) malloc( sizeof(GLUTwindow) );
   if( ptr == NULL )                       
       return (GLUTwindow *) NULL;       
   else {                                
       ptr->num = num;
       ptr->type = type;
       ptr->next=NULL;                   
       return ptr;                       
   }
}

void add( GLUTwindow *new ) 
{
   if( head == NULL ) {  
     head = new;        
     end = new;
     new->next=NULL;
   } else {
     new->next = end->next;  
     end->next = new;        
     end = new;
   }              
}

void delnode( GLUTwindow *node )
{
  GLUTwindow *temp;
  if(node == head) {
    head=node->next;
    node->next=end->next;
    end->next=node;
  } else {
    temp=head;
    while(temp->next!=node)
      temp=temp->next;
    temp->next=node->next;
    node->next=end->next;
    end->next=node;
  }
}

GLUTwindow *search( GLUTwindow *ptr, int type, int num )
{
    while(1) { 
       if( ptr == NULL )  break;                                
       if( num == ptr->num && type == ptr->type) break;
       if( ptr == end ) {
	 ptr=NULL;
	 return ptr;  
        }
       ptr = ptr->next;                                
    }
    return ptr;                                 
}                                               

void
__glutSetTimerFCB(int which, void *func, int value, int *num)
{
  int type;
  GLUTwindow *__glutCurrentWindow;

  type=which;  
  *num=__timernum++;
  __glutCurrentWindow=search(head,type, *num);
  if(__glutCurrentWindow==NULL) {
   __glutCurrentWindow=initnode(type, *num);
   add(__glutCurrentWindow);
  }
  __glutCurrentWindow->value=value;
  __glutCurrentWindow->ffunc = (GLUTtimerFCB) func;  
}


/*void APIENTRY*/
void
__glutSetFCB(int which, void *func)
{
  int type,num;
  GLUTwindow *__glutCurrentWindow;


  type=which;
  switch (which) {
   case GLUT_FCB_SELECT:
    num=glutGetMenu();
   break;
   default:
    num=glutGetWindow();
   } 
  __glutCurrentWindow=search(head,type,num);
  if(__glutCurrentWindow==NULL) {
   __glutCurrentWindow=initnode(type,num);
   add(__glutCurrentWindow);
  }
   
  switch (which) {
  case GLUT_FCB_DISPLAY:
    __glutCurrentWindow->fdisplay = (GLUTdisplayFCB) func;
    break;
  case GLUT_FCB_RESHAPE:
    __glutCurrentWindow->freshape = (GLUTreshapeFCB) func;
    break;
  case GLUT_FCB_MOUSE:
    __glutCurrentWindow->fmouse = (GLUTmouseFCB) func;
    break;
  case GLUT_FCB_MOTION:
    __glutCurrentWindow->fmotion = (GLUTmotionFCB) func;
    break;
  case GLUT_FCB_PASSIVE:
    __glutCurrentWindow->fpassive = (GLUTpassiveFCB) func;
    break;
  case GLUT_FCB_ENTRY:
    __glutCurrentWindow->fentry = (GLUTentryFCB) func;
    break;
  case GLUT_FCB_KEYBOARD:
    __glutCurrentWindow->fkeyboard = (GLUTkeyboardFCB) func;
    break;
  case GLUT_FCB_KEYBOARD_UP:
    __glutCurrentWindow->fkeyboardUp = (GLUTkeyboardFCB) func;
    break;
  case GLUT_FCB_WINDOW_STATUS:
    __glutCurrentWindow->fwindowStatus = (GLUTwindowStatusFCB) func;
    break;
  case GLUT_FCB_VISIBILITY:
    __glutCurrentWindow->fvisibility = (GLUTvisibilityFCB) func;
    break;
  case GLUT_FCB_SPECIAL:
    __glutCurrentWindow->fspecial = (GLUTspecialFCB) func;
    break;
  case GLUT_FCB_SPECIAL_UP:
    __glutCurrentWindow->fspecialUp = (GLUTspecialFCB) func;
    break;
  case GLUT_FCB_BUTTON_BOX:
    __glutCurrentWindow->fbuttonBox = (GLUTbuttonBoxFCB) func;
    break;
  case GLUT_FCB_DIALS:
    __glutCurrentWindow->fdials = (GLUTdialsFCB) func;
    break;
  case GLUT_FCB_SPACE_MOTION:
    __glutCurrentWindow->fspaceMotion = (GLUTspaceMotionFCB) func;
    break;
  case GLUT_FCB_SPACE_ROTATE:
    __glutCurrentWindow->fspaceRotate = (GLUTspaceRotateFCB) func;
    break;
  case GLUT_FCB_SPACE_BUTTON:
    __glutCurrentWindow->fspaceButton = (GLUTspaceButtonFCB) func;
    break;
  case GLUT_FCB_TABLET_MOTION:
    __glutCurrentWindow->ftabletMotion = (GLUTtabletMotionFCB) func;
    break;
  case GLUT_FCB_TABLET_BUTTON:
    __glutCurrentWindow->ftabletButton = (GLUTtabletButtonFCB) func;
    break;
#ifdef _WIN32
  case GLUT_FCB_JOYSTICK:
    __glutCurrentWindow->fjoystick = (GLUTjoystickFCB) func;
    break;
#endif
  case GLUT_FCB_OVERLAY_DISPLAY:
/*    __glutCurrentWindow->overlay->fdisplay = (GLUTdisplayFCB) func; */
    __glutCurrentWindow->fdisplay_ov = (GLUTdisplayFCB) func;
    break;
  case GLUT_FCB_SELECT:
    __glutCurrentWindow->fselect = (GLUTselectFCB) func;
    break;
/*  case GLUT_FCB_TIMER:
    __glutNewTimer->ffunc = (GLUTtimerFCB) func;
    break;
*/
  }
}

/* Get a Fortran callback function. */

void*
__glutGetTimerFCB(int which,int num, int *value)
{
  int type;
  GLUTwindow *__glutCurrentWindow;

  type=which;  
  *value=0;
  __glutCurrentWindow=search(head,type,num);
  if(__glutCurrentWindow==NULL) return NULL;
  *value=__glutCurrentWindow->value;
  return (void *) __glutCurrentWindow->ffunc;
}

/* void* APIENTRY*/
void*
__glutGetFCB(int which)
{
  int num;
  int type;
  GLUTwindow *__glutCurrentWindow;

  type=which;
  switch (which) {
   case GLUT_FCB_SELECT:
    num=glutGetMenu();
   break;
   default:
    num=glutGetWindow();
   } 

  __glutCurrentWindow=search(head,type,num);
  if(__glutCurrentWindow==NULL) return NULL;

  switch (which) {
  case GLUT_FCB_DISPLAY:
    return (void *) __glutCurrentWindow->fdisplay;
  case GLUT_FCB_RESHAPE:
    return (void *) __glutCurrentWindow->freshape;
  case GLUT_FCB_MOUSE:
    return (void *) __glutCurrentWindow->fmouse;
  case GLUT_FCB_MOTION:
    return (void *) __glutCurrentWindow->fmotion;
  case GLUT_FCB_PASSIVE:
    return (void *) __glutCurrentWindow->fpassive;
  case GLUT_FCB_ENTRY:
    return (void *) __glutCurrentWindow->fentry;
  case GLUT_FCB_KEYBOARD:
    return (void *) __glutCurrentWindow->fkeyboard;
  case GLUT_FCB_KEYBOARD_UP:
    return (void *) __glutCurrentWindow->fkeyboardUp;
  case GLUT_FCB_WINDOW_STATUS:
    return (void *) __glutCurrentWindow->fwindowStatus;
  case GLUT_FCB_VISIBILITY:
    return (void *) __glutCurrentWindow->fvisibility;
  case GLUT_FCB_SPECIAL:
    return (void *) __glutCurrentWindow->fspecial;
  case GLUT_FCB_SPECIAL_UP:
    return (void *) __glutCurrentWindow->fspecialUp;
  case GLUT_FCB_BUTTON_BOX:
    return (void *) __glutCurrentWindow->fbuttonBox;
  case GLUT_FCB_DIALS:
    return (void *) __glutCurrentWindow->fdials;
  case GLUT_FCB_SPACE_MOTION:
    return (void *) __glutCurrentWindow->fspaceMotion;
  case GLUT_FCB_SPACE_ROTATE:
    return (void *) __glutCurrentWindow->fspaceRotate;
  case GLUT_FCB_SPACE_BUTTON:
    return (void *) __glutCurrentWindow->fspaceButton;
  case GLUT_FCB_TABLET_MOTION:
    return (void *) __glutCurrentWindow->ftabletMotion;
  case GLUT_FCB_TABLET_BUTTON:
    return (void *) __glutCurrentWindow->ftabletButton;
  case GLUT_FCB_JOYSTICK:
#ifdef _WIN32
    return (void *) __glutCurrentWindow->fjoystick;
#else
    return NULL;
#endif
  case GLUT_FCB_OVERLAY_DISPLAY:
/*    return (void *) __glutCurrentWindow->overlay->fdisplay; */
    return (void *) __glutCurrentWindow->fdisplay_ov;
  case GLUT_FCB_SELECT:
    return (void *) __glutCurrentWindow->fselect;
/*
  case GLUT_FCB_TIMER:
    return (void *) __glutTimerList->ffunc;
*/
  default:
    return NULL;
  }
}

#endif
