//Wireframe Robot
//  *********   Key inputs: *********************
//
//  Recording Animation for robot involved in football/soccer
//
//	b - begin recording
//
//  e - end recording
//
//  s - save to playback file
//
//  l - load playback file and reset transformation
//
//  p - playback the recording
//
//  r - reset transformation
//
//  ***********************************************
//
//  left arrow key - move object[robot/ball] left
//
//  right arrow key - move object[robot/ball] right
//
//  up arrow key - move ball/robot up
//
//  down arrow key - move ball/robot down
//
//  *************  Mouse functions *****************
//
//  middle - to pick/select elements
//
//  left - left rotation 
//
//  right - right rotation
//
//  ***************************************************

#include "stdafx.h"
#include <gl/glut.h>

/* Height and Radius of different parts of robot and the ball*/
#define TORSO_HEIGHT 5.0
#define TORSO_RADIUS 1.0

#define UPPER_ARM_HEIGHT 2.5
#define LOWER_ARM_HEIGHT 1.0
#define UPPER_ARM_RADIUS  0.5
#define LOWER_ARM_RADIUS  0.35
#define HAND_RADIUS 0.35
#define HAND_HEIGHT 1.0

#define LOWER_LEG_RADIUS  0.35
#define LOWER_LEG_HEIGHT 2.0
#define UPPER_LEG_HEIGHT 3.0
#define UPPER_LEG_RADIUS  0.5
#define ANKLE_RADIUS 0.25
#define ANKLE_HEIGHT 1.0
#define FEET_RADIUS 0.25
#define FEET_HEIGHT 1.0

#define HEAD_HEIGHT 1.5
#define HEAD_RADIUS 1.0
#define BALL_RADIUS 1.0

/* Playback */
#define MAXEVENTS 1000
#define RECORDSIZE 4

//Each Record contains an object id, transformation type = translation(1)/rotation(2), transformation value
int obj_id, trans_type,trans_value1,trans_value2;
// Buffer contains the saved user inputs. Each event is of RECORDSIZE.  
int event_buffer[MAXEVENTS*RECORDSIZE];

// event_ptr is for recording into the event_buffer array.
int event_ptr = 0;

// playback_ptr is for reading/playing back from the event_buffer array.
int playback_ptr=0;

//recordMode and playbackMode are flags used for recording and playback.
int recordMode = 0;
int playbackMode = 0;

//Recorded events are saved to <first nameLastInitial>
FILE *jFile = NULL;
char *fileName = "manasan.txt";
/*End of Playback*/

/*width and height of window*/
int max_w = 500; int max_h = 500;

typedef float point[3];

/* initial joint angles */
static GLfloat theta[19] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
            180.0,0.0,0.0,0.0,180.0,0.0,0.0,0.0,0.0}; 

/*Quadric objects*/
GLUquadricObj *t, *h, *lua, *lla, *rua, *rla, *lll, *rll, *rul, *lul, *rt, *ra, *la, *lh, *rh, *rf, *lf, *ball;

int xpos1 = -5.0,ypos1 = 0.0,zpos1=0.0; /*Robot coordinates*/
int xpos2 = 0.0, ypos2 = 0.0 + BALL_RADIUS /*align ball along with robot*/,zpos2 = 0.0; /*Ball coordinates*/

static int selected = 0;
GLfloat aspectRatio = 0;

void timerFunc(int val);

/*Draw the body*/
void torso()
{
   if (selected == 1){ //if picked change the color
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(1.5, 0.0, 0.0);
   }
   glPushMatrix();
   glRotatef(-90.0, 1.0, 0.0, 0.0);
   gluCylinder(t,TORSO_RADIUS, TORSO_RADIUS, TORSO_HEIGHT,10,10);
   glPopMatrix();
}

/*Draw the head*/
void head()
{
   if (selected == 2){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(1.0, 0.5, 0.0);
   }
   glPushMatrix();
   glTranslatef(0.0, 0.5*HEAD_HEIGHT,0.0);
   glScalef(HEAD_RADIUS, HEAD_HEIGHT, HEAD_RADIUS);
   gluSphere(h,1.0,10,10);
   glPopMatrix();
}

/*Draw the left upper arm*/
void left_upper_arm()
{
   if (selected == 3){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(0.0, 1.5, 1.5);
   }
   glPushMatrix();
   glRotatef(90.0, 1.0, 0.0, 0.0);
   gluCylinder(lua,UPPER_ARM_RADIUS, UPPER_ARM_RADIUS, UPPER_ARM_HEIGHT,10,10);
   glPopMatrix();
}

/*Draw the left lower arm*/
void left_lower_arm()
{
   if (selected == 4){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(1.0, 0.5, 0.0);
   }
   glPushMatrix();
   glRotatef(90.0, 1.0, 0.0, 0.0);
   gluCylinder(lla,LOWER_ARM_RADIUS, LOWER_ARM_RADIUS, LOWER_ARM_HEIGHT,10,10);
   glPopMatrix();
}

/*Draw the left hand*/
void left_hand()
{
   if (selected == 5){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(1.0, 0.5, 0.0);
   }
   glPushMatrix();
   glRotatef(90.0, 1.0, 0.0, 0.0);
   gluSphere(lh,HAND_RADIUS,10,10);
   glPopMatrix();
}

/*Draw the right upper arm*/
void right_upper_arm()
{
   if (selected == 6){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(0.0, 1.5, 1.5);
   }
   glPushMatrix();
   glRotatef(90.0, 1.0, 0.0, 0.0);
   gluCylinder(rua,UPPER_ARM_RADIUS, UPPER_ARM_RADIUS, UPPER_ARM_HEIGHT,10,10);
   glPopMatrix();
}

/*Draw the right lower arm*/
void right_lower_arm()
{
   if (selected == 7){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(1.0, 0.5, 0.0);
   }
   glPushMatrix();
   glRotatef(90.0, 1.0, 0.0, 0.0);
   gluCylinder(rla,LOWER_ARM_RADIUS, LOWER_ARM_RADIUS, LOWER_ARM_HEIGHT,10,10);
   glPopMatrix();
}

/*Draw the right hand*/
void right_hand()
{
   if (selected == 8){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(1.0, 0.5, 0.0);
   }
   glPushMatrix();
   glRotatef(90.0, 1.0, 0.0, 0.0);
   gluSphere(rh,HAND_RADIUS,10,10);
   glPopMatrix();
}

/*Draw the left upper leg*/
void left_upper_leg()
{
   if (selected == 9){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(2.0, 0.0, 0.0);
   }
   glPushMatrix();
   glRotatef(-90.0, 1.0, 0.0, 0.0);
   gluCylinder(lul,UPPER_LEG_RADIUS, UPPER_LEG_RADIUS, UPPER_LEG_HEIGHT,10,10);
   glPopMatrix();
}

/*Draw the left lower leg*/
void left_lower_leg()
{
   if (selected == 10){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(1.0, 0.5, 0.0);
   }
   glPushMatrix();
   glRotatef(-90.0, 1.0, 0.0, 0.0);
   gluCylinder(lll,LOWER_LEG_RADIUS, LOWER_LEG_RADIUS, LOWER_LEG_HEIGHT,10,10);
   glPopMatrix();
}

/*Draw the left ankle*/
void left_ankle()
{
   if (selected == 11){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(1.0, 0.5, 0.0);
   }
   glPushMatrix();
   glRotatef(-90.0, 1.0, 0.0, 0.0);
   gluCylinder(la,ANKLE_RADIUS, ANKLE_RADIUS, ANKLE_HEIGHT,10,10);
   glPopMatrix();
}

/*Draw the left feet*/
void left_feet()
{
   if (selected == 12){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(0.0, 0.0, 1.0);
   }
   glPushMatrix();
   glRotatef(-90.0, 1.0, 0.0, 0.0);
   gluSphere(lf,FEET_RADIUS,10,10);
   glPopMatrix();
}

/*Draw the right upper leg*/
void right_upper_leg()
{
   if (selected == 13){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(2.0, 0.0, 0.0);
   }
   glPushMatrix();
   glRotatef(-90.0, 1.0, 0.0, 0.0);
   gluCylinder(rul,UPPER_LEG_RADIUS, UPPER_LEG_RADIUS, UPPER_LEG_HEIGHT,10,10);
   glPopMatrix();
}

/*Draw the right lower leg*/
void right_lower_leg()
{
   if (selected == 14){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(1.0, 0.5, 0.0);
   }
   glPushMatrix();
   glRotatef(-90.0, 1.0, 0.0, 0.0);
   gluCylinder(rll,LOWER_LEG_RADIUS, LOWER_LEG_RADIUS, LOWER_LEG_HEIGHT,10,10);
   glPopMatrix();
}

/*Draw the right ankle*/
void right_ankle()
{
   if (selected == 15){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(1.0, 0.5, 0.0);
   }
   glPushMatrix();
   glRotatef(-90.0, 1.0, 0.0, 0.0);
   gluCylinder(ra,ANKLE_RADIUS, ANKLE_RADIUS, ANKLE_HEIGHT,10,10);
   glPopMatrix();
}

/*Draw the right feet*/
void right_feet()
{
   if (selected == 16){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(0.0, 0.0, 1.0);
   }
   glPushMatrix();
   glRotatef(-90.0, 1.0, 0.0, 0.0);
   gluSphere(rf,FEET_RADIUS,10,10);
   glPopMatrix();
}

/*Draw the ball*/
void drawBall()
{
   if (selected == 17){
    glColor3f(1.0, 1.0, 1.0);
   }else{
    glColor3f(1.0, 0.0, 1.0);
   }
   glPushMatrix();
   glRotatef(-90.0, 1.0, 0.0, 0.0);
   gluSphere(ball,BALL_RADIUS,10,10);
   glPopMatrix();
}

/*function to draw the objects*/
void draw(GLenum mode)
{
	if(mode == GL_SELECT) glLoadName(1);
    glLoadIdentity();

	glTranslatef(xpos1,ypos1,zpos1);
	glRotatef(theta[1], 0.0, 1.0, 0.0);
    torso();
    glPushMatrix();

	if(mode == GL_SELECT) glLoadName(2);
    glTranslatef(0.0, TORSO_HEIGHT+0.5*HEAD_HEIGHT, 0.0);
	glRotatef(theta[2], 0.0, 1.0, 0.0);
    head();

    glPopMatrix();
    glPushMatrix();
	if(mode == GL_SELECT) glLoadName(3);
    glTranslatef(-(TORSO_RADIUS+UPPER_ARM_RADIUS), 0.9*TORSO_HEIGHT, 0.0);
	glRotatef(theta[3], 1.0, 0.0, 0.0);
    left_upper_arm();

	if(mode == GL_SELECT) glLoadName(4);
    glTranslatef(0.0, -UPPER_ARM_HEIGHT, 0.0);
	glRotatef(theta[4], 0.0, 1.0, 0.0);
    left_lower_arm();

	if(mode == GL_SELECT) glLoadName(5);
	glTranslatef(0.0, -LOWER_ARM_HEIGHT, 0.0);
	glRotatef(theta[5], 0.0, 1.0, 0.0);
    left_hand();

    glPopMatrix();
    glPushMatrix();
	if(mode == GL_SELECT) glLoadName(6);
    glTranslatef(TORSO_RADIUS+UPPER_ARM_RADIUS, 0.9*TORSO_HEIGHT, 0.0);
	glRotatef(theta[6], 1.0,0.0, 0.0);
    right_upper_arm();

	if(mode == GL_SELECT) glLoadName(7);
    glTranslatef(0.0, -UPPER_ARM_HEIGHT, 0.0);
	glRotatef(theta[7], 0.0, 1.0, 0.0);
    right_lower_arm();

	if(mode == GL_SELECT) glLoadName(8);
	glTranslatef(0.0, -LOWER_ARM_HEIGHT, 0.0);
	glRotatef(theta[8], 0.0, 1.0, 0.0);
    right_hand();

    glPopMatrix();
    glPushMatrix();
	if(mode == GL_SELECT) glLoadName(9);
    glTranslatef(-(TORSO_RADIUS+UPPER_LEG_RADIUS), 0.1*UPPER_LEG_HEIGHT, 0.0);
	glRotatef(theta[9], 1.0, 0.0, 0.0);
    left_upper_leg();

	if(mode == GL_SELECT) glLoadName(10);
    glTranslatef(0.0, UPPER_LEG_HEIGHT, 0.0);
	glRotatef(theta[10], 0.0, 1.0, 0.0);
    left_lower_leg();
	
	if(mode == GL_SELECT) glLoadName(11);
    glTranslatef(0.0, LOWER_LEG_HEIGHT, 0.0);
	glRotatef(theta[11], 0.0, 1.0, 0.0);
    left_ankle();

	if(mode == GL_SELECT) glLoadName(12);
	glTranslatef(0.0, ANKLE_HEIGHT, 0.0);
	glRotatef(theta[12], 0.0, 1.0, 0.0);
    left_feet();

    glPopMatrix();
    glPushMatrix();
	if(mode == GL_SELECT) glLoadName(13);
    glTranslatef(TORSO_RADIUS+UPPER_LEG_RADIUS, 0.1*UPPER_LEG_HEIGHT, 0.0);
	glRotatef(theta[13], 1.0, 0.0, 0.0);
    right_upper_leg();

	if(mode == GL_SELECT) glLoadName(14);
    glTranslatef(0.0, UPPER_LEG_HEIGHT, 0.0);
	glRotatef(theta[14], 0.0, 1.0, 0.0);
    right_lower_leg();
	
	if(mode == GL_SELECT) glLoadName(15);
    glTranslatef(0.0, LOWER_LEG_HEIGHT, 0.0);
	glRotatef(theta[15], 0.0, 1.0, 0.0);
    right_ankle();
	
	if(mode == GL_SELECT) glLoadName(16);
	glTranslatef(0.0, ANKLE_HEIGHT, 0.0);
	glRotatef(theta[16], 0.0, 1.0, 0.0);
    right_feet();

    glPopMatrix();

	glLoadIdentity();
	glPushMatrix();
	if(mode == GL_SELECT) glLoadName(17);
	glTranslatef(/*(TORSO_RADIUS+UPPER_LEG_RADIUS*/0.0,-(UPPER_LEG_HEIGHT+LOWER_LEG_HEIGHT+ANKLE_HEIGHT), 0.0);
	glRotatef(theta[17], 0.0, 1.0, 0.0);
	glTranslatef(xpos2,ypos2,zpos2);
    drawBall();

    glPopMatrix();

	if(mode == GL_SELECT) glLoadName(18);
	glPushMatrix();
	glColor3f(1.0, 1.0, 1.0);
	glRectf(9.5,-(UPPER_LEG_HEIGHT+LOWER_LEG_HEIGHT+ANKLE_HEIGHT),10.0,7.0);
	glPopMatrix();

    glFlush();
    glutSwapBuffers();
}

/*display objects on screen*/
void display()
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    draw(GL_RENDER);
}

/*  processHits prints out the contents of the 
 *  selection array.
 */
void processHits(GLint hits, GLuint buffer[]){
  unsigned int i,j;
  GLuint names, *ptr;

  printf("hits = %d\n",hits);
  ptr = (GLuint *) buffer;
  if ( hits == 0 ) selected = 0;
  for (i=0; i<hits; i++){
    names = *ptr;
    ptr+=3;
    selected = *ptr;
	obj_id = selected;
    printf("now selected %d\n",selected);
    for (j=0; j<names; j++){
      if (*ptr==1) {
	     printf("torso\n");
      }
      if (*ptr==2) {
	     printf("head\n");
      }
      if (*ptr==3) {
	     printf("upper arm(left)\n");
      }
	  if (*ptr==4) {
	     printf("lower arm(left)\n");
      }
	  if (*ptr==5) {
	     printf("hand(left)\n");
      }
	  if (*ptr==6) {
	     printf("upper arm(right)\n");
      }
	  if (*ptr==7) {
	     printf("lower arm(right)\n");
      }
	  if (*ptr==8) {
	     printf("hand(right)\n");
      }
	  if (*ptr==9) {
	     printf("upper leg(left)\n");
      }
	  if (*ptr==10) {
	     printf("lower leg(left)\n");
      }
	  if (*ptr==11) {
	     printf("ankle(left)\n");
      }
	  if (*ptr==12) {
	     printf("feet(left)\n");
      }
	  if (*ptr==13) {
	     printf("upper leg(right)\n");
      }
	  if (*ptr==14) {
	     printf("lower leg(right)\n");
      }
	  if (*ptr==15) {
	     printf("ankle(right)\n");
      }
	  if (*ptr==16) {
	     printf("feet(right)\n");
      }
	  if (*ptr==17) {
	     printf("ball\n");
      }
      ptr++;
    }
  }
}

void handleBtnClick(int click)
{
	if(click == 0)
	{
		theta[selected] += 5.0;
        if( theta[selected] > 360.0 ) theta[selected] -= 360.0;
	}
	else
	{
		theta[selected] -= 5.0;
        if( theta[selected] < 360.0 ) theta[selected] += 360.0;
	}    
}

#define SIZE 512

void mouse(int btn, int state, int x, int y)
{
   GLuint selectBuf[SIZE];
   GLint hits;
   GLint viewport[4];

   if (btn == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) 
   {
	    glGetIntegerv (GL_VIEWPORT, viewport);

        glSelectBuffer (SIZE, selectBuf);
        glRenderMode(GL_SELECT);

        glInitNames();
        glPushName(0);

        glMatrixMode (GL_PROJECTION);
        glPushMatrix ();
        glLoadIdentity ();
        /*  create 10x10 pixel picking region near cursor location	*/
        gluPickMatrix ((GLdouble) x, (GLdouble) (viewport[3] - y), 
                  10.0, 10.0, viewport);
		gluOrtho2D (-10.0, 10.0, -10.0*aspectRatio, 10.0*aspectRatio);
		glMatrixMode(GL_MODELVIEW);
        draw(GL_SELECT);
        glMatrixMode (GL_PROJECTION);
        glPopMatrix ();
		glMatrixMode(GL_MODELVIEW);
        glFlush ();
		glutSwapBuffers();
        hits = glRenderMode (GL_RENDER);
        processHits (hits, selectBuf);
        glutPostRedisplay();
   }  
   else if (btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
   {
	   handleBtnClick(0);
	   if(recordMode == 1)
	   {		
			event_buffer[event_ptr++] = obj_id;
			event_buffer[event_ptr++] = 2;
			event_buffer[event_ptr++] = theta[obj_id];
			event_buffer[event_ptr++] = 0;
	   }
	   glutPostRedisplay();
   }
   else if (btn == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
   {
	   handleBtnClick(1);
	   glutPostRedisplay();
	   if(recordMode == 1)
	   {		
			event_buffer[event_ptr++] = obj_id;
			event_buffer[event_ptr++] = 2;
			event_buffer[event_ptr++] = theta[obj_id];
			event_buffer[event_ptr++] = 0;
	   }
   }
}

void myReshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (w <= h)
	{   aspectRatio = (GLfloat) h / (GLfloat) w;
        glOrtho(-10.0, 10.0, -10.0 * aspectRatio,
            10.0 * aspectRatio, -10.0, 10.0);
	}
    else
	{
	   aspectRatio = (GLfloat) w/ (GLfloat) h;
       glOrtho(-10.0 * aspectRatio,
            10.0 * aspectRatio, -10.0, 10.0, -10.0, 10.0);
	}
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

//Reset resets all the transformation values and update the screen.
void reset()
{
	xpos1 = -5;
	xpos2 = 0.0;
	ypos2 = 0.0 + BALL_RADIUS;
	ypos1 = 0.0;
	zpos1 = 0.0;
	zpos2 = 0.0;

	for(int i = 1; i <= 17; i++)
	{
		if(i == 9 || i == 13)
			theta[i] = 180.0;
		else
			theta[i] = 0.0;
	}
	display();
}

void myinit()
{
       // glClearColor(1.0, 1.0, 1.0, 1.0);
		glClearColor (0.0, 0.0, 0.0, 0.0);
        glColor3f(0.0, 0.0, 0.0);
		selected = 0;
/* allocate quadrics with filled drawing style */

        h=gluNewQuadric();
        gluQuadricDrawStyle(h, GLU_LINE);
        t=gluNewQuadric();
        gluQuadricDrawStyle(t, GLU_LINE);
        lua=gluNewQuadric();
        gluQuadricDrawStyle(lua, GLU_LINE);
        lla=gluNewQuadric();
        gluQuadricDrawStyle(lla, GLU_LINE);
        rua=gluNewQuadric();
        gluQuadricDrawStyle(rua, GLU_LINE);
        rla=gluNewQuadric();
        gluQuadricDrawStyle(rla, GLU_LINE);
        lul=gluNewQuadric();
        gluQuadricDrawStyle(lul, GLU_SILHOUETTE);
        lll=gluNewQuadric();
        gluQuadricDrawStyle(lll, GLU_LINE);
        rul=gluNewQuadric();
        gluQuadricDrawStyle(rul, GLU_SILHOUETTE);
        rll=gluNewQuadric();
        gluQuadricDrawStyle(rll, GLU_LINE);
        ra=gluNewQuadric();
        gluQuadricDrawStyle(ra, GLU_LINE);
		la=gluNewQuadric();
        gluQuadricDrawStyle(la, GLU_LINE);	
        rh=gluNewQuadric();
        gluQuadricDrawStyle(rh, GLU_LINE);
		lh=gluNewQuadric();
        gluQuadricDrawStyle(lh, GLU_LINE);
		rf=gluNewQuadric();
        gluQuadricDrawStyle(rf, GLU_LINE);
		lf=gluNewQuadric();
        gluQuadricDrawStyle(lf, GLU_LINE);
		ball=gluNewQuadric();
        gluQuadricDrawStyle(ball, GLU_LINE);
		glEnable(GL_DEPTH_TEST);
}

void keyboard(unsigned char key, int x, int y){
  switch (key){
  case 'q':
  case 'Q':
  case 27: exit(0);
 /*Playback*/
 case 'b': //begin recording. Set recordMode
			recordMode = 1;
			playbackMode = 0;
			printf("Recording...  Press 'e' to end\n");
			event_ptr=0;
			break;

		case 'e': //stop recording. Reset recordMode.
			if(recordMode == 1)
			{
				recordMode = 0;
				printf("Recording stopped.  Press 's' to save to playback file.\n");
			}			
			break;

		case 'l': //Load file. Reset everything and load the contents of the file into the buffer.
			recordMode = 0;
			playbackMode = 0;

			event_ptr=0;
			playback_ptr=0;
			reset();

			printf("Loading file %s\n", fileName);
			
			jFile = fopen(fileName, "r");
			if ( jFile == NULL ) {
					printf("Warning: Could not open %s\n", fileName);		
					playbackMode = 0;
			}
			else {
					// Store the events to event_buffer
					while((fscanf(jFile, "%d ", &event_buffer[event_ptr])) != EOF)
					{
						event_ptr++;
					}
					fclose(jFile);
					playbackMode = 1;
			}
			break;

		case 'r': //Reset everything.
			recordMode = 0;
			playbackMode = 0;

			event_ptr=0;
			playback_ptr=0;
			reset();
			break;

		case 'p': //Playback 
			if(playbackMode==1)
			{
				reset();
				glutTimerFunc(4,timerFunc,1);		
			}
			break;

		case 's': //Save file.
			recordMode = 0;
			playbackMode = 0;

			jFile = fopen(fileName, "w");
			if (jFile == NULL) 
			{
				printf("Warning: Could not open %s\n", fileName);
			}
			else {
				for(int j=0;j<event_ptr;j++)
					fprintf(jFile, "%d ", event_buffer[j]);
				fclose(jFile);
				printf("\nEvents saved in %s\n", fileName);
			}
			playback_ptr=0;
			break;


		case 'L':
			handleBtnClick(0);
			if(recordMode == 1)
			{		
			event_buffer[event_ptr++] = obj_id;
			event_buffer[event_ptr++] = 2;
			event_buffer[event_ptr++] = theta[obj_id];
			event_buffer[event_ptr++] = 0;
			}
			glutPostRedisplay();
			break;
		case 'R':
			handleBtnClick(1);
			if(recordMode == 1)
			{		
			event_buffer[event_ptr++] = obj_id;
			event_buffer[event_ptr++] = 2;
			event_buffer[event_ptr++] = theta[obj_id];
			event_buffer[event_ptr++] = 0;
			}
			glutPostRedisplay();
			break;
  }
}

void special_key(int key, int x, int y){
  switch (key){
  case GLUT_KEY_RIGHT: //Move the Ball or Robot to Right
	  if(selected == 17){
			if(xpos2 > 8.0) 
				xpos2 = 8.0;
			xpos2 += 1.0;
		}
		else{
			if(xpos1 > 5.0) 
				xpos1 = 5.0;
			xpos1 += 1.0;
	   }		
	   
	   if(recordMode == 1)
	   {		
			event_buffer[event_ptr++] = obj_id;
			event_buffer[event_ptr++] = 1;
			if(selected == 17)
			{
				event_buffer[event_ptr++] = xpos2;
				event_buffer[event_ptr++] = ypos2;
			}
			else
			{
				event_buffer[event_ptr++] = xpos1;
				event_buffer[event_ptr++] = ypos1;
			}
	   }
	   glutPostRedisplay();
	   break;
  case GLUT_KEY_LEFT: //Move the Ball or Robot to Left
	  if(selected == 17){
			if(xpos2 < -8.0) 
				xpos2 = -8.0;
			xpos2 -= 1.0;
		}
		else{
			if(xpos1 < -6.0) 
				xpos1 = -6.0;
			xpos1 -= 1.0;
		}
		if(recordMode == 1)
	    {		
			event_buffer[event_ptr++] = obj_id;
			event_buffer[event_ptr++] = 1;
			if(selected == 17)
			{
				event_buffer[event_ptr++] = xpos2;
				event_buffer[event_ptr++] = ypos2;
			}
			else
			{
				event_buffer[event_ptr++] = ypos1;
				event_buffer[event_ptr++] = xpos1;
			}
	    }
		glutPostRedisplay();
	break;
  case GLUT_KEY_DOWN:  //Move the Ball Down
	    if(selected == 17){
			if(ypos2 < -2.0) 
				ypos2 = -2.0;
			ypos2 -= 1.0;
		}
		else
		{
			if(ypos1 < -2.0)
				ypos1 = -2.0;
			ypos1 -= 1.0;
		}
		if(recordMode == 1)
	    {		
			event_buffer[event_ptr++] = obj_id;
			event_buffer[event_ptr++] = 1;
			if(selected == 17)
			{
				event_buffer[event_ptr++] = xpos2;
				event_buffer[event_ptr++] = ypos2;
			}
			else
			{
				event_buffer[event_ptr++] = xpos1;
				event_buffer[event_ptr++] = ypos1;
			}
	     }	
        glutPostRedisplay();
    break;
  case GLUT_KEY_UP: //Move the Ball Up
	   if(selected == 17){
			if(ypos2 > 14.0) 
				ypos2 = 14.0;
			ypos2 += 1.0;
	   }
	   else
	   {
		   if(ypos1 > 1.5)
			   ypos1 = 1.5;
		   ypos1 += 1.0;
	   }
       if(recordMode == 1)
	   {		
			event_buffer[event_ptr++] = obj_id;
			event_buffer[event_ptr++] = 1;
			if(selected == 17)
			{
				event_buffer[event_ptr++] = xpos2;
				event_buffer[event_ptr++] = ypos2;
			}
			else
			{
			   event_buffer[event_ptr++] = xpos1;
               event_buffer[event_ptr++] = ypos1;
			}
	   }
       glutPostRedisplay();
    break;
  }
}

//The timer function
void timerFunc(int val)
{
			// Check if playback_ptr has reached the last event
			if(playback_ptr<event_ptr)
			{	
				obj_id = event_buffer[playback_ptr++];
				trans_type = event_buffer[playback_ptr++];
				trans_value1 = event_buffer[playback_ptr++];
				trans_value2 = event_buffer[playback_ptr++];
				
				// Update the object's transformation value
				switch(obj_id)
				{
				case 1:
                case 2:
				case 3:
				case 4:
				case 5:
				case 6:
				case 7:
				case 8:
				case 9:
				case 10:
				case 11:
				case 12:
				case 13:
				case 14:
				case 15:
				case 16:
					if(trans_type == 2)
						theta[obj_id] = trans_value1;
					else
						xpos1 = trans_value1;
						ypos1 = trans_value2;
					break;
				case 17:
					if(trans_type == 2)
						theta[obj_id] = trans_value1;
					else
						xpos2 = trans_value1;
						ypos2 = trans_value2;
					break;
				default:
					break;
				}

				// Update the screen with the current transformation retrieved from the event_buffer.
				display();

				// Call TimerFunc to read another event
				glutTimerFunc(50, timerFunc, 1);		
			}
			else
			{
				playback_ptr=0;
			}
}

void main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(max_w, max_h);
	glutInitWindowPosition(0,0);	// position of top-left window corner 
    glutCreateWindow("robot");
    myinit();
    glutReshapeFunc(myReshape);
    glutDisplayFunc(display);
    glutMouseFunc(mouse);
	glutKeyboardFunc(keyboard);
    glutSpecialFunc(special_key);

    glutMainLoop();
}
