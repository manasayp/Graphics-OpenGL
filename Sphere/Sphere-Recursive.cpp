// Sphere-Recursive.cpp : Defines the entry point for the console application.
//
// This program generates a sphere by recursive subdivision
// of a tetrahedron.


#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gl/glut.h>
//#include "C:\Users\David\Downloads\glut-3.7.6-bin\glut-3.7.6-bin\glut.h"
#include <math.h>

void myinit();
void display(void);

void triangle(GLfloat *a, GLfloat *b, GLfloat *c);
void tetra(GLfloat *a, GLfloat *b, GLfloat *c, GLfloat *d);
void divide_tetra(GLfloat *a, GLfloat *b, GLfloat *c, GLfloat *d,int m);
void divide_triangle(GLfloat *a, GLfloat *b, GLfloat *c, int m);

void myReshape(int w, int h);

/*int mountain = 0;
float min_dist, max_dist;
float scale = 0.1;*/

/* Number of subdivision  */
int n = 6;

/*  Vertex definitions  */

GLfloat v[4][3] = {{0.0, 0.0, 1.0}, {0.0, 0.942809, -0.33333}, {-0.816497, -0.471405, -0.333333}, {0.816497, -0.471405, -0.333333}};
GLfloat colors[4][3] = {{1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0}, {1.0,1.0,0.0}};

/*  Update the viewport when it the window is resized  */
void myReshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if (w <= h)
        glOrtho(-1.25, 1.25, -1.25 * (GLfloat) h / (GLfloat) w,
            1.25 * (GLfloat) h / (GLfloat) w, -10.0, 10.0);
    else
        glOrtho(-1.25 * (GLfloat) w / (GLfloat) h,
		1.25 * (GLfloat) w / (GLfloat) h, -1.25, 1.25, -10.0, 10.0);
	//glOrtho(-2.0,2.0,-2.0,2.0 ,-10.0,10.0);

    glMatrixMode(GL_MODELVIEW);
    glutPostRedisplay();
}

/*  Initialize the global settings.  */
void myinit()
{
   glClearColor(1.0, 1.0, 1.0, 1.0);
   glColor3f(1.0, 0.0, 0.0);
}

/*  Main  */
void main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(500,500);
	glutCreateWindow("Wireframe Sphere");
	myinit();

	glutReshapeFunc(myReshape);
	glutDisplayFunc(display);

	glEnable(GL_DEPTH_TEST);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glutMainLoop();
}

/*  Draw the current object.  */

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glBegin(GL_LINE_LOOP);
	    divide_tetra(v[0],v[1],v[2],v[3],n);
	glEnd();
	glFlush();
    glutSwapBuffers();
}

/*  Draw the given triangle.  */

void triangle(GLfloat *a, GLfloat *b, GLfloat *c)
{
	glColor3fv(colors[0]);
	glVertex3fv(a);
	glColor3fv(colors[1]);
	glVertex3fv(b);
	glColor3fv(colors[2]);
	glVertex3fv(c);
}

/*Caculate the normal(unit vector)*/
void normalize(GLfloat *v)
{
    double d=0.0;
    int i;
    for(i=0;i<3;i++) d+=v[i]*v[i];
    d=sqrt(d);
    if(d>0.0) for(i=0;i<3;i++) v[i]/=d;
}

void divide_triangle(GLfloat *a, GLfloat *b, GLfloat *c, int m) 
{
    /* triangle subdivision using vertex numbers */
    GLfloat mid[3][3];     
	int j;     
	if(m>0)    
	{         
		for(j=0; j<3; j++) mid[0][j]=(a[j]+b[j])/2; normalize(mid[0]);        
		for(j=0; j<3; j++) mid[1][j]=(a[j]+c[j])/2; normalize(mid[1]);        
		for(j=0; j<3; j++) mid[2][j]=(b[j]+c[j])/2; normalize(mid[2]);    
		divide_triangle(a, mid[1], mid[0], m-1);         
		divide_triangle(c, mid[2], mid[1], m-1);         
		divide_triangle(b, mid[0], mid[2], m-1); 
		divide_triangle(mid[0],mid[1],mid[2], m-1);
	}     
	else triangle(a,b,c); /* draw triangle at end of recursion */
}

/* Subdivide the given tetrahedron faces recursively into multiple triangles */
void divide_tetra(GLfloat *a, GLfloat *b, GLfloat *c, GLfloat *d,int m)
{
	divide_triangle(a,b,c,m);
	divide_triangle(a,d,b,m);
	divide_triangle(d,c,b,m);
	divide_triangle(a,c,d,m);
}
