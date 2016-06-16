// BSpline_Curve.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
using namespace std;
/*
**
**  Simple interactive curve drawing program.
**
**  The following keyboard/menu commands are used to control the
**  program:
**  "Ctrl Pt ON/OFF", 			'0',
**  "Ctrl Polygon ON/OFF", 		'1',
**  "BSpline ON/OFF", 			'2',
**	"Surface ON/OFF", 			'7',	
**	"lines/polygons",			'9',
**	"Flat/Smooth",				'f',
**	"Texture ON/OFF",			't',
**	"Lightning ON/OFF",			'l',
**  "Select ON/OFF", 			'3',
**  "Delete ON/OFF", 			'4',
**   "Save to File", 			'5',
**  "Load from File", 			'6',
**  "Clear", 					'8',
**   "quit", 					'q',
	
**
*/
#include "gl/glut.h"

typedef enum
{
   BSPLINE
} curveType;

void keyboard(unsigned char key, int x, int y);

/* Colors to draw them in */
GLfloat colors[][4] = {{1.0,0.5,1.0,0.5},{1.0,0.0,0.0,0.5},{1.0,1.0,0.0,0.5}};

#define MAX_CPTS  75            /* Fixed maximum number of control points */
#define MAX_KNOTS MAX_CPTS+5    /*Maximum Knots*/
#define MAX_BPTS MAX_KNOTS*50   /*Maximum BSpline points on the curve*/
#define PI 3.14159

GLfloat knot[MAX_KNOTS];        /*computed knot values*/
GLfloat cpts[MAX_BPTS][3];      /*control points entered by user*/
GLfloat BSpline[MAX_BPTS][3];   /*Holds the BSpline Curve points*/

/*Stores the vertex and triangle information when the curve is 
rotated from 0 to 360 degree to compute the surface*/
float vertex_normal[(MAX_BPTS+1)*36][3];
float tri_normal[(MAX_BPTS+1)*36*2][3];
float vertex[(MAX_BPTS+1)*36][3];
int triangle[(MAX_BPTS+1)*36*2][3];

int ncpts = 0;      /*Number of control points*/
static int width = 500, height = 500;           /* Window width and height */

/* Flags to Enable/Disable features*/
int ctrlPt_On = 0;						/*Control Point*/
int ctrlPoly_On = 0;					/*Control Polygon*/
int BSpline_On = 0;						/*Curve*/
int Selection_On = 0;					/*Select Control Point*/
int Delete_On = 0;						/*Delete Ctrl Point*/
int BSurface_On = 0;					/*Surface Mesh*/
int drawLines = 0;							/*Wireframe/Polygon*/
int lightingEnabled = 0;				/*Enable Light*/
int smoothEnabled = 0;					/*Flat/Smooth Shading*/
int texEnabled = 0;						/*Texture*/
int fogEnabled = 0;
int depthEnabled = 0;					/*Hidden Surface Removal*/
int texture = 0;

// remember the moving control pt name
int SelectedCtrlPt = -1;

// set up a flag to ensure only one control point is selected
int choice = 0;

int m , num_knots;
int numBSplinePts = -1;
int num_vertices = -1;
int num_traingles = -1;

FILE *jFile = NULL;
char *fileName = "bspline.txt";

/*----------------------------------------------------------------------*/
/* 
** Materials setup.
*/
typedef struct materialStruct {
   GLfloat ambient[4];
   GLfloat diffuse[4];
   GLfloat specular[4];
   GLfloat shininess;
} materialStruct;

materialStruct colorCubeMaterials = {
	{ 0.2F, 0.0F, 0.0F, 1.0F },
    { 0.8F, 0.0F, 0.0F, 1.0F },
    { 1.0F, 1.0F, 1.0F, 1.0F },
    {100.0F}
};

materialStruct *currentMaterials = &colorCubeMaterials;

void materials( materialStruct *materials)
{
    /* define material proerties for front face of all polygons */
    glMaterialfv(GL_FRONT, GL_AMBIENT, materials->ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, materials->diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, materials->specular);
    glMaterialf(GL_FRONT, GL_SHININESS, materials->shininess);
}

/*----------------------------------------------------------------------*/
/*
** Lighting setup.
*/
typedef struct lightingStruct {
    GLfloat ambient[4];
    GLfloat diffuse[4];
    GLfloat specular[4];
} lightingStruct;

lightingStruct whiteLighting = {
    {0.0F, 0.0F, 0.0F, 1.0F},
    {1.0F, 1.0F, 1.0F, 1.0F},
    {1.0F, 1.0F, 1.0F, 1.0F}
};
lightingStruct colorCubeLighting = {
    {0.2F, 0.0F, 0.0F, 1.0F},
    {0.0F, 1.0F, 0.0F, 1.0F},
    {0.0F, 0.0F, 1.0F, 1.0F}
};

lightingStruct *currentLighting = &colorCubeLighting;

void lighting(lightingStruct *lightSettings)
{
    /* set up light 0 ambient, diffuse, specular, and spotlight */
    glLightfv(GL_LIGHT0, GL_AMBIENT, lightSettings->ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightSettings->diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSettings->specular);

    glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 0.0F);
    glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 180.0F);
    glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1.0F);
    glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.0F);
    glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0.0F);

    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
}

/*----------------------------------------------------------------------*/
/*
** Light position setup.
*/
void setLightPos()
{
    //GLfloat light0_pos[4] = { 0.90F, 0.90F, 2.25F, 0.00F };
	GLfloat light0_pos[] = {0.0, 0.0, 10.0, 1.0};
	GLfloat light0_spotDir[3] = {0.0F, 0.0F, -1.0F};

    glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light0_spotDir);
}

/*----------------------------------------------------------------------*/
/*
** Initial texture settings.
*/

/*Computing Texture*/
#define imageWidth 64
#define imageHeight 64
GLubyte checkImage[imageHeight][imageWidth][3];
GLubyte image[3*imageWidth*imageHeight];
GLfloat texpts[2][2][2] = {{{0.0, 0.0}, {0.0, 1.0}}, 
                        {{1.0, 0.0}, {1.0, 1.0}}};

void makeCheckImage(void)
{
   int i, j, c;
    
   for (i = 0; i < imageHeight; i++) {
      for (j = 0; j < imageWidth; j++) {
         c = (((((i&0x8)==0)^((j&0x8))==0))*255);
         checkImage[i][j][0] = (GLubyte) c;
         checkImage[i][j][1] = (GLubyte) c;
         checkImage[i][j][2] = (GLubyte) c;
      }
   }
}

void makeImage(void)
{
   int i, j;
   float ti, tj;
   
   for (i = 0; i < imageWidth; i++) {
      ti = 2.0*3.14159265*i/imageWidth;
      for (j = 0; j < imageHeight; j++) {
         tj = 2.0*3.14159265*j/imageHeight;
         image[3*(imageHeight*i+j)] = 
              (GLubyte) 127*(1.0+sin(ti));
         image[3*(imageHeight*i+j)+1] = 
              (GLubyte) 127*(1.0+cos(2*tj));
         image[3*(imageHeight*i+j)+2] = 
              (GLubyte) 127*(1.0+cos(ti+tj));
      }
   }
}

void initTexture()
{
	if(texture)
	{
		makeImage();
	    glTexImage2D(GL_TEXTURE_2D,0,3,64,64,0,GL_RGB,GL_UNSIGNED_BYTE, image); 
	}
	else
	{
		makeCheckImage();
		glTexImage2D(GL_TEXTURE_2D,0,3,64,64,0,GL_RGB,GL_UNSIGNED_BYTE, checkImage);
		//gluBuild2DMipmaps(GL_TEXTURE_2D,3,64,64,GL_RGB,GL_UNSIGNED_BYTE, image);
	}
    glPixelStorei(GL_UNPACK_ALIGNMENT,1); 
    glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);    
	glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);    
	glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);    
	glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);    
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
}

/*----------------------------------------------------------------------*/
/*
** Set state according to user interaction.
*/
void userSettings(void) 
{
    lighting(currentLighting);
    materials(currentMaterials);
	
	if (lightingEnabled) {
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
    } else {
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
    }

    if (smoothEnabled) {
		glShadeModel(GL_SMOOTH);
    } else {
		glShadeModel(GL_FLAT);
    }

    if (texEnabled) {
		glEnable(GL_TEXTURE_2D);
    } else {
		glDisable(GL_TEXTURE_2D);
    }

    if (depthEnabled) {
	glEnable(GL_DEPTH_TEST); /* Enable hidden--surface--removal */
    } else {
	glDisable(GL_DEPTH_TEST);
    }
}

/*Calculate the knots which determines shape of curve*/
static void calculateKnot()
{
	if(ncpts > 4)
	{
		m = ncpts-1;
		num_knots = m+4;

		for(int i = 0; i <= 3; i++)
		{
			knot[i] = 0;
		}
		for(int i = 4; i <= m; i++)
		{
			knot[i] = i-3;
		}
		for(int i = m+1; i <= m+4; i++)
		{
			knot[i] = m-2;
		}
		cout << endl;
	}
}

/*Blending Function Calculation*/
float CoxdeBoor(int i , int p, float t)
{
	float left,right;

	if(p == 1){
		if((knot[i] < knot[i+1]) && (knot[i] <= t) && (t < knot[i+1]))
			return 1.0;
		else
			return 0.0;
	}
	else
	{
		if(knot[i+p-1] - knot[i] != 0.0)
			left = (CoxdeBoor(i, p-1, t) * (t - knot[i])) / (knot[i+p-1] - knot[i]);
		else
			left = 0.0;

		if(knot[i+p] - knot[i+1] != 0.0)
			right = (CoxdeBoor(i+1, p-1, t) * (knot[i+p] - t))/ (knot[i+p] - knot[i+1]);
		else
			right = 0.0;
		return left+right;
	}
}

/*Computes the BSpline curve*/
static void calculateBSpline()
{
	GLfloat B0,B1,B2,B3;

	numBSplinePts = -1;
	for(int i = 3; i < num_knots-3 ; i++)
	{
		cout << cpts[i][0] << " " << cpts[i][1] << endl;
		for(float t = knot[i] ; t < knot[i+1]; t+= 0.02)
		{
			B0 = CoxdeBoor(i , 4, t);
			B1 = CoxdeBoor(i-1 , 4, t);
			B2 = CoxdeBoor(i-2 , 4, t);
			B3 = CoxdeBoor(i-3 , 4, t);
			numBSplinePts += 1;
			BSpline[numBSplinePts][0] = cpts[i][0]*B0 + cpts[i-1][0]*B1 + cpts[i-2][0]*B2 + cpts[i-3][0]*B3;
			BSpline[numBSplinePts][1] = cpts[i][1]*B0 + cpts[i-1][1]*B1 + cpts[i-2][1]*B2 + cpts[i-3][1]*B3;
		}
	}
	numBSplinePts += 1;
	BSpline[numBSplinePts][0] = cpts[ncpts-1][0];
	BSpline[numBSplinePts][1] = cpts[ncpts-1][1];
}

/*Draws the Curve*/
void drawSpline()
{
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINE_STRIP);
	for(int i = 0; i <= numBSplinePts; i++)
	{
         glVertex2fv(BSpline[i]);
	}
	glEnd();
}

/*Drwas the surface*/
void drawSurface()
{
	int V1, V2, V3;
	if(drawLines) //Filled
	{
		for(int i = 0; i <= num_traingles; i++)
		{
			//glNormal3fv(tri_normal[i]);
			glBegin(GL_POLYGON);
			V1 = triangle[i][0];
			V2 = triangle[i][1];
			V3 = triangle[i][2];
			glColor3fv(colors[0]);
			glTexCoord2f(0.0,0.0);
			glNormal3fv(vertex_normal[V1]);
			glVertex3fv(vertex[V1]);
		    glColor3fv(colors[1]);
			glTexCoord2f(0.0,1.0);
			glNormal3fv(vertex_normal[V2]);
			glVertex3fv(vertex[V2]);
			glColor3fv(colors[2]);
			glTexCoord2f(1.0,1.0);
			glNormal3fv(vertex_normal[V3]);
			glVertex3fv(vertex[V3]);
			glEnd();
		}
	}
	else //Wireframe
	{
		for(int i = 0; i <= num_traingles; i++)
		{
			//glNormal3fv(tri_normal[i]);
			glBegin(GL_LINE_STRIP);
			glColor3fv(colors[0]);
			V1 = triangle[i][0];
			glColor3fv(colors[1]);
			V2 = triangle[i][1];
			glColor3fv(colors[2]);
			V3 = triangle[i][2];
			glVertex3fv(vertex[V1]);
			glVertex3fv(vertex[V2]);
			glVertex3fv(vertex[V3]);
			glEnd();
		}
	}
}

// Add the normal from each triangle.
static void Accumulate(int triangle_num, int vertex1, int vertex2, int vertex3)
{
	vertex_normal[vertex1][0] += tri_normal[triangle_num][0];
	vertex_normal[vertex1][1] += tri_normal[triangle_num][1];
	vertex_normal[vertex1][2] += tri_normal[triangle_num][2];

	vertex_normal[vertex2][0] += tri_normal[triangle_num][0];
	vertex_normal[vertex2][1] += tri_normal[triangle_num][1];
	vertex_normal[vertex2][2] += tri_normal[triangle_num][2];

	vertex_normal[vertex3][0] += tri_normal[triangle_num][0];
	vertex_normal[vertex3][1] += tri_normal[triangle_num][1];
	vertex_normal[vertex3][2] += tri_normal[triangle_num][2];
}

static void Compute_Triangle_Normals()
{

	int counter_clock_wise = 1;  // May need to invert this for your surface. 

	int	i, vertex1, vertex2, vertex3;
	float	ax, ay, az, bx, by, bz, norm;

    /*tri_normal =  (float **)malloc((num_traingles+1)*sizeof(float  *));
    for(int i = 0; i <= num_traingles ; i++)
	   tri_normal[i] = (float *)malloc(3*sizeof(float));*/

	for (i=0; i<=num_traingles; i++) {
		vertex1 = triangle[i][0];
		vertex2 = triangle[i][1];
		vertex3 = triangle[i][2];

	if (counter_clock_wise) {
			ax = vertex[vertex2][0] - vertex[vertex1][0];
			ay = vertex[vertex2][1] - vertex[vertex1][1];
			az = vertex[vertex2][2]- vertex[vertex1][2];

			bx = vertex[vertex3][0] - vertex[vertex1][0];
			by = vertex[vertex3][1] - vertex[vertex1][1];
			bz = vertex[vertex3][2] - vertex[vertex1][2];
		}
	else {
			ax = vertex[vertex3][0] - vertex[vertex1][0];
			ay = vertex[vertex3][1] - vertex[vertex1][1];
			az = vertex[vertex3][2] - vertex[vertex1][2];

			bx = vertex[vertex2][0] - vertex[vertex1][0];
			by = vertex[vertex2][1] - vertex[vertex1][1];
			bz = vertex[vertex2][2] - vertex[vertex1][2];
		} 

        //  Compute the cross product. 

		tri_normal[i][0] = ay*bz - by*az;
		tri_normal[i][1] = bx*az - ax*bz;
		tri_normal[i][2] = ax*by - bx*ay;

		norm = sqrt(tri_normal[i][0]*tri_normal[i][0] + 
			tri_normal[i][1]*tri_normal[i][1] + 
			tri_normal[i][2]*tri_normal[i][2]);

		if (norm != 0.0) {
			tri_normal[i][0] = tri_normal[i][0]/norm;
			tri_normal[i][1] = tri_normal[i][1]/norm;
			tri_normal[i][2] = tri_normal[i][2]/norm;
		}
	}
}

// Compute the surface normal at each vertex. 
static void Compute_Vertex_Normals()
{
	int     i;
	float	norm;

    /*vertex_normal =  (float **)malloc((num_vertices+1)*sizeof(float  *));
    for(int i = 0; i <= num_vertices ; i++)
	   vertex_normal[i] = (float *)malloc(3*sizeof(float));*/

	for (i=0; i<=num_vertices; i++) {
		vertex_normal[i][0] = 0.0;
		vertex_normal[i][1] = 0.0;
		vertex_normal[i][2] = 0.0;
	}

//  Add the normal at each triangle to the triangle's vertices.
	for (i=0; i<=num_traingles; i++)
		Accumulate(i, triangle[i][0], triangle[i][1], triangle[i][2]);

// Normalize the normal at each vertex. 
	for (i=0; i<=num_vertices; i++) {
		norm = sqrt(vertex_normal[i][0]*vertex_normal[i][0] +
				vertex_normal[i][1]*vertex_normal[i][1] +
				vertex_normal[i][2]*vertex_normal[i][2]);

		if (norm != 0.0) {
			vertex_normal[i][0] /= norm;
			vertex_normal[i][1] /= norm;
			vertex_normal[i][2] /= norm;
		}
	}
}

//Calculate the surface mesh
static void calculateSurface()
{
   float radian;
   float costheta,sintheta;
   float thethaincr = 10.0;
   num_vertices = -1;
   num_traingles = -1;

   /*int Vsize = (numBSplinePts+1)*(360/thethaincr);
   int Tsize = (numBSplinePts+1)*(360/thethaincr)*2;

   vertex =  (float **)malloc(Vsize*sizeof(float  *));
   for(int i = 0; i < Vsize; i++)
		vertex[i] = (float *)malloc(3*sizeof(float));

   triangle =  (int **)malloc(Tsize*sizeof(int  *));
   for(int i = 0; i < Tsize; i++)
		triangle[i] =  (int *)malloc(3*sizeof(int));*/

   for(int i = 0; i <= numBSplinePts; i++)
   {
		vertex[++num_vertices][0] = BSpline[i][0];
		vertex[num_vertices][1] = BSpline[i][1];
		vertex[num_vertices][2] = 0;
   }

   int offset = num_vertices + 1;

   for(float thetha = thethaincr; thetha < 360.0; thetha += thethaincr)
   {
       radian = (thetha*PI)/180;
	   costheta = cos(radian);
	   sintheta = sin(radian);

	   for(int i = 0; i <= numBSplinePts; i++)
	   {
	       vertex[++num_vertices][0] = BSpline[i][0] * costheta;
		   vertex[num_vertices][1] = BSpline[i][1];
		   vertex[num_vertices][2] = -BSpline[i][0] * sintheta;

		   if(i != 0){
		       triangle[++num_traingles][0] = num_vertices;
			   triangle[num_traingles][1] = num_vertices-offset-1;
			   triangle[num_traingles][2] = num_vertices-1;
			   
		       triangle[++num_traingles][0] = num_vertices;
			   triangle[num_traingles][1] = num_vertices-offset;
			   triangle[num_traingles][2] = num_vertices-offset-1;
			   }
        }
	}

	int last_tri_strip = num_traingles+1;
	int prev_vertex = num_vertices-offset+1;
	for(int i = 1; i <= numBSplinePts; i++)
	{
	    prev_vertex++;
		triangle[++num_traingles][0] = i;
		triangle[num_traingles][1] = prev_vertex-1;
		triangle[num_traingles][2] = i-1;

		triangle[++num_traingles][0] = i;
		triangle[num_traingles][1] = prev_vertex;
		triangle[num_traingles][2] = prev_vertex-1;
	}
	Compute_Triangle_Normals();
	Compute_Vertex_Normals();
}

/* Draw the indicated curves using the current control points. */
static void drawCurves(curveType type)
{
    /* Draw the curves */
    calculateKnot();
    calculateBSpline();
	drawSpline();
}

/*Display Control Polygon*/
void displayCtrlPolygon()
{
	glLineWidth(1.0);
	glColor3f(1.0f, 0.0f, 0.0f);
  
	glBegin(GL_LINE_STRIP);
	for (int i = 0; i < ncpts; i++)
	{
      glVertex2fv(cpts[i]);
	}
    glEnd();
}

/*Display Control Points*/
void displayCtrlPoint()
{
	glColor3f(0.0, 0.0, 0.0);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (int i = 0; i < ncpts; i++)
        glVertex3fv(cpts[i]);
    glEnd();
}

//void dealloc_memory()
//{
//	for(int i = 0; i <= num_vertices; i++)
//	{
//		free(vertex_normal[i]);
//	}
//	free(vertex_normal);
//	vertex_normal = NULL;
//	for(int i = 0; i <= num_traingles; i++)
//	{
//		free(triangle[i]);
//	}
//	free(triangle);
//	triangle = NULL;
//	for(int i = 0; i <= num_traingles; i++)
//	{
//		free(tri_normal[i]);
//	}
//	free(tri_normal);
//	tri_normal = NULL;
//	for(int i = 0; i <= num_vertices; i++)
//	{
//		free(vertex[i]);
//	}
//	free(vertex);
//	vertex = NULL;
//}

/* This routine displays the control points */
static void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if(ctrlPoly_On)
		displayCtrlPolygon();
    
	displayCtrlPoint();

	if(BSpline_On)
		drawCurves(BSPLINE);

	if(BSpline_On && BSurface_On)
	{
		calculateSurface();
		drawSurface();
		//dealloc_memory();
	}
    glutSwapBuffers();
}


/* This routine inputs new control points */
static void mouse(int button, int state, int x, int y)
{
    float wx, wy;

    /* We are only interested in left clicks */
    if (button != GLUT_LEFT_BUTTON || state != GLUT_DOWN)
	{
        return;
	}

    /* Translate back to our coordinate system */
    wx = (2.0 * x) / (float)(width - 1) - 1.0;
    wy = (2.0 * (height - 1 - y)) / (float)(height - 1) - 1.0;


	if(ctrlPt_On)
	{
		/* See if we have room for any more control points */
		if (ncpts == MAX_CPTS)
			return;

		/* Save the point */
		cpts[ncpts][0] = wx;
		cpts[ncpts][1] = wy;
		cpts[ncpts][2] = 0.0;
		ncpts++;
	
		/* Draw the point */
		glColor3f(0.0, 0.0, 0.0);
		glPointSize(5.0);
		glBegin(GL_POINTS);
		glVertex3f(wx, wy, 0.0);
		glEnd();
		glutPostRedisplay();
	}

	if(Selection_On && !choice)
	{
		// determine which control point is picked
		float min_dist = FLT_MAX;
		for (int i = 0; i < ncpts; i++)
		{
			// Compute the distance between the user clicked position and control points
			float dist = sqrt((cpts[i][0] - wx)*(cpts[i][0] - wx) 
				+ (cpts[i][1] - wy)*(cpts[i][1] - wy));
			if(dist < min_dist)
			{
				min_dist = dist;
				SelectedCtrlPt = i;
				choice = 1;
				printf("Control point %d was selected.\n",SelectedCtrlPt);
			}
		}

		/* Draw the point */
		glColor3f(0.0, 1.0, 0.0);
		glPointSize(5.0);
		glBegin(GL_POINTS);
		glVertex3f(cpts[SelectedCtrlPt][0], cpts[SelectedCtrlPt][1], 0.0);
		glEnd();
		glutSwapBuffers();
	}

	if(Delete_On)
	{
		int DeletePoint = -1;
		// determine which control point is picked
		float min_dist = FLT_MAX;
		for (int i = 0; i < ncpts; i++)
		{
			float dist = sqrt((cpts[i][0] - wx)*(cpts[i][0] - wx) 
				+ (cpts[i][1] - wy)*(cpts[i][1] - wy));
			if(dist < min_dist)
			{
				min_dist = dist;
				DeletePoint = i;
				choice = 1;
				printf("Control point %d was picked Distance %f.\n",DeletePoint,dist);
			}
		}

		if(DeletePoint != -1)
		{
			for(int i = DeletePoint; i < ncpts-1; i++)
			{
				cpts[i][0] = cpts[i+1][0];
				cpts[i][1] = cpts[i+1][1];
				cpts[i][2] = cpts[i+1][2];
			}
			cpts[ncpts-1][0] = 0.0;
			cpts[ncpts-1][1] = 0.0;
			cpts[ncpts-1][2] = 0.0;
			ncpts--;
			glutPostRedisplay();
		}
	}
	
	if(!Selection_On)
	{
		SelectedCtrlPt = -1;
		choice = 0;
	}
}

// mouse motion function
void myMouseMove(int x, int y)
{
	if ((SelectedCtrlPt > -1) && Selection_On)
	{
	 float wx,wy;
	 wx = (2.0 * x) / (float)(width - 1) - 1.0;
     wy = (2.0 * (height - 1 - y)) / (float)(height - 1) - 1.0;

	 printf("New control point %d location: x = %d   y = %d\n", SelectedCtrlPt, wx, wy);

	 cpts[SelectedCtrlPt][0] = wx;
	 cpts[SelectedCtrlPt][1] = wy;
	 cpts[SelectedCtrlPt][2] = 0.0;

	 //SelectedCtrlPt = -1;
	 //Selection_On = 0;
	 choice = 0;
	 glutPostRedisplay();
	}
}

void userEventAction(int key)
{
   /* menu selects which angle to change or whether to quit */
   switch(key)
   {
		case '0': //Control Point ON/OFF
			ctrlPt_On = !ctrlPt_On;
			Selection_On = 0;
			Delete_On = 0;
			break;

		case '1': //Control Polygon ON/OFF
			ctrlPoly_On = !ctrlPoly_On;
			break;

		case '2': //BSpline Curve ON/OFF
			BSpline_On = !BSpline_On;
			break;

		case '3': //Control Point Selection ON/OFF
			Selection_On = !Selection_On;
			Delete_On = 0;
			ctrlPt_On = 0;
			break;

		case '4': //Delete ON/OFF
			Delete_On = !Delete_On;
			Selection_On = 0;
			ctrlPt_On = 0;
			break;

		case '5': //Save to file
			jFile = fopen(fileName, "w");
			if (jFile == NULL) 
			{
				printf("Warning: Could not open %s\n", fileName);
			}
			else 
			{
				for(int j=0;j < ncpts;j++)
				{
					fprintf(jFile," %f %f", cpts[j][0],cpts[j][1]);
				}
				printf("\nControl Points saved in %s\n", fileName);
			}
			if(jFile != NULL)
			fclose(jFile);
			break;

		case '6': //BSpline Curve ON/OFF
			printf("Loading file %s\n", fileName);
	
			ncpts = 0;
			jFile = fopen(fileName, "r");
			if ( jFile == NULL ) {
				printf("Warning: Could not open %s\n", fileName);		
			}
			else 
			{
		        // Store the control points to event_buffer
				while (!feof(jFile)) {
					fscanf(jFile, "%f %f", &cpts[ncpts][0], &cpts[ncpts][1]);
					//printf("%f %f\n", cpts[ncpts][0], cpts[ncpts][1]);
					ncpts++;
				}
			}
			if(jFile != NULL)
			fclose(jFile);
			break;

		case '7': //BSpline Surface ON/OFF
			BSurface_On = !BSurface_On;
			break;

		case '8': //Clear
			for(int i = 0; i < ncpts; i++)
			{
				cpts[i][0] = 0.0;
				cpts[i][1] = 0.0;
				cpts[i][2] = 0.0;
			}
			ncpts = 0;
			BSurface_On = 0;
			BSpline_On = 0;
			ctrlPoly_On = 0;
			Delete_On = 0;
			Selection_On = 0;
			num_traingles = -1;
			num_vertices = -1;
			numBSplinePts = -1;
			break;

		case '9': //Wireframe ON/OFF
			drawLines = !drawLines;
			break;

		case 'f': //Smooth Shading ON/OFF
			smoothEnabled = !smoothEnabled;
			break;

		case 't': //Texture ON/OFF
			texEnabled = !texEnabled;
			break;

		case 'l': //Lightning ON/OFF
			lightingEnabled = !lightingEnabled;
			break;

		case 'q': //Exit
			/*for(int i = 0; i <= num_vertices; i++)
			{
				free(vertex_normal[i]);
			}
			free(vertex_normal);
			for(int i = 0; i <= num_traingles; i++)
			{
				free(triangle[i]);
			}
			free(triangle);
			for(int i = 0; i <= num_traingles; i++)
			{
				free(tri_normal[i]);
			}
			free(tri_normal);
			for(int i = 0; i <= num_vertices; i++)
			{
				free(vertex[i]);
			}
			free(vertex);*/
			exit(0);
			break;

		default:break;
   }
   userSettings();
   glutPostRedisplay();
}

/* This routine handles keystroke commands */
void keyboard(unsigned char key, int x, int y)
{
   userEventAction(key);
}

/* This routine handles window resizes */
void reshape(int w, int h)
{
    width = w;
    height = h;

    /* Set the transformations */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glViewport(0, 0, w, h);
}


/*----------------------------------------------------------------------*/

typedef struct menuEntryStruct {
    char *label;
    char key;
} menuEntryStruct;

static menuEntryStruct mainMenu[] = {
    "Ctrl Pt ON/OFF", 			'0',
    "Ctrl Polygon ON/OFF", 		'1',
    "BSpline ON/OFF", 			'2',
	"Surface ON/OFF", 			'7',	
	"lines/polygons",			'9',
	"Flat/Smooth",				'f',
	"Texture ON/OFF",			't',
	"Lightning ON/OFF",			'l',
    "Select ON/OFF", 			'3',
    "Delete ON/OFF", 			'4',
    "Save to File", 			'5',
    "Load from File", 			'6',
    "Clear", 					'8',
    "quit", 					'q',
	
};
int mainMenuEntries = sizeof(mainMenu)/sizeof(menuEntryStruct);

void selectMain(int choice)
{
    userEventAction(mainMenu[choice].key);
}

void setMenuEntries(bool init)
{
    int i;

    if (init) {
	glutCreateMenu(selectMain);
	for (i=0; i < mainMenuEntries; i++) {
	    glutAddMenuEntry(mainMenu[i].label, i);
	}
	glutAttachMenu(GLUT_RIGHT_BUTTON);
    } else {
    }
}

/*----------Init---------------------*/
void init(void)
{
   glClearColor (0.0, 0.0, 0.0, 0.0);
   glEnable(GL_DEPTH_TEST);
   setMenuEntries(true);
   glColor3f (0.0, 0.0, 0.0);
  // glShadeModel(GL_FLAT);
}

void main(int argc, char **argv)
{
    /* Intialize the program */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(width, height);
    glutCreateWindow("curves");
	initTexture();
	init();
	userSettings();

    /* Register the callbacks */
    glutDisplayFunc(display);
    glutMouseFunc(mouse);
    glutKeyboardFunc(keyboard);
	glutMotionFunc(myMouseMove);
    glutReshapeFunc(reshape);
	glClearColor(1.0, 1.0, 1.0, 1.0);

    glutMainLoop();
}