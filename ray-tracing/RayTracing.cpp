/*******************************************************************************

  Purpose:	This program renders three spheres and three planes using 
			non-recursive ray tracing.

			The purpose of this program is to demonstrate a basic ray tracer. 

*******************************************************************************/

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gl/glut.h>

#define DEBUG   0

/*  Define some constants.  */
#define	PI	3.1415926

/* Max image size allowed. */

#define MAX_SIZE 512

/**/
#define MAX_DEPTH 1
/*  Define some structures.  */

#define NUM_SPHERES 3
#define NUM_PLANES 3

struct	points	{
		   float   x, y, z;
};

typedef struct	rgb_struct	{
		   float   r, g, b;
} rgb;

/*  Viewing parameters.  */

struct	points	from, at, up;
float	VXR, VXL, VYB, VYT;
float	ax, ay, az, bx, by, bz, cx, cy, cz; /* line of sight and others for matrix*/
float	viewangle, angle, tanv2;
float	xinterval, yinterval;

/*  Illumination parameters.  */

float	lsx, lsy, lsz;
rgb	il, ia;
rgb ka_S[NUM_SPHERES];
rgb ka_P[NUM_PLANES];

rgb kd_S[NUM_SPHERES];
rgb kd_P[NUM_PLANES];

rgb ks_S[NUM_SPHERES];
rgb ks_P[NUM_PLANES];

float ns_S[NUM_SPHERES];
float ns_P[NUM_PLANES];

float nt_S[NUM_SPHERES];
float nt_P[NUM_PLANES];

float krg_S[NUM_SPHERES];
float ktg_S[NUM_SPHERES];

float krg_P[NUM_PLANES];
float ktg_P[NUM_PLANES];

rgb	tka1, tkd1, tka2, tkd2, tka3, tkd3;
rgb	tka4, tkd4, tka5, tkd5, tka6, tkd6;

/*  Image parameters.  */

int		xmax_pixel, ymax_pixel;

/* Image buffer.  A more efficient approach is to use one single array (texture_RGB)
   rather than using the three rays.  */

float *texture_R;
float *texture_G;
float *texture_B;

/*  Object parameters.  */

float	c[NUM_SPHERES][4]; //Sphere coordinates
float	p[NUM_PLANES][4];  //Plane coordinates
float	pminmax[NUM_PLANES][3][2];  //Plane coordinates

/* Image output file. */
FILE	*outpfile;

/*******************************************************************************

  Title:	Read_Information

  Purpose:	This function reads in the information about the objects (three spheres
            and three planes) to be rendered.  The information is assumed to be 
			stored in an ASCII text file called "ray.dat".	

Object Parameters: [ 3 Spheres and 3 Planes]
		Sphere l m n r ns nt krg ktg ka kd ks
		Planes a b c d ns nt krg ktg ka kd ks
			   xmin ymin zmin xmax ymax zmax
*******************************************************************************/

int Read_Information()
{
	char	str[132];
	FILE	*inpfile;

	if ( (inpfile = fopen("ray.dat","r"))==NULL) {
		printf("ERROR: Could not open ray.dat for read!\n");
		return(0);
	}

    /*Read in viewing information.*/
    fscanf(inpfile,"%s %f %f %f\n",str, &from.x, &from.y, &from.z);
	fscanf(inpfile,"%s %f %f %f\n",str, &at.x, &at.y, &at.z);
	fscanf(inpfile,"%s %f %f %f\n",str, &up.x, &up.y, &up.z);
	fscanf(inpfile,"%s %f\n",str, &viewangle);
	angle = viewangle * PI/180.0;
	tanv2 = tan(angle/2.0);
	printf("tanv2 = %f\n", tanv2);

    /*Read in the viewport information */
    fscanf(inpfile,"%s %f %f\n", str, &VXL, &VXR);
	fscanf(inpfile,"%s %f %f\n", str, &VYB, &VYT);

    /*Read in the light vector (L).  Note: To specify a light source position, the light 
    vector will need to be computed in Compute_Color().*/
    fscanf(inpfile,"%s %f %f %f\n",str, &lsx, &lsy, &lsz);

	printf("From: %f %f %f\n", from.x, from.y, from.z);
	printf("At: %f %f %f\n", at.x, at.y, at.z);
	printf("Up: %f %f %f  View Angle: %f \n", up.x, up.y, up.z, viewangle);

    /*Read in spheres' information.  */
	for(int i = 0; i < NUM_SPHERES; i++)
	{
		fscanf(inpfile,"%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
			str, &c[i][0], &c[i][1], &c[i][2], &c[i][3],
			&ns_S[i], &nt_S[i], &krg_S[i], &ktg_S[i],
			&ka_S[i].r, &ka_S[i].g, &ka_S[i].b,
			&kd_S[i].r, &kd_S[i].g, &kd_S[i].b,
			&ks_S[i].r, &ks_S[i].g, &ks_S[i].b);

		/*printf("%s l = %f m = %f n = %f r = %f \n ns = %f  nt = %f krg = %f ktg = %f \n ka = (%f %f %f) \n kd = (%f %f %f) \n kd = (%f %f %f)\n",
			str, c[i][0], c[i][1], c[i][2], c[i][3], ns_S[i], nt_S[i], krg_S[i], ktg_S[i],
			ka_S[i].r,ka_S[i].g,ka_S[i].b,
			kd_S[i].r,kd_S[i].g,kd_S[i].b,
			ks_S[i].r,ks_S[i].g,ks_S[i].b);
		printf("\n");*/
	}

    /*Read in checker boards' (planes') information.  */
	for(int i = 0; i < NUM_SPHERES; i++)
	{
		fscanf(inpfile,"%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",str, &p[i][0], &p[i][1], &p[i][2], &p[i][3], 
			   &ns_P[i], &nt_P[i], &krg_P[i], &ktg_P[i],
		       &ka_P[i].r, &ka_P[i].g, &ka_P[i].b,
		       &kd_P[i].r, &kd_P[i].g, &kd_P[i].b,
		       &ks_P[i].r, &ks_P[i].g, &ks_P[i].b);

		/*printf("%s a = %f  b= %f c = %f d = %f \n ns = %f nt = %f krg = %f ktg = %f\n ka = (%f %f %f)\n kd = (%f %f %f) \n ks = (%f %f %f)\n",str, p[0][0], p[0][1], p[0][2], p[0][3],
			   ns_P[i], nt_P[i], krg_P[i], ktg_P[i],
		       ka_P[i].r, ka_P[i].g, ka_P[i].b,
		       kd_P[i].r, kd_P[i].g, kd_P[i].b,
		       ks_P[i].r, ks_P[i].g, ks_P[i].b);
		printf("\n");*/
	}

	fscanf(inpfile,"%s %f %f %f %f %f %f\n", str, &pminmax[0][0][0], &pminmax[0][1][0],
		&pminmax[0][2][0], &pminmax[0][0][1], &pminmax[0][1][1],
		&pminmax[0][2][1]);
	fscanf(inpfile,"%s %f %f %f %f %f %f\n", str, &pminmax[1][0][0], &pminmax[1][1][0],
		&pminmax[1][2][0], &pminmax[1][0][1], &pminmax[1][1][1],
		&pminmax[1][2][1]);
	fscanf(inpfile,"%s %f %f %f %f %f %f\n", str, &pminmax[2][0][0], &pminmax[2][1][0],
		&pminmax[2][2][0], &pminmax[2][0][1], &pminmax[2][1][1],
		&pminmax[2][2][1]);

    /* Read in the image size.  */
    fscanf(inpfile,"%s %d %d\n",str, &xmax_pixel, &ymax_pixel);

	/* Make sure the image size does not exceed MAX_SIZE x MAX_SIZE.  */
	if (xmax_pixel > MAX_SIZE || ymax_pixel > MAX_SIZE) {
		printf("Error: Exceeded max image size %d x %d\n", xmax_pixel, ymax_pixel);
		printf("Reset to max image size: %d x %d\n", MAX_SIZE, MAX_SIZE);
		xmax_pixel = MAX_SIZE - 1;
		ymax_pixel = MAX_SIZE - 1;
	}

	fclose(inpfile);

    /*Open an output file to store the intensity values of the output image.  */
	if ((outpfile = fopen("image.out","wb")) == NULL) {
		printf("ERROR:  cannot open image.out for write.\n");
		return(0);
	}

    /* Allocate memory for the image buffer.  */
    texture_R = new float [xmax_pixel * ymax_pixel ];
	texture_G = new float [xmax_pixel * ymax_pixel ];
	texture_B = new float [xmax_pixel * ymax_pixel ];
	printf("image_buf allocated.  Image size %d x %d\n", xmax_pixel, ymax_pixel);

	return(1);
}


/*******************************************************************************

  Title:	Normalize

  Purpose:	This function normalizes the given vector.

*******************************************************************************/

void Normalize(float *x,float *y,float *z)
{
	float	norm;

	norm = sqrt( *x * *x + *y * *y + *z * *z );
	if (norm != 0.0) {
		*x = *x / norm;
		*y = *y / norm;
		*z = *z / norm;
	}
}

/*******************************************************************************

  Title:	Power

  Purpose:	This function computes the power of the given base and
		    exponent.

*******************************************************************************/

float 	Power(float base,int exp)
{
	int	i;
	float	value;

	value = 1.0;
	for (i=1; i<=exp; i++)
		value *= base;

	return( value );
}


/*******************************************************************************

  Title:	Compute_M

  Purpose:	This function computes the transformation matrix to be used
		    in the perspective viewing model.

*******************************************************************************/

void Compute_M()
{

/*  Compute the line-of-sight vector, c.  */

	cx = at.x - from.x;
	cy = at.y - from.y;
	cz = at.z - from.z;
	Normalize(&cx, &cy, &cz);

/*  Compute the cross product of vector c and the up vector.  */

	ax = cy*up.z - up.y*cz;
	ay = up.x*cz - cx*up.z;
	az = cx*up.y - up.x*cy;
	Normalize(&ax, &ay, &az);

/*  Compute the cross product of vector a and c.  */

	bx = ay*cz - cy*az;
	by = cx*az - ax*cz;
	bz = ax*cy - cx*ay;
}

/*******************************************************************************

  Title:	Setup_Parameters

  Purpose:	This function sets up the necessary parameters for 
		    performing the ray trace.  It first computes the 
		    transformation matrix for the perspective viewing model, then
		    sets up the default illumination parameters.

*******************************************************************************/

void Setup_Parameters()
{

	/*  Compute the transformation matrix for converting world coordinates to eye
    coordinates.  */

	Compute_M(); 

	/*  Normalized the given directional light vector. 

	Note:  DO NOT normalize this vector, if the position of the light source is given.
	       The light vector (L) would need to computed for each intersction point,
		   then normalized in Compute_Color().   */

	Normalize(&lsx, &lsy, &lsz);
	printf("light %f %f %f\n", lsx, lsy, lsz);

    /*  Set up the conversion factors for converting from pixel coordinates to 
    view port coordinates.  */

	xinterval = (VXR - VXL) / xmax_pixel;
	yinterval = (VYT - VYB) / ymax_pixel;

    /*  Set up default illumination (Phong lighting) parameters.  */

	il.r = 1.0;	il.g = 1.0;	il.b = 1.0;
	ia.r = 1.0;	ia.g = 1.0;	ia.b = 1.0;

    /* Phong lighting parameters for the three spheres.  */

	tka1.r = 0.2;	tka1.g = 0.0;	tka1.b = 0.0;
	tkd1.r = 0.2;	tkd1.g = 0.0;	tkd1.b = 0.0;

	tka2.r = 0.0;	tka2.g = 0.2;	tka2.b = 0.0;
	tkd2.r = 0.0;	tkd2.g = 0.2;	tkd2.b = 0.0;

	tka3.r = 0.0;	tka3.g = 0.0;	tka3.b = 0.2;
	tkd3.r = 0.0;	tkd3.g = 0.0;	tkd3.b = 0.2;

    /* Phong lighting parameters for the three planes (not shown).  */

	tka4.r = 0.4;	tka4.g = 0.0;	tka4.b = 0.0;
	tkd4.r = 0.7;	tkd4.g = 0.0;	tkd4.b = 0.0;

	tka5.r = 0.0;	tka5.g = 0.0;	tka5.b = 0.8;
	tkd5.r = 0.0;	tkd5.g = 0.0;	tkd5.b = 0.7;

	tka6.r = 0.0;	tka6.g = 0.1;	tka6.b = 0.0;
	tkd6.r = 0.0;	tkd6.g = 0.7;	tkd6.b = 0.0;
}

/*******************************************************************************

  Title:	Bump Map Textures
  Purpose:	functions below computes the bump map texture for each point

*******************************************************************************/

int TABLE_SIZE_1 = 60;
int TABLE_SIZE = 60;
float noise_tabl[65][65];

float calc_noise(float iu,float iv, int direction)
{
	int i,j,x,y,left,right,val;
	float noise,u,v,w,ul,vl,wl;

    //Filling noise table with random noise
    for (i = 0; i < 65; i++)
	{
        for (j = 0; j < 65; j++)
		{
            val = rand() % 255;
            noise_tabl[i][j] = val / 256;
		}
	}

	i = (int) iu;
	j = (int) iv;
	x = i % TABLE_SIZE;
	y = j % TABLE_SIZE;

	if(direction == 1)
	{
		if(x <= 0)
			left = 0;
		else
			left = x-1;
		if(x >= TABLE_SIZE_1)
			right = TABLE_SIZE_1;
		else
			right = x+1;
		noise = (noise_tabl[right][y] - noise_tabl[left][y])/2.0;
	}
	else
	{
		if(y <= 0)
			left = 0;
		else
			left = y -1;
		if(y >= TABLE_SIZE_1)
			right = TABLE_SIZE_1;
		else
			right = y+1;
		noise = (noise_tabl[x][right] - noise_tabl[x][left])/2.0;
	}
	return noise;
}

void Bump_map(float x, float y, float z, float xc, float yc, float zc, float r, float *nx, float *ny , float *nz, int noise_method)
{
	float xp,yp,zp,iu,iv,xu,yu,zu,xv,yv,zv;
	float fu,fv,a,dx,dy,dz,u,v,nnx,nny,nnz;

	/*Translate to origin*/
	xp = (x - xc)/r;
	yp = (y - yc)/r;
	zp = (z - zc)/r;

	/*convert to (u,v) coordinates*/
	u = asin(zp);
	v = atan2(yp, xp);

	/*convert to integer (iu,iv) coordinates*/
	iu = (u + PI) / (2 * PI) * (TABLE_SIZE_1);
	iv = (u + 0.5*PI) / PI * (TABLE_SIZE_1);

	/*get the partials*/
	xu = -r * cos(v) * sin(v);
	xv = -r * sin(v) * cos(u);
	yu =  r * cos(v) * cos(u);
	yv = -r * sin(v) * sin(u);
	zu = 0.0;
	zv = r * cos(v);

	/*calculate the pertubations*/
	fu = calc_noise(iu,iv,1);
	fv = calc_noise(iu,iv,2);

	/*compute D*/
	if(noise_method == 1)
	{
		dx = fu*xu + fv*xv;
		dy = fu*yu + fv*yv;
		dz = fu*zu + fv*zv;
	}
	else
	{
		nnx = *nx;
		nny = *ny;
		nnz = *nz;
		Normalize(&nnx, &nny, &nnz);

		/*Compute the cross product of Pu * N */

		dx = fv*(yu*nnz - nny*zu);
		dy = fv*(nnx*zu - xu*nnz);
	    dz = fv*(xu*nny - nnx*yu);

		/*Compute the cross product of Pv * N */

		dx += fu*(yv*nnz - nny*zv);
		dy += fu*(nnx*zv - xv*nnz);
	    dz += fu*(xv*nny - nnx*yv);
	}

	/*Normalize and Scale D*/
	Normalize(&dx,&dy,&dz);

	a = sqrt(fu*fu + fv*fv);
	dx *= a;
	dy *= a;
	dz *= a;
	*nx += dx;
	*ny += dy;
	*nz += dz;
}

/*******************************************************************************

  Title:	Wood Grain Textures

  Purpose:	functions below computes the wood grain texture for each point

*******************************************************************************/

int round(float number)
{
    return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

int wood_grain(float u,float v,float w)
{
	float radius,angle;
	int grain;

	radius = sqrt(u*u + w*w);
	if ( w == 0 ) 
		angle = PI/2;
	else
		angle = atan(u/w);
	radius = radius + 2*sin(20*angle+v/4);
	grain = round(radius) % 6;
	return grain;
}

/*******************************************************************************

  Title:	Marble Textures

  Purpose:	functions below computes the marble texture for each point

*******************************************************************************/

void get_marble_color(float value, rgb *color)
{
	float x;
    x = sqrt( value + 1.0)*0.7071;
	color->g = 0.3 + 0.4*x;
	x = sqrt( x );
	color->r = 0.3 + 0.2*x;
	color->b = 0.6 + 0.7*x;
}

/* Evaluate the turbulence function based on the given point.  */

float turbulence(float x,float y,float z,float pixel_size)
{

	float t, scale, new_point[3];

	t = 0;
	new_point[0] = x;
	new_point[1] = y;
	new_point[2] = z;

	for (scale=1.0; scale>pixel_size; scale/=2.0) {
		new_point[0] /= scale;
		new_point[1] /= scale;
		new_point[2] /= scale;
		t += calc_noise(new_point[0],new_point[1],2)*scale;
	}
	return( t );
}


/*******************************************************************************

  Title:	Check_Sphere

  Purpose:	This function determines if the given ray intercepts the given
		    sphere.

*******************************************************************************/

void Check_Sphere(float px,float py,float pz,float dx,float dy,float dz,float xc,float yc,float zc,float r,float *t1, float *t2)
{
	float	a, b, c, xdiff, ydiff, zdiff, discr;

	xdiff = px-xc;
	ydiff = py-yc;
	zdiff = pz-zc;
	a = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff - r*r;
	b = 2.0*( dx*xdiff + dy*ydiff + dz*zdiff );
	c = dx*dx + dy*dy + dz*dz;

/*  Check if there are any intersections.  */

	discr = b*b - 4.0*a*c;
	if (discr < 0.0) {
		*t1 = -1.0;
		*t2 = -1.0;
	}
	else if (discr == 0.0) {
		*t1 = -b / (2.0*c);
		*t2 = -1.0;
	}
	else {
		discr = sqrt(discr);
		*t1 = (-b + discr) / (2.0*c);
		*t2 = (-b - discr) / (2.0*c);
	}
}

/*******************************************************************************

  Title:	Check_Plane

  Purpose:	This function checks if the given ray intercepts the given
		    plane.

*******************************************************************************/

void Check_Plane(float px,float py,float pz,float dx,float dy,float dz,float a,float b,float c,float d,float *t1)
{
	*t1 = (-a*px - b*py - c*pz - d) / (a*dx + b*dy + c*dz);
}


/*******************************************************************************

  Title:	Compute_Intersection

  Purpose:	This function computes the intersection of ray with an
		    object.  The intersection point is given by a parametric value
		    t, where ray = p + d*t, d = the direction of the ray, and p is
		    the starting point of the ray.

*******************************************************************************/

void Compute_Intersection(float px,float py,float pz,float dx,float dy, float dz,float t,float *newx,float *newy,float *newz)
{
	*newx = px + t*dx;
	*newy = py + t*dy;
	*newz = pz + t*dz;
}

/*******************************************************************************

  Title:	Compute_Color

  Purpose:	This function computes the intensity of the color for the
		    given location based on the Phong lighting model.

*******************************************************************************/

rgb Compute_Color(int shadow_flag, float ipx,float ipy,float  ipz,float  nx,float  ny,float  nz,
				   rgb ia,rgb ka,rgb kd, rgb ks,int n,int krg, int ktg, int depth)
{
	float	vx, vy, vz, rx, ry, rz;
	float	ndotl, vdotr, cosalphapower;
	rgb RGB_color;

/*  Compute the view vector.  */

	vx = from.x - ipx;
	vy = from.y - ipy;
	vz = from.z - ipz;
	Normalize(&vx, &vy, &vz);

/*  Compute the R (reflection) vector.  */

	ndotl = nx*lsx + ny*lsy + nz*lsz;
	rx = 2.0*ndotl*nx - lsx;
	ry = 2.0*ndotl*ny - lsy;
	rz = 2.0*ndotl*nz - lsz;

	/* Compute the V (view) vector. */

	vdotr = vx*rx + vy*ry + vz*rz;

	/* Compute Ia * Ka.  */

	RGB_color.r = ia.r * ka.r;
	RGB_color.g = ia.g * ka.g;
	RGB_color.b = ia.b * ka.b;

	/* Compute diffuse reflection. */

	if (ndotl >= 0.0 && shadow_flag==0) { 

		/*  diffuse reflection = kd * N dot L * Il  */

		RGB_color.r = RGB_color.r + kd.r*ndotl*il.r;
		RGB_color.g = RGB_color.g + kd.g*ndotl*il.g;
		RGB_color.b = RGB_color.b + kd.b*ndotl*il.b;

		if (vdotr >= 0.0) {

			/*  specular reflection = ks * cos(alpha)**K^n * Il */

			cosalphapower = Power(vdotr, n);

			RGB_color.r = RGB_color.r + ks.r*cosalphapower*il.r;
			RGB_color.g = RGB_color.g + ks.g*cosalphapower*il.g;
			RGB_color.b = RGB_color.b + ks.b*cosalphapower*il.b;
		}
	}
	/*  Make sure that the color is within range.  */

	if (RGB_color.r > 1.0) RGB_color.r = 1.0;
	if (RGB_color.g > 1.0) RGB_color.g = 1.0;
	if (RGB_color.b > 1.0) RGB_color.b = 1.0;

	return RGB_color;
}

/*******************************************************************************
  Title:	compute_refracted_ray

  Purpose:	This function generates refracted ray of a given ray.

*******************************************************************************/

void compute_refracted_ray(float dx, float dy, float dz, float nx, float ny, float nz,
    float *rfx, float *rfy, float *rfz, float nt){
    
	/*[T = ηt (N · I) – √1 – ηt*ηt (1 – (N · I)*(N . I) ) ] N – ηt I */
	
	float ndoti;
	ndoti = (nx*dx)+(ny*dy)+(nz*dz);
	*rfx = nt*ndoti - sqrt(1 - (nt*nt*(1 - ndoti*ndoti)))*nx - nt*dx;
	*rfy = nt*ndoti - sqrt(1 - (nt*nt*(1 - ndoti*ndoti)))*ny - nt*dy;
	*rfz = nt*ndoti - sqrt(1 - (nt*nt*(1 - ndoti*ndoti)))*nz - nt*dz;
}


/*******************************************************************************
  Title:	compute_reflected_ray

  Purpose:	This function generates reflected ray of a given ray.

*******************************************************************************/

void compute_reflected_ray(float dx, float dy, float dz, float nx, float ny, float nz,
    float *rlx, float *rly, float *rlz){
    /*
    c1 = -dot_product( N, V )
    Rl = V + (2 * N * c1 )
    */
    float c1 = -(nx*dx + ny*dy + nz*dz);
    *rlx = dx + ( 2 * nx * c1);
    *rly = dy + ( 2 * ny * c1);
    *rlz = dz + ( 2 * nz * c1);

}

/*******************************************************************************
  Title:	Check_Shadow

  Purpose:	This function checks if the point is in shadow.

*******************************************************************************/

int Check_Shadow(float ipx, float ipy, float ipz, int obj_num)
{
	float t1, t2;
	for(int i = 1; i <= NUM_SPHERES; i++)
	{
		Check_Sphere(ipx, ipy, ipz, lsx, lsy, lsz, c[i-1][0], 
				c[i-1][1], c[i-1][2], c[i-1][3], &t1, &t2);

		if ((t1>=0.00001 || t2>=0.00001) && i != obj_num){
			return 1;
		}

	}

	for(int i = 1; i <= NUM_PLANES; i++)
	{
		Check_Plane(ipx, ipy, ipz, lsx, lsy, lsz, p[i-1][0], p[i-1][1], p[i-1][2], p[i-1][3], &t1);

		if (t1 >= 0.00001 && i != (NUM_SPHERES+1+obj_num))
			return 1;
	}
	return 0;
}

/*******************************************************************************

  Title:	Ray_Trace	

  Purpose:	This function performs simple ray tracing.  This is a non-recursive
            ray tracing without any reflection and refraction rays. 

*******************************************************************************/

rgb Ray_Trace(int obj_num,float dx, float dy,float dz, points from, int depth)
{
	int	    shadow_flag;
	int	    texture, buf_ptr  = 0;
	int     num_image = 1;
	float	nx, ny, nz;
	float	t_min, t1, t2, ipx = 0, ipy = 0, ipz = 0;
	float   rlx, rly, rlz, rfx,rfy,rfz;
	rgb		rgb_color;

	/*  Check if the current ray intercepts spheres */
	t_min = 999.0;
	//obj_num = 0;
	texture = 0;

	for(int i = 1; i <= NUM_SPHERES; i++)
	{
		Check_Sphere(from.x, from.y, from.z, dx, dy, dz, c[i-1][0], 
				c[i-1][1], c[i-1][2], c[i-1][3], &t1, &t2);

		if (t1>=0.00001) {
			t_min = t1;
			obj_num = i;
			Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t1, &ipx, &ipy, &ipz);
		}

		if (t2>=0.00001 && t2<t_min) {
			t_min = t2;
			obj_num = i;
			Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t2, &ipx, &ipy, &ipz);
		}

	}
	
	/*  Check if the current ray intercepts planes.  */
	for(int i = 1; i <= NUM_PLANES; i++)
	{
		Check_Plane(from.x, from.y, from.z, dx, dy, dz, p[i-1][0], p[i-1][1], p[i-1][2], p[i-1][3], &t1);

		if (t1 >= 0.00001 && t1<t_min) {
			/*  Check if the intersection point is inside the min/max values. */

			Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t1, &ipx, &ipy, &ipz);

			if (ipx >= pminmax[i-1][0][0] && ipx <= pminmax[i-1][0][1] && 
				ipy >= pminmax[i-1][1][0]  && ipy <= pminmax[i-1][1][1] && 
				ipz >=  pminmax[i-1][2][0] && ipz <= pminmax[i-1][2][1] ) { 

					t_min = t1;
					obj_num = i+NUM_SPHERES;
			}
		}
	}

	int grain;
	float u,v,w,t,value,x;

	/* Compute the intensity to use at the current pixel.  */
    switch (obj_num) {

	/*  The current ray does not intersect any of the objects.  */
    case 0 :rgb_color.r = 0.1;
			rgb_color.g = 0.1;
			rgb_color.b = 0.1;
			break;

	/*  The current ray intercept sphere #1.  */
	case 1 : 
			nx = ipx - c[0][0];
			ny = ipy - c[0][1];
			nz = ipz - c[0][2];
			Normalize(&nx, &ny, &nz);
					 
			shadow_flag = 0;
			shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
			texture = 0;

			
			rgb_color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka_S[0], kd_S[0], ks_S[0], ns_S[0], krg_S[0], ktg_S[0], 1);
			
		    break;

	/*  The current ray intercepts sphere #2.  */
	case 2 : 
			nx = ipx - c[1][0];
			ny = ipy - c[1][1];
			nz = ipz - c[1][2];
			Normalize(&nx, &ny, &nz);
			shadow_flag = 0;
			shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
					
		    texture = rand() % 2;
			
			if(texture == 1)
			{
				Bump_map(ipx,ipy,ipz,dx,dy,dz,1,&nx,&ny,&nz,2);
				rgb_color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka_S[1], kd_S[1], ks_S[1], ns_S[1], krg_S[1], ktg_S[1], 1);
			}
			else
				rgb_color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka_S[1], kd_S[1], ks_S[1], ns_S[1], krg_S[1], ktg_S[1], 1);
		    break;

	/*  The current ray intercepts sphere #3.  */
	case 3 : 
			nx = ipx - c[2][0];
			ny = ipy - c[2][1];
			nz = ipz - c[2][2];
			Normalize(&nx, &ny, &nz);
			shadow_flag = 0;
			shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);	 
			texture = rand() % 2;

		    if (texture==1) // Compute texture. */
			{
				Bump_map(ipx,ipy,ipz,dx,dy,dz,1,&nx,&ny,&nz,1);
				rgb_color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka3, tkd3, ks_S[2], ns_S[2], krg_S[2], ktg_S[2], 1);
			}
			else
				rgb_color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka_S[2], kd_S[2], ks_S[2], ns_S[2], krg_S[2], ktg_S[2], 1);
			break;
    
	case 4 :  /*  The current ray intercepts plane #1.  */
			nx = p[0][0];
		    ny = p[0][1];
			nz = p[0][2];
			shadow_flag = 0;
			shadow_flag = Check_Shadow(ipx, ipy, ipz,obj_num);	

			grain = wood_grain(ipx,ipy,ipz);

			if(grain < 3)
				rgb_color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka4 ,tkd4 , ks_P[0], ns_P[0] ,krg_P[0], ktg_P[0], 1);
			else
				rgb_color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka_P[0], kd_P[0], ks_P[0], ns_P[0], krg_P[0], ktg_P[0], 1);
			break;
			 
   case 5 :  /*  The current ray intercepts plane #2.  */
			nx = p[1][0];
			ny = p[1][1];
			nz = p[1][2];
			shadow_flag = 0;
			shadow_flag = Check_Shadow(ipx, ipy, ipz,obj_num);	
					
			if (ipz < 2.0 || (ipz>=4.0 && ipz<6.0)) {
				if ((ipy>=2.0 && ipy<4.0) || (ipy>=6.0))
					texture = 1;
				else
					texture = 0;
			}
			else {
				if ((ipy<2.0) || (ipy>=4.0 && ipy<6.0)) 
					texture = 1;
				else
					texture = 0;
			}
			if (texture == 1)
				rgb_color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka5, tkd5, ks_P[1], ns_P[1], krg_P[1], ktg_P[1], 1);
			else
				rgb_color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka_P[1], kd_P[1], ks_P[1], ns_P[1], krg_P[1], ktg_P[1], 1);
			break;

	case 6 : /*  The current ray intercepts plane #3.  */
			nx = p[2][0];
			ny = p[2][1];
		    nz = p[2][2];

			shadow_flag = 0;
		    shadow_flag = Check_Shadow(ipx, ipy, ipz,obj_num);	
	
			if (ipx < 2.0 || (ipx >= 4.0 && ipx < 6.0)) {
			     if ((ipz >= 2.0 && ipz < 4.0) || (ipz >=6.0))
					texture = 1;
				else
					texture = 0;
			}
			else {
				 if ((ipz<2.0) || (ipz>=4.0 && ipz<6.0)) 
					texture = 1;
				 else
					texture = 0;
			}
			t = turbulence(ipx,ipy,ipz,0.5);
			value = ipx+t;
			x = sqrt( value + 1.0)*0.7071;
			get_marble_color( (float)sin(value*3.1415926), &rgb_color);
			get_marble_color( (float)sin(value*3.1415926), &rgb_color);
			get_marble_color( (float)sin(value*3.1415926), &rgb_color);
			
			/*if (texture == 1)
				rgb_color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka6, tkd6, ks_P[2], ns_P[2], krg_P[2], ktg_P[2], 1);
			else
				rgb_color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka_P[2], kd_P[2], ks_P[2], ns_P[2], krg_P[2], ktg_P[2], 1);
			*/
			break;

		default:
			break;
	}

	if(depth < MAX_DEPTH)
	{
		if(obj_num > 0 && obj_num !=4 && nx >= 0 && ny >= 0 && nz >= 0)
		{
			compute_reflected_ray(dx, dy, dz, nx, ny, nz, &rlx, &rly, &rlz);

			points newFrom;
			newFrom.x = ipx;
			newFrom.y = ipy;
			newFrom.z = ipz;
			rgb rcolor = Ray_Trace(obj_num,rlx,rly,rlz,newFrom,depth+1);
			switch(obj_num)
			{
				case 1:
					rgb_color.r += krg_S[0]*rcolor.r;
					rgb_color.g += krg_S[0]*rcolor.g;
					rgb_color.b += krg_S[0]*rcolor.b;
					break;

				case 2:
					rgb_color.r += krg_S[1]*rcolor.r;
					rgb_color.g += krg_S[1]*rcolor.g;
					rgb_color.b += krg_S[1]*rcolor.b;
					break;

				case 3:
					rgb_color.r += krg_S[2]*rcolor.r;
					rgb_color.g += krg_S[2]*rcolor.g;
					rgb_color.b += krg_S[2]*rcolor.b;
					break;
			  
				case 4:
					rgb_color.r += krg_P[0]*rcolor.r;
					rgb_color.g += krg_P[0]*rcolor.g;
					rgb_color.b += krg_P[0]*rcolor.b;
					break;

				case 5:
					rgb_color.r += krg_P[1]*rcolor.r;
					rgb_color.g += krg_P[1]*rcolor.g;
					rgb_color.b += krg_P[1]*rcolor.b;
					break;
			
				case 6:
					rgb_color.r += krg_P[2]*rcolor.r;
					rgb_color.g += krg_P[2]*rcolor.g;
					rgb_color.b += krg_P[2]*rcolor.b;
					break;

				default:
					break;
			}
	    }

		float n1 = 1.0003; /*refractive index air*/
		/*nt -  refractive index material*/
		float thetha = 1000;

		if(obj_num > 0 && obj_num < 7)
		{
			if(obj_num < 4)
			{
				thetha = asin((float)n1/(float)nt_S[obj_num-1]);
				
				if(thetha < 41.2) /* Check if object is semi transparent */
				{
					points newFrom;
					newFrom.x = ipx;
					newFrom.y = ipy;
					newFrom.z = ipz;
			
					compute_refracted_ray(dx,dy,dz,nx,ny,nz,&rfx,&rfy,&rfz,nt_S[obj_num-1]); /*Compute refracted ray*/
					rgb rcolor = Ray_Trace(obj_num,rfx,rfy,rfz,newFrom,depth+1);
					rgb_color.r += ktg_P[0]*rcolor.r;
					rgb_color.g += ktg_P[0]*rcolor.g;
					rgb_color.b += ktg_P[0]*rcolor.b;
				}
			}
			else
			{
				int i = obj_num-NUM_SPHERES-1;
				thetha = asin((float)n1/(float)nt_P[i]);

				if(thetha < 41.2) /* Check if object is semi transparent */
				{
					points newFrom;
					newFrom.x = ipx;
					newFrom.y = ipy;
					newFrom.z = ipz;
			
					compute_refracted_ray(dx,dy,dz,nx,ny,nz,&rfx,&rfy,&rfz,nt_P[i]); /*Compute refracted ray*/
					rgb rcolor = Ray_Trace(obj_num,rfx,rfy,rfz,newFrom,depth+1);
					rgb_color.r += ktg_P[0]*rcolor.r;
					rgb_color.g += ktg_P[0]*rcolor.g;
					rgb_color.b += ktg_P[0]*rcolor.b;
				}
			}
		}
	}
	return rgb_color; 
}

/*****************************************************************************************

  Title:	RayTraceUtil

  Purpose:	Utility function to calculate reflected,refracted and shadow rays at each pixel

********************************************************************************************/

void RayTraceUtil()
{
	int	    xp, yp;
	float	xv, yv, dx, dy, dz;
	float   u, v;

	/*  Generate a ray for each pixel in the desired image.  */

	printf("Ray tracing...\n");
	for (xp=0; xp<xmax_pixel; xp++) 
	{
		u = (float)xp/xmax_pixel;

		for (yp=0; yp<ymax_pixel; yp++) 
		{
			v = (float)yp/ymax_pixel;

			/*  Compute the corresponding view port coordinates.  */

			xv = VXL + xp * xinterval;
			yv = VYB + yp * yinterval;

			/*  Compute the direction of the current ray from the "From" point to the 
			current position on the image.  */

			dx = ax*xv*tanv2 + bx*yv*tanv2 + cx;
			dy = ay*xv*tanv2 + by*yv*tanv2 + cy;
			dz = az*xv*tanv2 + bz*yv*tanv2 + cz;

			rgb rgb_color = Ray_Trace(0,dx,dy,dz,from,1);

			
			/* Save the computed color intensity to the image buffer. */

			texture_R[xp + xmax_pixel * yp] = rgb_color.r;
			texture_G[xp + xmax_pixel * yp] = rgb_color.g;
			texture_B[xp + xmax_pixel * yp] = rgb_color.b;

		}
	}
	/*  Write the image to the output file.  */

	printf("Writing to image...\n");
	fwrite(&xmax_pixel, sizeof(int), 1, outpfile);
	fwrite(&ymax_pixel, sizeof(int), 1, outpfile);
	
	fwrite(texture_R, sizeof(float), xmax_pixel*ymax_pixel, outpfile);
	fwrite(texture_G, sizeof(float), xmax_pixel*ymax_pixel, outpfile);
	fwrite(texture_B, sizeof(float), xmax_pixel*ymax_pixel, outpfile);

	fclose(outpfile);

}

/* Initialize the projection matrix.  */

void myinit(void)
{
	/* attributes */

    glClearColor(1.0, 1.0, 1.0, 1.0); /* white background */

	/* set up viewing */
	/* 512 x 512 window with origin lower left */

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 512.0, 0.0, 512.0);
    glMatrixMode(GL_MODELVIEW);
}

/* Display the ray traced image.   A more efficient method is 
   to use glDrawPixels(). */

void display( void )
{
	int s, t;
	float  r, g, b;

	glClear(GL_COLOR_BUFFER_BIT);  /*clear the window */

	for(t = 0; t < ymax_pixel; t++) {
		for(s = 0; s < xmax_pixel; s++) {

			r = texture_R[s + xmax_pixel * t];
			g = texture_G[s + xmax_pixel * t];
			b = texture_B[s + xmax_pixel * t];

			glColor3f(r, g, b);
			glBegin(GL_POINTS);
               glVertex2f(s,t); 
			glEnd();
		}
	 }

     glFlush(); /* clear buffers */
}

/*  Main routine.  */

int main(int argc, char**argv)
{

	Read_Information();
	Setup_Parameters();
	RayTraceUtil();

/* Standard GLUT initialization */

    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB); /* default, not needed */
    glutInitWindowSize(500,500); /* 500 x 500 pixel window */
    glutInitWindowPosition(0,0); /* place window top left on display */
    glutCreateWindow("Ray Trace"); /* window title */
    glutDisplayFunc(display); /* display callback invoked when window opened */

    myinit(); /* set attributes */
    glutMainLoop(); /* enter event loop */

	return(0);
}