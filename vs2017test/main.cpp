#define _CRT_SECURE_NO_WARNINGS
#include "glut.h"
#include <math.h>
#include <time.h>
#include <vector>
#include <windows.h>
#include <list>

using namespace std;

const int W = 600;
const int H = 600;
const double PI = 3.14;
const int GSZ = 150;

// texture matrices
const int  TSZ = 512;
unsigned char tx0[TSZ][TSZ][3]; // 3 stands for RGB
unsigned char tx1[1024][1024][3]; // 3 stands for RGB
unsigned char tx2[600][600][3]; // 3 stands for RGB
unsigned char tx3[1200][797][3];
unsigned char tx4[512][1024][3];
unsigned char tx5[1024][1024][3];


double ground[GSZ][GSZ] = { 0 };
double offset = 0;

unsigned char* bmp; // array of pixels

double eyex=-1, eyey=19, eyez=15;

double dx=0, dy=0, dz=0;
// ego-motion
double speed = 0, angular_speed = 0;
double sight_angle = PI; // initial sight
double pitch=-0.2;
double dir[3]; // x,y,z
// airplane
double tspeed = 0, turn_speed = 0;
double yaw = PI;
double apitch = 0;
double adir[3] = { sin(yaw),sin(apitch),cos(yaw) };
double ax = 0, ay = 12, az = 0;

bool is_building = true;
double hpitch = 0;
double vpitch = 0;

bool wasShot = false;
bool follow = true;
bool explosion = false;
double xi, yi, zi;
double vx = 50, vy=15, vz=50;
double xf =0 , yf = 0, zf = 0;
double yawi;
double hpitchi;
double vpitchi;
double t = 0, tExplostion = 0;
double a = -20.8;
int numOfTrees;


typedef struct {
	double x;
	double y;
	double z;
} TreeCoor ;

std::list<TreeCoor> treeList;

void DrawTree(double x, double y, double z);
void DrawCylinder(int n, double topr, double bottomr);
void UpdateTerrain1();
void UpdateTerrain2();
void Smooth();
void SetTexture(int);
void ReadPicture(char*);
void SetUpTrees(int);

void PrepareTextures()
{
	// texture definitions
	SetTexture(0); // brick wall
	glBindTexture(GL_TEXTURE_2D, 0); // setting texture #0
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, TSZ, TSZ, 0, GL_RGB, GL_UNSIGNED_BYTE, tx0);

	char name[30] = "camo_1.bmp";
	ReadPicture(name);
	SetTexture(1); // camo_2
	glBindTexture(GL_TEXTURE_2D, 1); // setting texture #1
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_UNSIGNED_BYTE, tx1);


	strcpy(name, "camo_2.bmp");
	ReadPicture(name);
	SetTexture(2); // camo_2
	glBindTexture(GL_TEXTURE_2D, 2); // setting texture #1
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 600, 600, 0, GL_RGB, GL_UNSIGNED_BYTE, tx2);


	strcpy(name, "camo_3.bmp");
	ReadPicture(name);
	SetTexture(3); // camo_3
	glBindTexture(GL_TEXTURE_2D, 3); // setting texture #3
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1200, 797, 0, GL_RGB, GL_UNSIGNED_BYTE, tx3);

	strcpy(name, "sky.bmp");
	ReadPicture(name);
	SetTexture(4); // sky
	glBindTexture(GL_TEXTURE_2D, 4); // setting texture #4
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, tx4);

	strcpy(name, "tree-1.bmp");
	ReadPicture(name);
	SetTexture(5); // sky
	glBindTexture(GL_TEXTURE_2D, 5); // setting texture #4
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_UNSIGNED_BYTE, tx5);


}
void ReadPicture(char* fname)
{
	FILE* pf;
	BITMAPFILEHEADER bf;
	BITMAPINFOHEADER bi;
	pf = fopen(fname, "rb");

	fread(&bf, sizeof(bf), 1, pf);
	fread(&bi, sizeof(bi), 1, pf);
	bmp = (unsigned char*)malloc(3 * bi.biWidth * bi.biHeight); // 3 stands for BGR
	fread(bmp, 1, 3 * bi.biWidth * bi.biHeight, pf);
	fclose(pf);
}

void SetTexture(int txnum)
{
	int i, j, k;
	switch (txnum)
	{
	case 0: // brick wall
		for (i = 0; i < TSZ; i++)
			for (j = 0; j < TSZ; j++)
			{
				k = rand() % 20;
				if (i % (TSZ / 2) < 20 ||  // horizontal separation
					i > TSZ / 2 && j % (TSZ / 2) < 20 || // upper part of bricks
					i <= TSZ / 2 && (j > TSZ / 4 && j < TSZ / 4 + 20 || j>3 * TSZ / 4 && j < 3 * TSZ / 4 + 20)  // lower paert of bricks
					)
				{    // separating layer
					tx0[i][j][0] = 150 + k; //red
					tx0[i][j][1] = 150 + k; //green 
					tx0[i][j][2] = 150 + k; //blue
				}
				else // brick
				{
					tx0[i][j][0] = 150 + k; //red
					tx0[i][j][1] = 50 + k; //green 
					tx0[i][j][2] = k; //blue
				}
			}
		break;
	case 1:
		for (i = 0; i < 1024; i++)
			for (j = 0; j < 1024; j++)
			{
				tx1[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx1[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx1[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;
	case 2:
		for (i = 0; i < 600; i++)
			for (j = 0; j < 600; j++)
			{
				tx2[i][j][0] = bmp[(i * 600 + j) * 3 + 2];  //red
				tx2[i][j][1] = bmp[(i * 600 + j) * 3 + 1];// green
				tx2[i][j][2] = bmp[(i * 600 + j) * 3];//blue
			}
		break;
	case 3:
		for (i = 0; i < 797; i++)
			for (j = 0; j < 1200; j++)
			{
				tx3[i][j][0] = bmp[(i * 1200 + j) * 3 + 2];  //red
				tx3[i][j][1] = bmp[(i * 1200 + j) * 3 + 1];// green
				tx3[i][j][2] = bmp[(i * 1200 + j) * 3];//blue
			}
		break;
	case 4:
		for (i = 0; i < 512; i++)
			for (j = 0; j < 1024; j++)
			{
				tx4[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx4[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx4[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;
	case 5:
		for (i = 0; i < 1024; i++)
			for (j = 0; j < 1024; j++)
			{	
				tx5[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx5[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx5[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;

	}
}

void DrawTextureCylinder(int n, double topr, double bottomr, int tnum, double hrepeats,
	double txbottom, double txtop);

void DrawPitchControlHorizontal() // 2D
{
	// background

	glBegin(GL_POLYGON);
	glColor3d(0.1, 1, 0.3);
	glVertex2d(-1, -1);
	glVertex2d(1, -1);
	glColor3d(0.3, 1, 0.8);
	glVertex2d(1, 1);
	glVertex2d(-1, 1);
	glEnd();

	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2d(-1, 0);
	glVertex2d(1, 0);
	glEnd();

	// slider
	glColor3d(0.5, 0.5, 0.5);
	glBegin(GL_POLYGON);
	glVertex2d(hpitch - 0.2, -0.2);
	glVertex2d(hpitch + 0.2, -0.2);
	glVertex2d(hpitch + 0.2, 0.2);
	glVertex2d(hpitch - 0.2, 0.2);
	glEnd();
	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2d(hpitch, -0.2);
	glVertex2d(hpitch, 0.2);
	glEnd();
}


void DrawPitchControlVertical() // 2D
{
	// background

	glBegin(GL_POLYGON);
	glColor3d(0.9, 0, 0);
	glVertex2d(-1, 1);
	glVertex2d(1, 1);
	glColor3d(0.9, 0.4, 0.6);
	glVertex2d(1, -1);
	glVertex2d(-1, -1);
	glEnd();

	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2d(0, -1);
	glVertex2d(0, 1);
	glEnd();

	// slider
	glColor3d(0.5, 0.5, 0.5);
	glBegin(GL_POLYGON);
	glVertex2d(-0.2, vpitch - 0.2);
	glVertex2d(-0.2, vpitch + 0.2);
	glVertex2d(0.2, vpitch + 0.2);
	glVertex2d(0.2, vpitch - 0.2);
	glEnd();
	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2d(-0.2, vpitch);
	glVertex2d(0.2, vpitch);
	glEnd();



}

// draw trapezoid
void DrawTrapezoid(int textureId) {
	int counter = 0;
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, textureId); // texture #0
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);// GL_REPLACE);// GL_MODULATE); 
	glPushMatrix();

	// front blate
	glBegin(GL_POLYGON);
	glTexCoord2d(0, 1);  glVertex3d(0, 0, -2);
	glTexCoord2d(4, 2);	 glVertex3d(5, 0, -2);
	glTexCoord2d(4, 4);  glVertex3d(5, 5, 0);
	glTexCoord2d(2, 4);  glVertex3d(0, 5, 0);
	glEnd();

	// back plate
	glBegin(GL_POLYGON);
	glTexCoord2d(1, 2);  glVertex3d(0, 0, 10);
	glTexCoord2d(4, 2);  glVertex3d(5, 0, 10);
	glTexCoord2d(4, 4);  glVertex3d(5, 5, 8);
	glTexCoord2d(2, 1);  glVertex3d(0, 5, 8);
	glEnd();


	// bottom plate
	glBegin(GL_POLYGON);
	glTexCoord2d(3, 6); glVertex3d(0, 0, -2);
	glTexCoord2d(6, 3); glVertex3d(5, 0, -2);
	glTexCoord2d(6, 6); glVertex3d(5, 0, 10);
	glTexCoord2d(2, 4); glVertex3d(0, 0, 10);
	glEnd();

	//upper plate
	glBegin(GL_POLYGON);
	glTexCoord2d(6, 3);  glVertex3d(0, 5, 0);
	glTexCoord2d(3, 6);  glVertex3d(5, 5, 0);
	glTexCoord2d(2, 4);  glVertex3d(5, 5, 8);
	glTexCoord2d(6, 6);  glVertex3d(0, 5, 8);
	glEnd();

	//fill
	//glColor3d(0, 0, 0);
	glBegin(GL_POLYGON);
	glTexCoord2d(4, 4);  glVertex3d(0.75, 0, -2); // x = 0 ?
	glTexCoord2d(2, 4);  glVertex3d(0.75, 0, 10);
	glTexCoord2d(2, 2);  glVertex3d(0.75, 5, 8);
	glTexCoord2d(4, 2);  glVertex3d(0.75, 5, 0);
	glEnd();


	glBegin(GL_POLYGON);
	glTexCoord2d(4, 4);  glVertex3d(4.75, 0, -2);  // x = 5 ?
	glTexCoord2d(2, 4);  glVertex3d(4.75, 0, 10);
	glTexCoord2d(2, 2);  glVertex3d(4.75, 5, 8);
	glTexCoord2d(4, 2);  glVertex3d(4.75, 5, 0);
	glEnd();

	glPopMatrix();
	glDisable(GL_TEXTURE_2D);
}


// draw 3d tank model
void DrawTank() {
	glColor3d(1, 1, 1);

	glPushMatrix();
		glPushMatrix();
			glRotated(-hpitch * 180 * 0.65 / PI, 0, 1, 0);


			// upper section of the tank
			glPushMatrix();
			glTranslated(0, 0, -5);
			
			DrawTrapezoid(2);
			glPopMatrix();

			glPushMatrix();
				glRotated(vpitch * 180 * 0.6 / PI, 1, 0, 0);

				// cannon tube
				glPushMatrix();
				glTranslated(2.5, 2, -14);
				glRotated(90, 1, 0, 0);
				glScaled(1, 14, 1);
				
				DrawTextureCylinder(30, 0.5, 0.8, 3, 1, 1, 1);

				glPopMatrix();


				// cannon end 
				glPushMatrix();
				glTranslated(2.5, 2, -15);
				glRotated(90, 1, 0, 0);
				DrawTextureCylinder(30, 0.8, 1.2, 3, 10, 0, 10);
				glPopMatrix();


				glPushMatrix();
				glTranslated(2.5, 2, -16);
				glRotated(90, 1, 0, 0);
				DrawTextureCylinder(30, 1.2, 1.2, 3, 10, 0, 10);
				glPopMatrix();
			glPopMatrix();
		glPopMatrix();
		// lower section of the tank
		glPushMatrix();
		glTranslated(-2, -0.1, 8.8);
		glScaled(2, 1, 1.8);
		glRotated(180, 1, 0, 0);
		glColor3d(1, 1, 1);
		DrawTrapezoid(1);
		glPopMatrix();

		// left side wheels
		for (int z = -1; z < 3; z++) {
			glPushMatrix();
			glTranslated(0, -3, z * 3.75);
			glRotated(90, 0, 0, 1);
			glScaled(2, 2, 2);
			DrawCylinder(30, 0.5, 1);
			glPopMatrix();
		}

		// right side wheels
		for (int z = -1; z < 3; z++) {
			glPushMatrix();
			glTranslated(6, -3, z * 3.75);
			glRotated(90, 0, 0, -1);
			glScaled(2, 2, 2);
			DrawCylinder(30, 0.5, 1);
			glPopMatrix();
		}
	glPopMatrix();

}


void init()
{
	int i,j;
	srand(time(0));
	numOfTrees = 2 + rand() % 3;
	//             red  green  blue
	glClearColor(0.0, 0.3, 0.6, 0);// color of window background
	glEnable(GL_DEPTH_TEST); // show the nearest object

	for (i = 0; i < GSZ; i++)
		for (j = 0; j < GSZ; j++)
			ground[i][j] = 0.1;//3*sin(j/3.0)+2*sin(i/2.0);

	for (i = 1; i <= 1000; i++)
		UpdateTerrain1();
	for (i = 1; i <= 1200;i++)
		UpdateTerrain2();
	Smooth();
	Smooth();
	
	
	glEnable(GL_NORMALIZE); // used for lighting
	PrepareTextures();

	SetUpTrees(numOfTrees);
}

// randomize trees coordiantes, num of trees determined in init
void SetUpTrees(int numOfTrees)
{
	for (int i = 0; i < numOfTrees; i++)
	{
		double x = -30 + rand() % 30;
		double y = 2 + rand() % 2;
		double z = -30 + rand() % 30;
		TreeCoor tc;
		tc.x = x; tc.y = y; tc.z = z;
		treeList.push_back(tc);
	}

	for (TreeCoor tc : treeList)
	{
		printf("%lf %lf %lf\n", tc.x, tc.y, tc.z);
	}
}


// draw all trees
void DrawTrees()
{
	for (TreeCoor tc : treeList)
	{
		DrawTree(tc.x, tc.y, tc.z);
	}
}

void UpdateTerrain()
{
	int i, j;
	double delta = 0.05;

	for (i = 0; i < GSZ; i++)
		for (j = 0; j < GSZ; j++)
		{
			if (rand() % 100 >= 50)
				delta = -delta;
			ground[i][j] += delta;
		}
}

void UpdateTerrain1()
{
	int i, j;
	double delta = 0.08,a,b;
	int x1, y1, x2, y2;

	if (rand() % 2 == 0)
		delta = -delta;

	x1 = rand() % GSZ;
	y1 = rand() % GSZ;
	x2 = rand() % GSZ;
	y2 = rand() % GSZ;

	if (x1 == x2) // exception
		return;

	a = (y2 - y1) / (double)(x2 - x1);
	b = y1 - a * x1;

	for (i = 0; i < GSZ; i++)
		for (j = 0; j < GSZ; j++)
		{
			if (i < j * a + b)
				ground[i][j] += delta;
			else
				ground[i][j] -= delta;
		}
}

// random walk
void UpdateTerrain2()
{
	int steps = 5000;
	double delta = 0.01;
	int x, y;
	if (rand() % 2 == 0)
		delta = -delta;

	x = rand() % GSZ;
	y = rand() % GSZ;

	for (int counter = 1; counter <= steps; counter++)
	{
		ground[y][x] += delta;
		// choose new direction
		switch (rand() % 4)
		{
		case 0: // up
			y++;
			break;
		case 1: // down
			y--;
			break;
		case 2: // right
			x++;
			break;
		case 3: // left
			x--;
			break;
		}
		// in case that we go out of matrix we use the formula:
		x = (x+GSZ)%GSZ; 
		y = (y + GSZ) % GSZ;
	}
}

void Smooth()
{
	int i, j;
	double tmp[GSZ][GSZ];

	for (i = 1; i < GSZ - 1; i++)
		for (j = 1; j < GSZ - 1; j++)
			tmp[i][j] = (	ground[i-1][j-1]+ 2*ground[i-1][j]+ground[i-1][j+1]+
								2*ground[i][j - 1] + 4*ground[i][j] + 2*ground[i][j + 1]+
								ground[i +1][j - 1] + 2*ground[i + 1][j] + ground[i + 1][j + 1]	) / 16.0;
// copy tmp to ground
	for (i = 1; i < GSZ - 1; i++)
		for (j = 1; j < GSZ - 1; j++)
			ground[i][j] = tmp[i][j];
}

void SetColor(double h)
{
	
//	glColor3d(fabs(h / 3), (h+3)/6, fabs(sin(h)));
	h = fabs(h);
	if(h<0.15) // sand
		glColor3d(0.8, 0.7, 0.5);
	else
	if(h<3) // grass
		glColor3d(0.1+h/11, 0.5-h/10, 0);
	else // snow
		glColor3d(h/5,h/5,h/4);
		
}

void SetNormal(int row, int col)
{
	// if  it was defined glEnable(GL_NORMALIZE) then we don't need
	// to normalize it manually
	glNormal3d(ground[row][col + 1] - ground[row][col], 1, 
												ground[row + 1][col] - ground[row][col]);
}


void DrawFloor()
{
	int i, j;

	glColor3d(0.8, 0.8, 1); // blue

	for (i = 1; i < GSZ - 2; i++)
		for (j = 1; j < GSZ - 2; j++)
		{
			glBegin(GL_POLYGON); // lines parallel to X
			SetColor(ground[i][j]);
			//			SetNormal(i, j);
			glVertex3d(j - GSZ / 2, ground[i][j], i - GSZ / 2);
			SetColor(ground[i - 1][j]);
			//			SetNormal(i - 1, j);
			glVertex3d(j - GSZ / 2, ground[i - 1][j], i - 1 - GSZ / 2);
			SetColor(ground[i - 1][j - 1]);
			//			SetNormal(i - 1, j - 1);
			glVertex3d(j - 1 - GSZ / 2, ground[i - 1][j - 1], i - 1 - GSZ / 2);
			SetColor(ground[i][j - 1]);
			//			SetNormal(i, j - 1);
			glVertex3d(j - 1 - GSZ / 2, ground[i][j - 1], i - GSZ / 2);
			glEnd();
		}


	// water
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4d(0.2, 0.4, 0.6, 0.7);
	glBegin(GL_POLYGON);
	glVertex3d(-GSZ / 2, 0, -GSZ / 2);
	glVertex3d(-GSZ / 2, 0, GSZ / 2);
	glVertex3d(GSZ / 2, 0, GSZ / 2);
	glVertex3d(GSZ / 2, 0, -GSZ / 2);
	glEnd();
	glDisable(GL_BLEND);
}

// topr and bottom r are the top and bottom radiuses
void DrawCylinder(int n, double topr, double bottomr)
{
	double alpha, teta = 2 * PI / n;

	for (alpha = 0; alpha <= 2 * PI; alpha += teta)
	{
		glBegin(GL_POLYGON);
		glColor3d(fabs(sin(alpha)/2) , fabs(sin(alpha)), fabs(cos(alpha))/2);
		glVertex3d(topr*sin(alpha), 1,topr* cos(alpha)); //1
	
		glVertex3d(topr*sin(alpha + teta), 1, topr*cos(alpha + teta));//2
		glColor3d(fabs(sin(alpha)/2), fabs(sin(alpha)), fabs(cos(alpha))/2);
		glVertex3d(bottomr*sin(alpha + teta), 0, bottomr * cos(alpha + teta));//3
		glVertex3d(bottomr * sin(alpha), 0, bottomr * cos(alpha));// 4
		glEnd();
	}
}
// topr and bottom r are the top and bottom radiuses
void DrawTexCylinder1(int n, double topr, double bottomr,int tnum,double hrepeats)
{
	double alpha, teta = 2 * PI / n;
	int side_counter;
	double tex_part = hrepeats / n; // part of a texture that covers one side
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1); // texture #0
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);// GL_REPLACE);// GL_MODULATE); 

	for (alpha = 0,side_counter=0; alpha <= 2 * PI; alpha += teta,side_counter++)
	{
		glBegin(GL_POLYGON);
		glColor3d(0.5+fabs(sin(alpha))/2 , 0.5 + fabs(sin(alpha))/2 , 0.5 +fabs( sin(alpha)/2 ));
		glTexCoord2d(side_counter* tex_part, 1);	glVertex3d(topr * sin(alpha), 1, topr * cos(alpha)); //1
//		glColor3d(0.2 + fabs(sin(alpha + teta) / 2), 0.2 + fabs(sin(alpha + teta) / 2), 0.2 +fabs( sin(alpha + teta) / 2));
		glTexCoord2d((side_counter+1) * tex_part, 1);		glVertex3d(topr * sin(alpha + teta), 1, topr * cos(alpha + teta));//2
		glTexCoord2d((side_counter + 1) * tex_part, 0);		glVertex3d(bottomr * sin(alpha + teta), 0, bottomr * cos(alpha + teta));//3
//		glColor3d(0.2 + fabs(sin(alpha) / 2), 0.2 +fabs( sin(alpha) / 2), 0.2 +fabs( sin(alpha) / 2));
		glTexCoord2d(side_counter * tex_part, 0);		glVertex3d(bottomr * sin(alpha), 0, bottomr * cos(alpha));// 4
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);

}

// topr and bottom r are the top and bottom radiuses
void DrawTextureCylinder(int n, double topr, double bottomr, int tnum, double hrepeats,
	double txbottom, double txtop)
{
	double alpha, teta = 2 * PI / n;
	int side_counter;
	double tex_part = hrepeats / n; // part of a texture that covers one side
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, tnum); // texture #0
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE,  GL_REPLACE);// GL_MODULATE); 

	for (alpha = 0, side_counter = 0; alpha <= 2 * PI; alpha += teta, side_counter++)
	{
		glBegin(GL_POLYGON);
		glColor3d(0.5 + fabs(sin(alpha)) / 2, 0.5 + fabs(sin(alpha)) / 2, 0.5 + fabs(sin(alpha) / 2));
		glTexCoord2d(side_counter * tex_part, txtop);	glVertex3d(topr * sin(alpha), 1, topr * cos(alpha)); //1
//		glColor3d(0.2 + fabs(sin(alpha + teta) / 2), 0.2 + fabs(sin(alpha + teta) / 2), 0.2 +fabs( sin(alpha + teta) / 2));
		glTexCoord2d((side_counter + 1) * tex_part, txtop);		glVertex3d(topr * sin(alpha + teta), 1, topr * cos(alpha + teta));//2
		glTexCoord2d((side_counter + 1) * tex_part, txbottom);		glVertex3d(bottomr * sin(alpha + teta), 0, bottomr * cos(alpha + teta));//3
//		glColor3d(0.2 + fabs(sin(alpha) / 2), 0.2 +fabs( sin(alpha) / 2), 0.2 +fabs( sin(alpha) / 2));
		glTexCoord2d(side_counter * tex_part, txbottom);		glVertex3d(bottomr * sin(alpha), 0, bottomr * cos(alpha));// 4
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);

}

// topr and bottom r are the top and bottom radiuses
void DrawLtCylinder1(int n, double topr, double bottomr)
{
	double alpha, teta = 2 * PI / n;

	for (alpha = 0; alpha <= 2 * PI; alpha += teta)
	{
		glBegin(GL_POLYGON);
		glNormal3d(sin(alpha),bottomr*(bottomr-topr) , cos(alpha));
		glVertex3d(topr * sin(alpha), 1, topr * cos(alpha)); //1
		glNormal3d(sin(alpha+teta), bottomr * (bottomr - topr), cos(alpha+teta));

		glVertex3d(topr * sin(alpha + teta), 1, topr * cos(alpha + teta));//2
		glVertex3d(bottomr * sin(alpha + teta), 0, bottomr * cos(alpha + teta));//3
		glNormal3d(sin(alpha), bottomr * (bottomr - topr), cos(alpha));
		glVertex3d(bottomr * sin(alpha), 0, bottomr * cos(alpha));// 4
		glEnd();
	}
}

void DrawCylinder(int n)
{
	double alpha, teta = 2 * PI / n;
	
	for (alpha = 0; alpha <= 2 * PI; alpha += teta)
	{
		glBegin(GL_POLYGON);
		glColor3d(fabs(sin(alpha))/2, (1+cos(alpha))/3, fabs(cos(alpha+PI/3))/2);
		glVertex3d(sin(alpha), 1, cos(alpha)); //1
		glColor3d(fabs(sin(alpha+teta)) / 2, (1 + cos(alpha+teta)) / 3, fabs(cos(alpha+teta + PI / 3)) / 2);
		glVertex3d(sin(alpha + teta), 1, cos(alpha + teta));//2
		glColor3d(fabs(sin(alpha)), (1 + cos(alpha)), fabs(cos(alpha)) );
		glVertex3d(sin(alpha +teta), 0, cos(alpha + teta));//3
		glVertex3d(sin(alpha  ), 0, cos(alpha ));// 4
		glEnd();
	}
}
// hrep are horizontal repeats of texture
// vrep are vertical repeats of texture

void DrawTexSphere(int slices, int stacks,int tnum,double hrep,double vrep)
{
	double beta = PI / slices, gamma; // gamma is running angle from -90 to 90
	int i;
	double vpart = vrep / stacks;

	for (gamma = -PI / 2, i = 0; gamma <= PI / 2; gamma += beta, i++)
	{
		glPushMatrix();
		glTranslated(0, sin(gamma), 0);
		glScaled(1, sin(gamma + beta) - sin(gamma), 1);
		DrawTextureCylinder(stacks, cos(gamma + beta), cos(gamma),tnum,hrep,vpart*i,vpart*(i+1));
		glPopMatrix();
	}
}

// draws skybox with predetermined texture
void DrawSkyBox() {
	glPushMatrix();
	glRotated(offset / 20, 0, 1, 0);
	glTranslated(0, 30, 0);
	glScaled(100, 100, 100);
	DrawTexSphere(200, 200, 4, 1, 1);
	glPopMatrix();
}

// handle shooting and explosions visuals
void HandleShooting() {
	if (wasShot && !explosion) {
		glPushMatrix();
		glColor3d(0, 0, 0);
		glTranslated(xf, yf, zf);

		glRotated(180, 0, 1, 0);
		
		glRotated(yawi * 180 / PI, 0, 1, 0);

		glTranslated(0, 0, -5);
		glTranslated(2.5, 2, -12);
		glutSolidSphere(1, 20, 20);
		glPopMatrix();
	}
	

	if (explosion) {
		glPushMatrix();
		glColor3d(0.9, cos(tExplostion/8), sin(tExplostion/8));
		glTranslated(xf, yf, zf);
		glRotated(180, 0, 1, 0);
		glRotated(yawi * 180 / PI, 0, 1, 0);
		glTranslated(0, 0, -5);
		glTranslated(2.5, 2, -12);
		glutSolidSphere(2 + tExplostion, 20, 20);
		glPopMatrix();
		

	}
}

// draws a single tree with textures
void DrawTree(double x, double y, double z)
{
	glPushMatrix();
	glTranslated(x, y, z);
	glScaled(1, 16, 1);
	DrawTextureCylinder(10, 2, 3.2, 5, 2, 0, 10);
	glScaled(1, 0.8, 1);
	glTranslated(0, 1, 0);
	DrawTextureCylinder(10, 2, 2, 5, 2, 0, 10);

	glTranslated(0, 1.1, 0);
	glScaled(1, 0.1, 1);
	glColor3d(0, 0.5, 0.3);
	glutSolidSphere(6, 30, 30);
	glPopMatrix();
}

// default display
void regularDisplay()
{
	// clean frame buffer and Z-buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, W, H );
	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity(); 
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(eyex,eyey,eyez,  // eye
		eyex+dir[0], eyey+sin(pitch), eyez+dir[2],// center or PointOfInterest
		 0, 1, 0); // up
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix

	DrawSkyBox();
	DrawFloor();
	glEnable(GL_DEPTH_TEST);
	glPushMatrix();
	glTranslated(ax, ay, az);
	glRotated(180, 0, 1, 0);
	glRotated(yaw * 180 / PI, 0, 1, 0);
	DrawTank();
	glPopMatrix();

	HandleShooting();

	DrawTrees();

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST); // last painted occludes preveous painted (2D)

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity(); // start transformations from the start
	glOrtho(-1, 1, -1, 1, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glViewport(0, 0, 100, 100);
	DrawPitchControlHorizontal();

	glViewport(W - 100, 0, 100, 100);
	DrawPitchControlVertical();

	
	glutSwapBuffers(); // show all
}

// top down view
void topDisplay()
{
	// clean frame buffer and Z-buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, W, H);

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(eyex, eyey+40, eyez,  // eye
		eyex , eyey  , eyez - 0.1,// center or PointOfInterest
		0, 1, 0); // up

	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix

	glEnable(GL_TEXTURE_2D);
	glEnable(GL_DEPTH_TEST);


	DrawSkyBox();
	DrawFloor();
	glEnable(GL_DEPTH_TEST);
	glPushMatrix();
	glTranslated(ax, ay, az);
	glRotated(180, 0, 1, 0);
	glRotated(yaw * 180 / PI, 0, 1, 0);
	glRotated(turn_speed, 0, 1, 0);
	DrawTank();
	glPopMatrix();

	HandleShooting();

	DrawTrees();

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST); // last painted occludes preveous painted (2D)

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity(); // start transformations from the start
	glOrtho(-1, 1, -1, 1, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glViewport(0, 0, 100, 100);
	DrawPitchControlHorizontal();

	glViewport(W - 100, 0, 100, 100);
	DrawPitchControlVertical();
	glutSwapBuffers(); // show all
}

// follows tank cannon view
void gunnerDisplay()
{
	
	double h_angle = -hpitch * 180 * 0.2 / PI;
	double tdir[3] = { sin(sight_angle), sin(vpitch), cos(sight_angle) }; // x
	// clean frame buffer and Z-buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, W, H);

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(ax + 3 * adir[0], ay + 6 + 3 * adir[1], az + 2 * adir[2],  // eye is now placed into airplane
		ax + 4 * adir[0] + hpitch, ay + 6 + 4 * adir[1] + 0.5 * vpitch , az + 3 * adir[2] + hpitch,// center or PointOfInterest
		0, 1, 0); // up

	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix

	glEnable(GL_TEXTURE_2D);
	glEnable(GL_DEPTH_TEST);

	DrawSkyBox();
	DrawFloor();
	glEnable(GL_DEPTH_TEST);
	glPushMatrix();

	glTranslated(ax, ay, az);
	glRotated(180, 0, 1, 0);
	glRotated(yaw * 180 / PI, 0, 1, 0);
	DrawTank();
	glPopMatrix();
	
	HandleShooting();

	DrawTrees();

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST); // last painted occludes preveous painted (2D)

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity(); // start transformations from the start
	glOrtho(-1, 1, -1, 1, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glViewport(0, 0, 100, 100);
	DrawPitchControlHorizontal();

	glViewport(W - 100, 0, 100, 100);
	DrawPitchControlVertical();
	glutSwapBuffers(); // show all
}


// all the views combined - each view type has its own viewport
void combinedDisplay()
{
	// clean frame buffer and Z-buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// regular view
	glViewport(0, H / 2, W / 2, H / 2);

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(eyex, eyey, eyez,  // eye
		eyex + dir[0], eyey + sin(pitch), eyez + dir[2],// center or PointOfInterest
		0, 1, 0); // up

	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix


	DrawSkyBox();
	DrawFloor();
	glEnable(GL_DEPTH_TEST);
	glPushMatrix();
	glTranslated(ax, ay, az);
	glRotated(180, 0, 1, 0);
	glRotated(yaw * 180 / PI, 0, 1, 0);
	DrawTank();
	glPopMatrix();
	HandleShooting();
	DrawTrees();

	// top view
	glViewport(W / 2, H / 2, W / 2, H / 2);
	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(eyex, eyey + 40, eyez,  // eye
		eyex, eyey, eyez - 0.1,// center or PointOfInterest
		0, 1, 0); // up

	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix

	DrawSkyBox();
	DrawFloor();
	glPushMatrix();
	glTranslated(ax, ay, az);
	glRotated(180, 0, 1, 0);
	glRotated(yaw * 180 / PI, 0, 1, 0);//yaw is in radians so we transform it to degrees
	glRotated(turn_speed, 0, 0, -1); // roll around main axis of airplane
	DrawTank();
	glPopMatrix();

	HandleShooting();
	DrawTrees();


	// cannon view
	glViewport(0, 0, W, H / 2);

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(ax + 3 * adir[0], ay + 6 + 3 * adir[1], az + 2 * adir[2],  // eye is now placed into airplane
		ax + 4 * adir[0] + hpitch, ay + 6 + 4 * adir[1] + 0.5 * vpitch, az + 3 * adir[2] + hpitch,// center or PointOfInterest
		0, 1, 0); // up

	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix

	DrawSkyBox();
	DrawFloor();

	glPushMatrix();
	glTranslated(ax, ay, az);
	glRotated(yaw * 180 / PI, 0, 1, 0);//yaw is in radians so we transform it to degrees
	glRotated(turn_speed, 0, 0, -1); // roll around main axis of airplane
	glRotated(180, 0, 1, 0);
	DrawTank();
	glPopMatrix();

	HandleShooting();
	DrawTrees();

	// cannon control horizontal 
	glViewport(0, 0, 100, 100);
	glDisable(GL_DEPTH_TEST); // last painted occludes preveous painted (2D)

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity(); // start transformations from the start
	glOrtho(-1, 1, -1, 1, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	DrawPitchControlHorizontal();

	// cannon control vertical 
	glViewport(W - 100, 0, 100, 100);
	DrawPitchControlVertical();


	glEnable(GL_DEPTH_TEST); // near occludes far (3D)
	glutSwapBuffers(); // show all
}


// logic process
void idle()
{
	t += 0.01;
	if (wasShot) {
		// calculates distance for x,y,z positions after t seconds of shot
		xf = xi + sin(yawi) * t * vx - hpitchi * a * t * t;
		yf = yi + (vy * t * sin(vpitchi)) + (a * t * t);
		zf = zi + cos(yawi) * t * vz;
		
		int i, j;

		// iterate over height terrain pixels and check if shot hit terrain height
		for (i = 0; i < GSZ;  i++) {
			bool collider = false;
			for (j = 0; j < GSZ; j++) {
				if (ground[i][j] >= yf) {
					wasShot = false;
					collider = true;
					break;
				}
			}

			// collision occured, trigger explosion
			if (collider) {
				explosion = true;
				tExplostion = 0;
			}		
		}		
	}

	// explosion timer
	if (explosion) {
		tExplostion += 0.2;
		if (tExplostion >= 10) {
			explosion = false;
		}
	}

	offset += 0.5;


	// ego motion
	sight_angle += angular_speed;

	dir[0] = sin(sight_angle); // x
	dir[1] = sin(pitch);   //y
	dir[2] = cos(sight_angle); // z

	eyex += speed*dir[0];
	eyey += speed*dir[1];
	eyez += speed*dir[2];

	yaw += turn_speed;
	adir[0] = sin(yaw);
	adir[1] = sin(apitch);
	adir[2] = cos(yaw);

	ax += tspeed * adir[0];
	ay += tspeed * adir[1];
	az += tspeed * adir[2];
	glutPostRedisplay(); // indirect call to display
}

void SpecialKey(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_UP:
		speed+=0.001;
		break;
	case GLUT_KEY_DOWN:
		speed-= 0.001;
		break;
	case GLUT_KEY_LEFT:
		angular_speed+=0.001;
		break;
	case GLUT_KEY_RIGHT:
		angular_speed -= 0.001;
		break;
	case GLUT_KEY_PAGE_UP:
		pitch += 0.1;
		break;
	case GLUT_KEY_PAGE_DOWN:
		pitch -= 0.1;
		break;
	}
}


void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'a':
		turn_speed += 0.001;
		break;
	case 'd':
		turn_speed -= 0.001;
		break;
	case 'w':
		tspeed += 0.006;
		break;
	case 's':
		tspeed -= 0.006;
		break;
	case ' ':
		if (!explosion && !wasShot) {
			t = 0;
			wasShot = true;
			yawi = yaw;
			hpitchi = hpitch;
			vpitchi = vpitch;
			xi = ax + 10 * hpitchi;
			yi = ay + 10 * vpitchi;
			zi = az;
		}
		break;
	}
}

void menu(int choice)
{
	switch (choice)
	{
	case 1:// regular view
		glutDisplayFunc(regularDisplay);
		break;
	case 2:// top view
		glutDisplayFunc(topDisplay);
		break;
	case 3:// cockpit view
		glutDisplayFunc(gunnerDisplay);
		break;
	case 4:// combined view
		glutDisplayFunc(combinedDisplay);
		break;

	}
}



void mouse_motion(int x, int y)
{
	double xh, yh;

	
	// for left pitch (horizontal)
	xh = (2.0 * x / 100) - 1; 
	yh = (2.0 * (H - y) / 100) - 1;
	if (-0.8 < xh && xh < 0.8 && hpitch - 0.8 < yh && yh < hpitch + 0.8)
		if (fabs(xh < 1))
			hpitch = xh;

	// for right pitch (vertical)
	double xv, yv;
	xv = (2.0 * x / W) - 1; // xxx is in range (-1,1)
	yv = (2.0 * (H - y) / 100) - 1;

	if (0.8 < xv && xv < 0.86 && 0 < yv && yv < 0.8)
		if (fabs(yv) < 1)
			vpitch = yv;
}

void main(int argc, char* argv[])
{
	glutInit(&argc, argv); // 
	// defines matrices: 1. Color matrix (frame buffer), 
	// 2. video buffer
	// 3. Z-buffer
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE|GLUT_DEPTH); 
	glutInitWindowSize(W, H); // physical size of window (in pixels)
	glutInitWindowPosition(200, 100);
	glutCreateWindow("3D Tank Model");

	glutDisplayFunc(regularDisplay); //  display is a refresh window function
	glutIdleFunc(idle); // kind of timer function
	//glutMouseFunc(mouse);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(SpecialKey);
	glutMotionFunc(mouse_motion);
	// menu
	glutCreateMenu(menu);
	glutAddMenuEntry("Regular View", 1);
	glutAddMenuEntry("Top View", 2);
	glutAddMenuEntry("Gunner View", 3);
	glutAddMenuEntry("Combined View", 4);
	glutAttachMenu(GLUT_RIGHT_BUTTON);


	init();

	glutMainLoop();
}