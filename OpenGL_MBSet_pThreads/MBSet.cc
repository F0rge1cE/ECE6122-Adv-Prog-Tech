// Calculate and display the Mandelbrot set
// ECE4893/8893 final project, Fall 2011

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#include <GL/glut.h>
#include <GL/glext.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "complex.h"

using namespace std;

#define WindowSize 512	// Use the macrodefinition in order to declare the Iteration_Count array.
#define Max_Depth 20 	// The maximum depth of expanding the MBSet figure.

// Min and max complex plane values
Complex  minC(-2.0, -1.2);
Complex  maxC( 1.0, 1.8);
int      maxIt = 2000;     // Max iterations for the set computations
// int WindowSize = 512;
int nThreads = 16;	// Use 16 threads to compute
int Real_per_thread =  WindowSize / nThreads;
int Imag_per_thread =  WindowSize;	
// Divide 512 real values into 16 parts. Each thread computes 32 reals * 512 imags
double Real_Step = (maxC.real - minC.real) / (WindowSize);
double Imag_Step = (maxC.imag - minC.imag) / (WindowSize);

int Iteration_Count[WindowSize * WindowSize];	// Store each pixel's iteration count in an array
int updateRate = 5;		// Update setting
int It_Threshold = 6;	
// Threshold for drawing the MBSet. Points iterated less than Threshold should be white. 

// Mutex variables 
pthread_mutex_t exitMutex; // Used for creating and waiting
pthread_cond_t  exitCond;
pthread_mutex_t ActiveCountMutex;
int      ActiveCount = 0;   // Counter of active threads



class Color
{
	public:
  		Color();
  		Color(double r0, double g0, double b0);
	public:
  		double r;
  		double g;
  		double b;
};
Color::Color()
	: r(0), g(0), b(0)
{
}
Color::Color(double r0, double g0, double b0)
    : r(r0), g(g0), b(b0)
{
}


Color* colors;	// Array for storing colors of caertain number of iterations
void Init_Colors()
{
	drand48();	
	// Initialize the psedo-random function to avoid the unexpected black area
	// Or the first randomly generated RGB will be too small (looks like black)

  	colors = new Color[maxIt + 1];
  	for (int i = 0; i < maxIt; ++i)
    {	
      colors[i]=Color(1, 1, 1);
      if (i < It_Threshold)
        { 
          	colors[i] = Color(1, 1, 1);
        }
      else
        {
          	colors[i] = Color(drand48(), drand48(), drand48());	// RGB values are between 0.0 to 1.0, non-negative
        }
    }
  	colors[maxIt] = Color(0, 0, 0);  // Black for Z in MBSet, maxIt = 2000
  	colors[It_Threshold].r = colors[It_Threshold].r * 10.0;
  	colors[It_Threshold].g = colors[It_Threshold].g * 10.0;
  	colors[It_Threshold].b = colors[It_Threshold].b * 1000000.0;
  	// to avoid the unexpected black area
	// Or the first randomly generated RGB will be too small (looks like black)
	
	// cout << colors[It_Threshold].r << ',' << colors[It_Threshold].g << ',' << colors[It_Threshold].b << endl;

}



// GLfloat Vertex[WindowSize * WindowSize][2];
// GLfloat Color_vertex[WindowSize * WindowSize][3];

// void Update_VertexAndColor(int x, int y, int num, int It)
// {
// 	Vertex[num][0] = x;
// 	Vertex[num][1] = y;
// 	Color_vertex[num][0] = colors[It].r;
// 	Color_vertex[num][1] = colors[It].g;
// 	Color_vertex[num][2] = colors[It].b;
// }

void Update_Steps()		// Update the length of steps after using the selection window
{
	Real_Step = (maxC.real - minC.real) / (WindowSize);
	Imag_Step = (maxC.imag - minC.imag) / (WindowSize);
}

int In_MBset(Complex c)		// Judge whether c is in the MBSet
{
	Complex Z(c.real, c.imag);
	int i = 0; 		// Iteration count

	// cout << "Z1=(" << Z.real << "," << Z.imag << ")" << endl;

	while(i < maxIt)
	{
		if( Z.Mag2() <= (2.0 * 2.0) )		// Z.Mag() returns a Complex type result, hard to be compared
		{
			Z = (Z * Z) + c;
			i++;
		}
		else
			// break;
			return i;
	}

	// Test only
	// cout << "i:" << i << endl;
	// cout << "Z2=(" << Z.real << "," << Z.imag << ")" << endl;

	return i;	// Return the iteration count. i = 2000 means the Z is in the MBSet
}

void* Compute_MBSet_nThreads(void* v)
// void Compute_MBSet_nThreads(int v)
{	
	unsigned long myID = (unsigned long)v;

	pthread_mutex_lock(&ActiveCountMutex);
	cout << "Thread ID = " << myID << " is working!" << endl;
	pthread_mutex_unlock(&ActiveCountMutex);

	unsigned long start_index = Real_per_thread * Imag_per_thread * myID;

	double Zr = minC.real + myID * Real_per_thread * Real_Step;
	double Zi = minC.imag;
	Complex c(Zr, Zi);	// Z0, the begining point

	int current_location,It = 0;

	// Real: rows; Imag: columns
	for (int i = 0; i < Real_per_thread; ++i)		// x coordinate
	{
		for (int j = 0; j < Imag_per_thread; ++j)	// y coordinate
		{
			current_location = start_index + WindowSize * i + j;
			// Iteration_Count[current_location] = In_MBset(c);

			It = In_MBset(c);
			Iteration_Count[current_location] = It;

			// Update_VertexAndColor(i, j, current_location, It);

			c.imag = c.imag + Imag_Step;
		}
		c.real = c.real + Real_Step;
		c.imag = Zi;	// Reset the imaginary section
	}

	pthread_mutex_lock(&ActiveCountMutex); //Judge whether to send the signal
  	ActiveCount--;
  	printf("Threads Remain: %d\n",ActiveCount); //Test, how many active threads now?
  	if (ActiveCount == 0){
    	pthread_mutex_unlock(&ActiveCountMutex);
    	pthread_mutex_lock(&exitMutex);  //Use this lock to confirm the main thread is ready to receive signal, or it will block
    	pthread_cond_signal(&exitCond);
    	pthread_mutex_unlock(&exitMutex);
  	}
  	else{
    	pthread_mutex_unlock(&ActiveCountMutex);
  	}

}



// void Draw_MBSet_Array()	// Use OpenGL's vertex array, try to get better performance
// {	// Temporarily unavailable
// 	glVertexPointer(2, GL_FLOAT, 2, Vertex);
// 	glColorPointer(3, GL_FLOAT, 3, Color_vertex);
// 	glDrawArrays(GL_POINTS, 0, WindowSize*WindowSize);
// }

void Draw_MBSet()
{	
	int color_index;

	for(int i = 0; i < WindowSize; i++)
  	{
    	for(int j = 0; j < WindowSize; j++)
    	{
    		color_index = Iteration_Count[i*WindowSize + j];
      		glColor3f(colors[color_index].r, colors[color_index].g, colors[color_index].b);
      		glVertex2d(i, j);
      		// if (glGetError()!=0)
      			// std::cout << glGetError() << std::endl;
    	}
 	 }
}

class Position
	{
	public:
  		Position();
  		Position(int x, int y);
	public:
		int x;
		int y;
	};
Position::Position()
	: x(0), y(0)
{
}
Position::Position(int x0, int y0)
    : x(x0), y(y0)
{
}


void Creat_nThreads()
{
	pthread_mutex_lock(&exitMutex);
  	ActiveCount = nThreads; // Initialize the ActiveCount
  	for (int i = 0; i < nThreads; ++i)
    {
      	pthread_t pt; // pThread variable (output param from create)
      	// Third param is the thread starting function
      	// Fourth param is passed to the thread starting function
      	pthread_create(&pt, 0, Compute_MBSet_nThreads, (void*)i);
    }
  	pthread_cond_wait(&exitCond, &exitMutex); // Wait for all threads complete
  	pthread_mutex_unlock(&exitMutex);
}


int Enlarge_count = 0;
Complex Recovery_max[Max_Depth];
Complex Recovery_min[Max_Depth];

Position start, end, new_start_position, new_end_position;
bool draw_square = false;

void mouse(int button, int state, int x, int y)
{ // Your mouse click processing here
  // state == 0 means pressed, state != 0 means released
  // Note that the x and y coordinates passed in are in
  // PIXELS, with y = 0 at the top.

	int width, height;
	bool Clicked = false;

	if ( state == 0 && button == GLUT_LEFT_BUTTON)	
	// Press the left button, choose a new start point
	{
		start = Position(x, y);
		cout << "Clicked" << endl;

		draw_square = true;
		// Push Max, Min into the stack here
		Recovery_max[Enlarge_count].real = maxC.real;
		Recovery_max[Enlarge_count].imag = maxC.imag;
		Recovery_min[Enlarge_count].real = minC.real;
		Recovery_min[Enlarge_count].imag = minC.imag;
		Enlarge_count++;

	}

	if ( state == 1 && button == GLUT_LEFT_BUTTON)	
	// Release the left button, choose a new end point
	{
		draw_square = false;

		end = Position(x, y);
		new_start_position = Position(min(start.x, end.x), min(start.y, end.y));
		new_end_position = Position(max(start.x, end.x), max(start.y, end.y));

		width = new_end_position.x - new_start_position.x;
		height = new_end_position.y - new_start_position.y;

		if (width !=0 && height !=0)
		{
		if (width != height){
			new_end_position.x = new_start_position.x + min(width, height);
			new_end_position.y = new_start_position.y + min(width, height);
		}

		minC.real += new_start_position.x * Real_Step;
		minC.imag += new_start_position.y * Imag_Step;
		maxC.real -= (WindowSize - new_end_position.x) * Real_Step;
		maxC.imag -= (WindowSize - new_end_position.y) * Imag_Step;
		Update_Steps();	
		}
		else
		cout << "Selecting one pixel makes no sense." << endl;	

		// count = 0;
		cout << "Released" << endl;		
	}

	if ( state == 1 && button == GLUT_LEFT_BUTTON)
		// Update the figure when releasing buttion
	{
		// for (int i = 0; i < nThreads; ++i)
  			// Compute_MBSet_nThreads(i);	// Need to be implemented as pthreads
		Creat_nThreads();

		glutPostRedisplay();

	}
}

void Draw_Selected_Region(int x, int y)
{	
	int x0, x1, y0, y1;	//(x0,y0) , (x1, y0), (x1, y1), (x0, y1)
	int width, height, size;

	x0 = start.x;
	y0 = start.y;
	width = abs(x - x0);
	height = abs(y - y0);

	size = min(width, height);
	if (x > x0)
		x1 = x0 + size;
	else 
		x1 = x0 - size;
	if (y > y0)
		y1 = y0 + size;
	else 
		y1 = y0 - size;

	glBegin(GL_LINE_LOOP);
	glColor3f(1.0, 1.0, 1.0);	// White lines
  	glLineWidth(4.0);
	glVertex2d(x0,y0);
	glVertex2d(x1, y0);
	glVertex2d(x1, y1);
	glVertex2d(x0, y1);
	glEnd();

}

int count = 0;
int x_drag, y_drag;

void motion(int x, int y)
{ // Your mouse motion here, x and y coordinates are as above
	
  	// Draw_Selected_Region(x, y);	
  	x_drag = x;
  	y_drag = y;
 	// cout << "Draging: " << count++ << endl;
 	glutPostRedisplay();
}

void keyboard(unsigned char c, int x, int y)
{ // Your keyboard processing here
	// if (c == 's')
	// 	cout << "Runtime" << runtime << endl;

	if (c == 'q')
	{	cout << "Quit!" << endl;
		exit(0);
	}

	if (c == 'b')	// Go back to upper level
	{
		// Pop last Max, Min into the stack here
		if (Enlarge_count == 0)
		{
			cout << "Already in the first layer!" << endl;
		}
		else
		{

		Enlarge_count--;
		maxC.real = Recovery_max[Enlarge_count].real;
		maxC.imag = Recovery_max[Enlarge_count].imag;
		minC.real = Recovery_min[Enlarge_count].real;
		minC.imag = Recovery_min[Enlarge_count].imag;
		Update_Steps();

		// for (int i = 0; i < nThreads; ++i)
  			// Compute_MBSet_nThreads(i);	// Need to be implemented as pthreads
		Creat_nThreads();

		glutPostRedisplay();
		cout << "Go back successfully! Current depth:" << Enlarge_count << " / " << Max_Depth << endl;
		
		}
		
	}
}

void display(void)
{ // Your OpenGL display code here

	glClearColor(0.0, 0.0, 0.0, 0.0); 
  	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  	glLoadIdentity();
  	// glTranslatef(0, 64, 0);
  
  	glBegin(GL_POINTS);
  	Draw_MBSet();
 	glEnd();

 	// Draw_MBSet_Array();

 	if (draw_square)
  		Draw_Selected_Region(x_drag, y_drag);

 	glutSwapBuffers();
}

void init()
{ // Your OpenGL initialization code here

	//select clearing (background) color
	// White
	// glClearColor(1.0, 1.0, 1.0, 1.0);
	glClearColor(0.0, 0.0, 0.0, 0.0); 

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);

  	// glShadeModel(GL_FLAT);
}

void reshape(int w, int h)
{ // Your OpenGL window reshape code here
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);                                            
  	glMatrixMode(GL_PROJECTION); 
  	glLoadIdentity();

  	gluOrtho2D(0, (GLdouble)w, (GLdouble)h, 0);
  	glMatrixMode(GL_MODELVIEW);
  	glLoadIdentity(); 
}

// }

void timer(int)
{
  glutPostRedisplay();
  glutTimerFunc(1000.0 / updateRate, timer, 0);
}


void Init_pThreads()
{
	pthread_mutex_init(&exitMutex, 0);
	pthread_cond_init(&exitCond, 0);
	pthread_mutex_init(&ActiveCountMutex, 0);
	ActiveCount = 0;
}


int main(int argc, char** argv)
{
  // Initialize OpenGL, but only on the "master" thread or process.
  // See the assignment writeup to determine which is "master" 
  // and which is slave.
	glutInit(&argc, argv);
  	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  	glutInitWindowSize(WindowSize, WindowSize);
  	glutInitWindowPosition(100, 100);
  	glutCreateWindow("MBSet");
  	init();

  	Init_Colors();
  	Update_Steps();

  	Init_pThreads();
  	// for (int i = 0; i < nThreads; ++i)
  	// 	Compute_MBSet_nThreads(i);	// Need to be implemented as pthreads
  	Creat_nThreads();


  	glutDisplayFunc(display);
  	glutReshapeFunc(reshape);

  	// glutTimerFunc(1000.0 / updateRate, timer, 0);

  	glutKeyboardFunc(keyboard);  
  	glutMouseFunc(mouse);
  	glutMotionFunc(motion);
  	glutMainLoop();
  return 0;
}

