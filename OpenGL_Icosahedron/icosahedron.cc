// Draw an Icosahedron
// ECE4893/8893 Project 4
// Xueyang Xu

#include <iostream>
#include <math.h>
#include <GL/glut.h>
#include <GL/glext.h>
#include <GL/gl.h>
#include <GL/glu.h>

using namespace std;

#define NFACE 20
#define NVERTEX 12

#define X .525731112119133606 
#define Z .850650808352039932

static int updateRate = 10;

// These are the 12 vertices for the icosahedron
static GLfloat vdata[NVERTEX][3] = {    
   {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},    
   {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},    
   {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0} 
};

// These are the 20 faces.  Each of the three entries for each 
// vertex gives the 3 vertices that make the face.
static GLint tindices[NFACE][3] = { 
   {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},    
   {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},    
   {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6}, 
   {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11} };

int testNumber; // Global variable indicating which test number is desired
int depth; // Recursion depth

// // Assign similar colors to sub-triangles in the same Basic triangles
// GLfloat Color_R;
// GLfloat Color_G;
// GLfloat Color_B;
 
void drawTriangle(GLfloat v1[3], GLfloat v2[3], GLfloat v3[3])
{
  // Draw 1 triangle and 3 edges each time called
  
  glBegin(GL_LINE_LOOP);
  glColor3f(1.0, 1.0, 1.0);
  glLineWidth(2.0);
  glVertex3fv(v1);
  glVertex3fv(v2);
  glVertex3fv(v3);
  glEnd();

  // glBegin(GL_LINES);
  // glColor3f(1.0, 1.0, 1.0);
  // glLineWidth(2.0);
  // glVertex3fv(v1);
  // glVertex3fv(v2);
  // glEnd();

  // glBegin(GL_LINES);
  // glColor3f(1.0, 1.0, 1.0);
  // glLineWidth(2.0);
  // glVertex3fv(v1);
  // glVertex3fv(v3);
  // glEnd();

  // glBegin(GL_LINES);
  // glColor3f(1.0, 1.0, 1.0);
  // glLineWidth(2.0);
  // glVertex3fv(v2);
  // glVertex3fv(v3);
  // glEnd();
  
  // Assign similar colors to sub-triangles in the same Basic triangles
  // Decide colors by vertices of each triangle
  // Map -1.37~+1.37 to 0.0~+1.0 roughly
  GLfloat Color_R = (v1[0]+v2[0]+v3[0])/3.0+0.5;
  GLfloat Color_G = (v1[1]+v2[1]+v3[1])/3.0+0.5;
  GLfloat Color_B = (v1[2]+v2[2]+v3[2])/3.0+0.5;

  glBegin(GL_TRIANGLES);
  glColor3f(Color_R, Color_G, Color_B);
  glVertex3fv(v1);
  glVertex3fv(v2);
  glVertex3fv(v3);
  glEnd();
}

// Normalize the norm of vector from the origin point to each middle point to 1 unit
// Required in Test3's requirement 
// Look like a sphere when recursion depth go higer...
void Normalization(GLfloat v[3]){
  GLfloat ratio = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  v[0] = v[0] / ratio;
  v[1] = v[1] / ratio;
  v[2] = v[2] / ratio;
}

// Calculate recursive triangles and DRAW IT!
void Divide_Recursive(GLfloat v1[3], GLfloat v2[3], GLfloat v3[3], int depth)
{
  GLfloat v12[3], v13[3], v23[3]; // The middle point of each edge
  for (int i = 0; i < 3; ++i){
    v12[i] = (v1[i] + v2[i]) / 2.0;
    v13[i] = (v1[i] + v3[i]) / 2.0;
    v23[i] = (v2[i] + v3[i]) / 2.0;
  }
  Normalization(v12);
  Normalization(v13);
  Normalization(v23);

  if (depth == 0)
    drawTriangle(v1, v2, v3);
  
  if (depth != 0){ 
    // Each time the divide operation will create 4 triangles
    Divide_Recursive(v12, v13, v23, (depth-1));
    Divide_Recursive(v1, v13, v12, (depth-1));
    Divide_Recursive(v2, v12, v23, (depth-1));
    Divide_Recursive(v3, v13, v23, (depth-1));
  }

}

void init()
{
  //select clearing (background) color
  glClearColor(0.0, 0.0, 0.0, 0.0);
  //glShadeModel(GL_FLAT);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS); // As default, is not necessary to set
}

void timer(int)
{
  glutPostRedisplay();
  glutTimerFunc(1000.0 / updateRate, timer, 0);
}

// Test cases.  Fill in your code for each test case
void Test1()
{
  // for (int i = 0; i < NFACE; ++i){
  //   drawTriangle( 
  //   // Addresss of first element in every vertex's coordinate
  //     &vdata[tindices[i][0]][0], 
  //     &vdata[tindices[i][1]][0],
  //     &vdata[tindices[i][2]][0]
  //     );

  // Use general function instead of above form
   for (int i = 0; i < NFACE; ++i){
      Divide_Recursive( 
    // Addresss of first element in every vertex's coordinate
      &vdata[tindices[i][0]][0], 
      &vdata[tindices[i][1]][0],
      &vdata[tindices[i][2]][0],
      0 // Depth = 0 in this test case
      );
  }
  
}

void Test2()
{ 
  static GLfloat rotX = 0.0;
  static GLfloat rotY = 0.0;
  glRotatef(rotX, 1.0, 0.0, 0.0);
  glRotatef(rotY, 0.0, 1.0, 0.0);
  rotX += 1.0;
  rotY -= 1.0;

  Test1();

}

void Test3()
{
  for (int i = 0; i < NFACE; ++i){
    Divide_Recursive( 
    // Addresss of first element in every vertex's coordinate
      &vdata[tindices[i][0]][0], 
      &vdata[tindices[i][1]][0],
      &vdata[tindices[i][2]][0],
      1 // Depth = 1 in this test case
      );
  }
}

void Test4()
{
  static GLfloat rotX = 0.0;
  static GLfloat rotY = 0.0;
  glRotatef(rotX, 1.0, 0.0, 0.0);
  glRotatef(rotY, 0.0, 1.0, 0.0);
  rotX += 1.0;
  rotY -= 1.0;

  Test3();
}

void Test5(int depth)
{
    for (int i = 0; i < NFACE; ++i){
      Divide_Recursive( 
    // Addresss of first element in every vertex's coordinate
      &vdata[tindices[i][0]][0], 
      &vdata[tindices[i][1]][0],
      &vdata[tindices[i][2]][0],
      depth // Depth = 1 in this test case
      );
  }
}

void Test6(int depth)
{
  static GLfloat rotX = 0.0;
  static GLfloat rotY = 0.0;
  glRotatef(rotX, 1.0, 0.0, 0.0);
  glRotatef(rotY, 0.0, 1.0, 0.0);
  rotX += 1.0;
  rotY -= 1.0;

  Test5(depth);
}


void reshape(int w, int h)
{
  glViewport(0,0, (GLsizei)w, (GLsizei)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  glOrtho(0.0, (GLdouble)w, (GLdouble)0.0, h, (GLdouble)-w, (GLdouble)w);
  
  // glOrtho(-1.0, (GLdouble)1.0, (GLdouble)-1.0, 1.0, (GLdouble)-w, (GLdouble)w);
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

}

void display(void)
{
  // static int pass;
  // cout << "Displaying pass " << ++pass << endl;

  // clear all
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // Clear the matrix
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  // Set the viewing transformation
  gluLookAt(0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

  glTranslatef(250, 250, 0);
  glScalef(250.0, 250.0, 250.0);
  
  if (testNumber == 1)
    Test1();
  else if (testNumber == 2)
    Test2();
  else if (testNumber == 3)
    Test3();
  else if (testNumber == 4)
    Test4();
  else if (testNumber == 5){
    if (depth >= 5){
      cout << "Warning: Depth too large" << endl;
      exit(1);
    }
    Test5(depth);
  }  
  else if (testNumber == 6){
    if (depth >= 5){
      cout << "Warning: Depth too large" << endl;
      exit(1);
    }
    Test6(depth);
  } 
  else{
    cout << "ERROR: Wrong test number" << endl;
    exit(1);
  }

  // Flush buffer
  // glFlush(); // If single buffering
  glutSwapBuffers(); // If double buffering
}


int main(int argc, char** argv)
{
  if (argc < 2)
    {
      std::cout << "Usage: icosahedron testnumber depth" << endl;
      exit(1);
    }
  // Set the global test number
  testNumber = atol(argv[1]);
  if (argv[2] != NULL)
    // Set the recursion depth
    depth = atol(argv[2]);

  // Initialize glut  and create your window here
  // Set your glut callbacks here
  // Enter the glut main loop here
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  // glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(500, 500);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("Icosahedron");
  init();
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutTimerFunc(1000.0 / updateRate, timer, 0);
  glutMainLoop();

  return 0;
}

