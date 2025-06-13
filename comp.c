#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#define  ORBIT_LENGTH 300000

const double G = 6.6743e-11;
const int WIDTH = 1000; 
const int HEIGHT = 800;
const double SCALE = 1e-10;

double h = 10000;
float camera_distance = 5e11 * SCALE;  
float camera_speed = 10;           

static SDL_Window* window;
static SDL_GLContext context;

typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    double m;
    float color[3];
} Planet;

typedef struct {
    double x[ORBIT_LENGTH];
    double y[ORBIT_LENGTH];
    double z[ORBIT_LENGTH];
    int count;
    int index;
} OrbitHistory;

Planet* k1;
Planet* k2;
Planet* k3;
Planet* k4;
Planet* temp;

void f(int n, Planet* planets) {
    for (int i = 0; i < n; i++) {
        planets[i].ax = 0;
        planets[i].ay = 0;
        planets[i].az = 0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double dx = (planets[j].x - planets[i].x);
                double dy = (planets[j].y - planets[i].y);
                double dz = (planets[j].z - planets[i].z);

                double r = sqrt(dx*dx + dy*dy + dz*dz);
                if(r == 0) continue;

                double a = (G * planets[j].m) / (r * r * r);
                planets[i].ax += a * dx;
                planets[i].ay += a * dy;
                planets[i].az += a * dz;
            }
        }
    }    
}


void rungekutta(int n, Planet* planets) {

    for (int i = 0; i < n; i++) temp[i] = planets[i];

    //  k1
    f(n, temp);
    for (int i = 0; i < n; i++) k1[i] = temp[i];

    // k2
    for (int i = 0; i < n; i++) {
        temp[i].x = planets[i].x + h * 0.5 * k1[i].vx;
        temp[i].vx = planets[i].vx + h * 0.5 * k1[i].ax;
        temp[i].y = planets[i].y + h * 0.5 * k1[i].vy;
        temp[i].vy = planets[i].vy + h * 0.5 * k1[i].ay;
        temp[i].z = planets[i].z + h * 0.5 * k1[i].vz;
        temp[i].vz = planets[i].vz + h * 0.5 * k1[i].az;
    }
    f(n, temp);
    for (int i = 0; i < n; i++) k2[i] = temp[i];

    // k3
    for (int i = 0; i < n; i++) {
        temp[i].x = planets[i].x + h * 0.5 * k2[i].vx;
        temp[i].vx = planets[i].vx + h * 0.5 * k2[i].ax;
        temp[i].y = planets[i].y + h * 0.5 * k2[i].vy;
        temp[i].vy = planets[i].vy + h * 0.5 * k2[i].ay;
        temp[i].z = planets[i].z + h * 0.5 * k2[i].vz;
        temp[i].vz = planets[i].vz + h * 0.5 * k2[i].az;
    }
    f(n, temp);
    for (int i = 0; i < n; i++) k3[i] = temp[i];

    // k4
    for (int i = 0; i < n; i++) {
        temp[i].x = planets[i].x + h * k3[i].vx;
        temp[i].vx = planets[i].vx + h * k3[i].ax;
        temp[i].y = planets[i].y + h * k3[i].vy;
        temp[i].vy = planets[i].vy + h * k3[i].ay;
        temp[i].z = planets[i].z + h * k3[i].vz;
        temp[i].vz = planets[i].vz + h * k3[i].az;
    }
    f(n, temp);
    for (int i = 0; i < n; i++) k4[i] = temp[i];

    for (int i = 0; i < n; i++) {
        planets[i].x += h / 6.0 * (k1[i].vx + 2*k2[i].vx + 2*k3[i].vx + k4[i].vx);
        planets[i].vx += h / 6.0 * (k1[i].ax + 2*k2[i].ax + 2*k3[i].ax + k4[i].ax);
        planets[i].y += h / 6.0 * (k1[i].vy + 2*k2[i].vy + 2*k3[i].vy + k4[i].vy);
        planets[i].vy += h / 6.0 * (k1[i].ay + 2*k2[i].ay + 2*k3[i].ay + k4[i].ay);
        planets[i].z += h / 6.0 * (k1[i].vz + 2*k2[i].vz + 2*k3[i].vz + k4[i].vz);
        planets[i].vz += h / 6.0 * (k1[i].az + 2*k2[i].az + 2*k3[i].az + k4[i].az);
    }

}


void initPlanet(Planet *p, double x, double y, double z, double vx, double vy, double vz, double m, float r, float g, float b) {
    p->x = x; p->y = y; p->z = z;
    p->vx = vx; p->vy = vy; p->vz = vz;
    p->ax = 0.0; p->ay = 0.0; p->az = 0.0; 
    p->m = m; 
    p->color[0] = r; p->color[1] = g; p->color[2] = b;
}

void initOrbits(OrbitHistory* orbits, int num_planets) {
    for (int i = 0; i < num_planets; i++) {
        orbits[i].index = 0;  
        memset(orbits[i].x, 0, ORBIT_LENGTH * sizeof(double));
        memset(orbits[i].y, 0, ORBIT_LENGTH * sizeof(double));
        memset(orbits[i].z, 0, ORBIT_LENGTH * sizeof(double));
    }
}

void calculateBarycenter(Planet* planets, int n, double* bx, double* by, double* bz) {
    double total_mass = 0;
    *bx = *by = *bz = 0;
    
    for(int i = 0; i < n; i++) {
        total_mass += planets[i].m;
        *bx += planets[i].m * planets[i].x;
        *by += planets[i].m * planets[i].y;
        *bz += planets[i].m * planets[i].z;
    }
    
    *bx /= total_mass;
    *by /= total_mass;
    *bz /= total_mass;
}

void calculateBarycenterVelocity(Planet* planets, int n, double* bvx, double* bvy, double* bvz) {
    double total_mass = 0;
    *bvx = *bvy = *bvz = 0;
    
    for(int i = 0; i < n; i++) {
        total_mass += planets[i].m;
        *bvx += planets[i].m * planets[i].vx;
        *bvy += planets[i].m * planets[i].vy;
        *bvz += planets[i].m * planets[i].vz;
    }
    
    *bvx /= total_mass;
    *bvy /= total_mass;
    *bvz /= total_mass;
}

void adjustToBarycenter(Planet* planets, int n) {
    double bx, by, bz;
    double bvx, bvy, bvz;
    
    calculateBarycenter(planets, n, &bx, &by, &bz);
    calculateBarycenterVelocity(planets, n, &bvx, &bvy, &bvz);
    
    for(int i = 0; i < n; i++) {
        planets[i].x -= bx;
        planets[i].y -= by;
        planets[i].z -= bz;
        planets[i].vx -= bvx;
        planets[i].vy -= bvy;
        planets[i].vz -= bvz;
    }
}

void calculateAngularMomentum(Planet* planets, int n, double* Lx, double* Ly, double* Lz) {
    *Lx = *Ly = *Lz = 0;
    
    for(int i = 0; i < n; i++) {
        *Lx += planets[i].m * (planets[i].y * planets[i].vz - planets[i].z * planets[i].vy);
        *Ly += planets[i].m * (planets[i].z * planets[i].vx - planets[i].x * planets[i].vz);
        *Lz += planets[i].m * (planets[i].x * planets[i].vy - planets[i].y * planets[i].vx);
    }
}

double calculateTotalEnergy(Planet* planets, int n) {
    double kinetic = 0;
    double potential = 0;
    
    // Кинетическая энергия
    for(int i = 0; i < n; i++) {
        double v2 = planets[i].vx*planets[i].vx + 
                    planets[i].vy*planets[i].vy + 
                    planets[i].vz*planets[i].vz;
        kinetic += 0.5 * planets[i].m * v2;
    }
    
    // Потенциальная энергия
    for(int i = 0; i < n; i++) {
        for(int j = i+1; j < n; j++) {
            double dx = planets[j].x - planets[i].x;
            double dy = planets[j].y - planets[i].y;
            double dz = planets[j].z - planets[i].z;
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            potential -= G * planets[i].m * planets[j].m / r;
        }
    }
    
    return kinetic + potential;
}

void logSystemState(int step, int num_planets, Planet* planets, double initial_energy) {
    double bx, by, bz;
    calculateBarycenter(planets, num_planets, &bx, &by, &bz);
    
    double bvx, bvy, bvz;
    calculateBarycenterVelocity(planets, num_planets, &bvx, &bvy, &bvz);
    
    double Lx, Ly, Lz;
    calculateAngularMomentum(planets, num_planets, &Lx, &Ly, &Lz);
    
    double energy = calculateTotalEnergy(planets, num_planets);
    
    printf("\nStep: %d\n", step);
    printf("Barycenter position: (%e, %e, %e) m\n", bx, by, bz);
    printf("Barycenter velocity: (%e, %e, %e) m/s\n", bvx, bvy, bvz);
    printf("Angular momentum: (%e, %e, %e) kg·m²/s\n", Lx, Ly, Lz);
    printf("Total energy: %e J\n", energy);
    printf("Relative energy error: %e %%\n", fabs((energy - initial_energy) / initial_energy * 100));
}

void addOrbitPoint(OrbitHistory* orbit, const Planet* planet) {
    // Сохраняем текущую позицию в массив по текущему индексу
    
    orbit->x[orbit->index] = planet->x;
    orbit->y[orbit->index] = planet->y;
    orbit->z[orbit->index] =planet->z;

    orbit->index = (orbit->index + 1) % ORBIT_LENGTH;
}
 
void drawOrbit(const OrbitHistory* orbit, float color[3]) {
    glDisable(GL_LIGHTING);
    glLineWidth(1.0f);
    glBegin(GL_LINE_STRIP);
    glColor3f(1,1,1);
    for(int i = 0; i < orbit->index; i++) {
        glVertex3d(orbit->x[i] * SCALE, 
                  orbit->y[i] * SCALE,
                  orbit->z[i] * SCALE);
    }
    glEnd();
    glEnable(GL_LIGHTING);
}
void drawPlanet(Planet p) {
    glPushMatrix();
    GLfloat mat_color[] = {p.color[0], p.color[1], p.color[2], 1.0f};
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_color);
    glTranslatef(p.x * SCALE, p.y * SCALE, p.z * SCALE);
    
        
    float radius = 5e9 + log(p.m/1e24) * 1e9;
    GLUquadric* quadric = gluNewQuadric();
    gluSphere(quadric, radius * SCALE, 20, 20);
    gluDeleteQuadric(quadric);
    glPopMatrix();
}

void draw(int num_planets, Planet* planets, OrbitHistory* orbits) {
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT ); 

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, (double)WIDTH/HEIGHT, 1e10*SCALE, 1e16*SCALE);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0, 0, camera_distance, 0, 0, 0, 0, 1, 0);
    
    glDisable(GL_LIGHTING);
    for (int i = 1; i < num_planets; i++) { 
        drawOrbit(&orbits[i], planets[i].color);
    }
    
    for (int i = 0; i < num_planets; i++) {
        drawPlanet(planets[i]);
    }

    SDL_GL_SwapWindow(window);
}

int get_input() {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        if (event.type == SDL_KEYDOWN) {
            switch(event.key.keysym.sym) {
                case SDLK_ESCAPE: 
                    return 0;
                case SDLK_PLUS:
                case SDLK_EQUALS: // + =
                    camera_distance -= camera_speed;
                    break;
                case SDLK_MINUS:   // -
                    camera_distance += camera_speed;
                    break;
            }
        }    }
    return 1;
}

int main(int argc, char** argv) {
    if (argc >= 2) {
        h = atof(argv[1]);
    }
    int num_planets = 10;
    Planet* planets = (Planet*)malloc(num_planets * sizeof(Planet));
    OrbitHistory* orbits = (OrbitHistory*)malloc(num_planets * sizeof(OrbitHistory));
    k1 = (Planet*)malloc(num_planets * sizeof(Planet));
    k2 = (Planet*)malloc(num_planets * sizeof(Planet));
    k3 = (Planet*)malloc(num_planets * sizeof(Planet));
    k4 = (Planet*)malloc(num_planets * sizeof(Planet));
    temp = (Planet*)malloc(num_planets * sizeof(Planet));

    SDL_Init(SDL_INIT_VIDEO);
    window = SDL_CreateWindow("Solar System", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 
                            WIDTH, HEIGHT, SDL_WINDOW_OPENGL);
    context = SDL_GL_CreateContext(window);
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);

    initOrbits(orbits, num_planets);
    // Солнце
    initPlanet(&planets[0], -7.6969e8, -7.1919e8, -2.8414e8, 12.594, -4.175, -2.043, 1.989e30, 1.0f, 0.9f, 0.0f);
    // Земля
    initPlanet(&planets[1], -1.496e11, -5.6097e9, -2.4026e9, 582, -2.7421e4, -1.1887e4, 5.972e24, 0.2f, 0.4f, 1.0f);
    // Юпитер
    initPlanet(&planets[2], 6.5912e10, 6.9977e11, 2.9834e11, -13167, 1491, 959, 1.8987e27, 0.5f, 0.5f, 0.5f);
    // Луна
    initPlanet(&planets[3], -1.4972e11, -5.9522e9, -2.5901e9, 1567, -27285, -11807, 7.36e22, 0.8f, 0.6f, 0.4f);
    // Марс
    initPlanet(&planets[4], -2.0488e11, 1.2637e11, 6.3516e10, -12913, -16333, -7143, 6.39e23, 0.7f, 0.2f, 0.1f);
    // Сатурн
    initPlanet(&planets[5], 1.4216e12, -1.6081e11, -1.2765e11, 798, 8837, 3615, 5.683e26, 0.9f, 0.8f, 0.6f);
    // Уран
    initPlanet(&planets[6], 1.6205e12, 2.2355e12, 9.5618e11, -5716, 3139, 1455, 8.6849e25, 0.4f, 0.8f, 0.8f);
    // Нептун
    initPlanet(&planets[7], 4.4469e12, -1.2220e10, -1.1627e11, 34, 5061, 2070, 1.024e26, 0.2f, 0.3f, 0.9f);
    // Венера
    initPlanet(&planets[8], -1.0827e11, -6.4209e9, 3.9515e9, 959, -32045, 14479, 4.867e24, 0.6f, 0.6f, 0.6f);
    // Меркурий
    initPlanet(&planets[9], -5.750e10, 7.6497e8, 6.3911e9, -13559, -41651, -20843, 3.285e23, 0.5f, 0.2f, 1.0f);
    
    //Корректировка относительно барицентра
    adjustToBarycenter(planets, num_planets);
    // Вычисление начальной энергии
    double initial_energy = calculateTotalEnergy(planets, num_planets);
    int cnt = 0;
    while (get_input()) {
        rungekutta(num_planets, planets);
        for (int i = 0; i < num_planets; i++) {
            addOrbitPoint(&orbits[i], &planets[i]); 
        }              
        draw(num_planets, planets, orbits);
        if (cnt % 100 == 0) logSystemState(cnt, num_planets, planets, initial_energy);
        cnt++;
    }
    
    SDL_GL_DeleteContext(context);
    SDL_DestroyWindow(window);
    SDL_Quit();
    free(planets);
    free(orbits);
    return 0;
}
