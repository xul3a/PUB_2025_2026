#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;


// g++ -o2 -std=c++17 rk_presa_predador.cpp -o rk_presa_predador // ./rk_presa_predador

//parametros do modelo Lotka–Volterra
const float a = 10.0; // crescimento presas
const float b = 10.0; // predacao
const float c = 1.0;  // mortalidade predadores
const float d = 1.0;  // ganho predadores por encontro

//calcula dx/dt e dy/dt do modelo
void fLotka(float t, float x, float y,
             float &dxdt, float &dydt)
{
    dxdt = a * x - b * x * y;
    dydt = -c * y + d * x * y;
}

//constante de movimento de Lotka–Volterra:
// V(x,y) = d x - c ln x + b y - a ln y
float V_LV(float x, float y)
{
    float resultado;

    resultado = d * x - c * logf(x) + b * y - a * logf(y);

    return resultado;
}

// Faz um passo do metodo RKF45 
void rkf45_step(float t, float x, float y, float h,
                float &x5, float &y5, float &erro)
{
    float k1x, k1y;
    float k2x, k2y;
    float k3x, k3y;
    float k4x, k4y;
    float k5x, k5y;
    float k6x, k6y;

    float dxdt, dydt;
    float x_aux, y_aux;

    // k1 = h * f(t, x, y)
    fLotka(t, x, y, dxdt, dydt);
    k1x = h * dxdt;  k1y = h * dydt;

    // k2 = h * f(t + h/4, x + (1/4)*k1x, y + (1/4)*k1y)
    x_aux = x + (1.0/4.0) * k1x;  y_aux = y + (1.0/4.0) * k1y;
    fLotka(t + h/4.0, x_aux, y_aux, dxdt, dydt);
    k2x = h * dxdt;  k2y = h * dydt;

    // k3 = h * f(t + 3h/8, x + 3/32 k1x + 9/32 k2x, y + 3/32 k1y + 9/32 k2y)
    x_aux = x + (3.0/32.0)*k1x + (9.0/32.0)*k2x;
    y_aux = y + (3.0/32.0)*k1y + (9.0/32.0)*k2y;
    fLotka(t + 3.0*h/8.0, x_aux, y_aux, dxdt, dydt);
    k3x = h * dxdt;  k3y = h * dydt;

    x_aux = x + (1932.0/2197.0)*k1x - (7200.0/2197.0)*k2x + (7296.0/2197.0)*k3x;
    y_aux = y + (1932.0/2197.0)*k1y - (7200.0/2197.0)*k2y + (7296.0/2197.0)*k3y;
    fLotka(t + 12.0*h/13.0, x_aux, y_aux, dxdt, dydt);
    k4x = h * dxdt;  k4y = h * dydt;

    x_aux = x + (439.0/216.0)*k1x - 8.0*k2x + (3680.0/513.0)*k3x - (845.0/4104.0)*k4x;
    y_aux = y + (439.0/216.0)*k1y - 8.0*k2y + (3680.0/513.0)*k3y - (845.0/4104.0)*k4y;
    fLotka(t + h, x_aux, y_aux, dxdt, dydt);
    k5x = h * dxdt;  k5y = h * dydt;


    x_aux = x - (8.0/27.0)*k1x + 2.0*k2x - (3544.0/2565.0)*k3x + (1859.0/4104.0)*k4x - (11.0/40.0)*k5x;
    y_aux = y - (8.0/27.0)*k1y + 2.0*k2y - (3544.0/2565.0)*k3y + (1859.0/4104.0)*k4y - (11.0/40.0)*k5y;
    fLotka(t + h/2.0, x_aux, y_aux, dxdt, dydt);
    k6x = h * dxdt;  k6y = h * dydt;

    // y4 (ordem 4)
    float x4, y4;
    x4 = x + (25.0/216.0)*k1x + (1408.0/2565.0)*k3x + (2197.0/4104.0)*k4x - (1.0/5.0)*k5x;
    y4 = y + (25.0/216.0)*k1y + (1408.0/2565.0)*k3y + (2197.0/4104.0)*k4y - (1.0/5.0)*k5y;

    // y5 (ordem 5)
    x5 = x + (16.0/135.0)*k1x + (6656.0/12825.0)*k3x + (28561.0/56430.0)*k4x - (9.0/50.0)*k5x + (2.0/55.0)*k6x;
    y5 = y + (16.0/135.0)*k1y + (6656.0/12825.0)*k3y + (28561.0/56430.0)*k4y - (9.0/50.0)*k5y + (2.0/55.0)*k6y;

    // erro = max(|x5 - x4|, |y5 - y4|)
    float erro_x = x5 - x4;
    float erro_y = y5 - y4;

    float abs_erro_x = fabsf(erro_x);
    float abs_erro_y = fabsf(erro_y);

    if (abs_erro_x > abs_erro_y) {
        erro = abs_erro_x;
    } else {
        erro = abs_erro_y;
    }
}
void integrar_LV(float t0, float tf,
                 float x0, float y0,
                 float h0, float tol)
{
    float t = t0;
    float x = x0;
    float y = y0;
    float h = h0;

    const float SAFETY = 0.9;
    const float TINY   = 1e-30;

    float K0 = V_LV(x0, y0);
    cout << fixed << setprecision(10);
    ofstream arquivo("t60.tsv");
    arquivo << fixed << setprecision(10);

    cout    << "# K0 = " << K0 << "\n";
    arquivo << "# K0 = " << K0 << "\n";

    float erro_relativo = 0.0;

    // cabeçalho
    cout    << "# t\tx(t)\ty(t)\th\tK(t)\terro_relativo\n";
    arquivo << "# t\tx(t)\ty(t)\th\tK(t)\terro_relativo\n";


    while (t < tf) {

        // ajusta o passo para nao ultrapassar tf
        if (t + h > tf) {
            h = tf - t;
        }

        float x5, y5;
        float erro_local;

        rkf45_step(t, x, y, h, x5, y5, erro_local);

        //calcula yscal_x e yscal_y
        float dxdt, dydt;
        fLotka(t, x, y, dxdt, dydt);

        float yscal_x = fabsf(x) + fabsf(dxdt) * h + TINY;
        float yscal_y = fabsf(y) + fabsf(dydt) * h + TINY;

        float yscal;
        if (yscal_x > yscal_y) {
            yscal = yscal_x;
        } else {
            yscal = yscal_y;
        }

        //testa se o erro_local é aceitavel
        if (erro_local > tol * yscal) {

            //rejeita o passo: diminui h e tenta de novo
            float razao = (tol * yscal) / (erro_local + TINY);
            float S = SAFETY * powf(razao, 0.25);

            if (S < 0.1) S = 0.1;
            if (S > 4.0) S = 4.0;

            h = h * S;
            continue;   
        }

        // aceita o passo:
        t = t + h;
        x = x5;
        y = y5;

        // calcula K(t) e o erro relativo da constante
        float K_atual = V_LV(x, y);
        erro_relativo = fabsf(K_atual - K0) / fabsf(K0);

        cout << t << "\t" << x << "\t" << y
            << "\t" << h << "\t" << K_atual
            << "\t" << erro_relativo << "\n";

        arquivo << t << "\t" << x << "\t" << y
                << "\t" << h << "\t" << K_atual
                << "\t" << erro_relativo << "\n";


         //justa h para o proximo passo 
        float razao2 = (tol * yscal) / (erro_local + TINY);
        float S2 = SAFETY * powf(razao2, 0.2);

        if (S2 < 0.1) S2 = 0.1;
        if (S2 > 4.0) S2 = 4.0;

        h = h * S2;
    }
}
int main()
{
    float t0  = 0.0;    // tempo 
    float tf  = 60.0;   // tempo 
    float h0  = 0.05;   // passo 
    float tol = 1e-6;   // tolerancia 

    float x0 = 3.0;     // presa 
    float y0 = 1.0;     // predador 

    integrar_LV(t0, tf, x0, y0, h0, tol);

    return 0;
}
