Projekt vypracovala skupina Boháč, Bubeník, Koudelka. Projekt řeší problém stacionárního 3D vedení tepla v pravidelném homogenním tělese s vnitřním tepelným zdrojem.
Cílem je numericky vyřešit Poissonovu rovnici:
<img width="844" height="101" alt="image" src="https://github.com/user-attachments/assets/62f78bc7-c4a9-4b20-b27c-d222dbd7d2b8" />

Řešení je realizováno pomocí Gauss–Seidelovy iterační metody na pravidelné 3D mřížce.

**Struktura kódu**
1) GEOMETRICKÉ A FYZIKÁLNÍ PARAMETRY
struct parametry {

    // Geometrie
    double Lx = 0.1, Ly = 0.06, Lz = 0.01;
    int Nx = 200, Ny = 100, Nz = 50;

    // Odvozené geometrické parametry
    double dx, dy, dz;

    // Fyzikální parametry
    double lambda = 0.4;    // [W/mK]
    double rho = 1200.0;    // [kg/m3]
    double cp = 1000.0;     // [J/kgK]
    double Tinf = 293.15;   // [K]
    double Q0 = 5e4;        // [W/m3]

    // Konvekční koeficienty
    double h_x = 2.0;
    double h_y = 1.0;
    double h_z = 3.0;

    // Konstruktor
    parametry() {
        dx = Lx / (Nx - 1);
        dy = Ly / (Ny - 1);
        dz = Lz / (Nz - 1);
    }
};
