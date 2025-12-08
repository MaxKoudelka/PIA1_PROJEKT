#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

     struct parametry {

// Geometrie
     double Lx = 0.1, Ly = 0.06, Lz = 0.01;
     int Nx = 200, Ny = 100, Nz = 50;

     // Odvozené geometrické parametry
     double dx, dy, dz;

// Fyzikální parametry
     double lambda = 0.4;    // [W/mK]
     double rho = 1200.0;    // [kg/m3]
     double cp = 1000.0;   // [J/kgK]
     double Tinf = 293.15;   // [K]
     double Q0 = 5e4;    // [W/m3]

// Konvekční koeficienty
     double h_x = 2.0;   // chlazení na bocích
     double h_y = 1.0;   // slabší chlazení
     double h_z = 3.0;   // silnější chlazení

// Konstruktor pro výpočet kroků sítě
     parametry() {
         dx = Lx / (Nx - 1);
         dy = Ly / (Ny - 1);
         dz = Lz / (Nz - 1);
     }
};

class HeatSolver {
private:
     parametry p;             // Uložené parametry
     std::vector<double> T;   // Teplotní pole
     std::vector<double> Q;   // Zdrojové pole

     // Pomocná funkce pro index
        int idx(int i, int j, int k) const {
        return i + p.Nx * (j + p.Ny * k);
        }

     void zdrojtepla() {
         for (int k = p.Nz / 3; k < 2 * p.Nz / 3; k++)
             for (int j = p.Ny / 3; j < 2 * p.Ny / 3; j++)
                 for (int i = p.Nx / 3; i < 2 * p.Nx / 3; i++)
                     Q[idx(i, j, k)] = p.Q0;
     }

public:
     HeatSolver(const parametry& params) : p(params) {
         T.resize(p.Nx * p.Ny * p.Nz, p.Tinf);
         Q.resize(p.Nx * p.Ny * p.Nz, 0.0);
         zdrojtepla();
     }

// --- Hlavní výpočetní smyčka ---
     void solve(int maxIter, double tol) {
         std::cout << "Zacinam vypocet..." << std::endl;

         for (int iter = 0; iter < maxIter; iter++) {
             double maxChange = 0.0;

// Vnitřní body
             for (int k = 1; k < p.Nz - 1; k++) {
                 for (int j = 1; j < p.Ny - 1; j++) {
                     for (int i = 1; i < p.Nx - 1; i++) {
                         int id = idx(i, j, k);

                         double Tnew =
                             (T[idx(i + 1, j, k)] + T[idx(i - 1, j, k)]) / (p.dx * p.dx) +
                             (T[idx(i, j + 1, k)] + T[idx(i, j - 1, k)]) / (p.dy * p.dy) +
                             (T[idx(i, j, k + 1)] + T[idx(i, j, k - 1)]) / (p.dz * p.dz) +
                             Q[id] / p.lambda;

                         double denom = 2.0 / (p.dx * p.dx) + 2.0 / (p.dy * p.dy) + 2.0 / (p.dz * p.dz);
                         Tnew /= denom;

                         double diff = std::abs(Tnew - T[id]);
                         if (diff > maxChange) maxChange = diff;

                         T[id] = Tnew;
                     }
                 }
             }

// Okrajové podmínky (Robin)
             okrajovapodm();

             if (iter % 100 == 0)
                 std::cout << "Iterace " << iter << ", maxChange=" << 
maxChange << std::endl;

             if (maxChange < tol) {
                 std::cout << "Konvergence po " << iter << "iteracich." << std::endl;
                 break;
             }
         }
     }

// --- Aplikace okrajových podmínek ---
     void okrajovapodm() {
         // Směr X
         for (int j = 0; j < p.Ny; j++) {
             for (int k = 0; k < p.Nz; k++) {
                 // Left
                 T[idx(0, j, k)] = (p.lambda / p.dx * T[idx(1, j, k)] + p.h_x * p.Tinf) / (p.lambda / p.dx + p.h_x);
                 // Right
                 T[idx(p.Nx - 1, j, k)] = (p.lambda / p.dx * T[idx(p.Nx - 2, j, k)] + p.h_x * p.Tinf) / (p.lambda / p.dx + p.h_x);
             }
         }
         // Směr Y
         for (int i = 0; i < p.Nx; i++) {
             for (int k = 0; k < p.Nz; k++) {
                 // Front
                 T[idx(i, 0, k)] = (p.lambda / p.dy * T[idx(i, 1, k)] + p.h_y * p.Tinf) / (p.lambda / p.dy + p.h_y);
                 // Back
                 T[idx(i, p.Ny - 1, k)] = (p.lambda / p.dy * T[idx(i, p.Ny - 2, k)] + p.h_y * p.Tinf) / (p.lambda / p.dy + p.h_y);
             }
         }
         // Směr Z
         for (int i = 0; i < p.Nx; i++) {
             for (int j = 0; j < p.Ny; j++) {
                 // Bottom
                 T[idx(i, j, 0)] = (p.lambda / p.dz * T[idx(i, j, 1)] + p.h_z * p.Tinf) / (p.lambda / p.dz + p.h_z);
                 // Top
                 T[idx(i, j, p.Nz - 1)] = (p.lambda / p.dz * T[idx(i, j, p.Nz - 2)] + p.h_z * p.Tinf) / (p.lambda / p.dz + p.h_z);
             }
         }
     }

// --- Uložení výsledků ---
     void saveVTK(const std::string& filename) const {
         std::ofstream vtk(filename);
         if (!vtk.is_open()) {
             std::cerr << "Chyba pri otevirani souboru!" << std::endl;
             return;
         }

         vtk << "# vtk DataFile Version 3.0\n";
         vtk << "Temperature field\n";
         vtk << "ASCII\n";
         vtk << "DATASET STRUCTURED_POINTS\n";
         vtk << "DIMENSIONS " << p.Nx << " " << p.Ny << " " << p.Nz << "\n";
         vtk << "ORIGIN 0 0 0\n";
         vtk << "SPACING " << p.dx << " " << p.dy << " " << p.dz << "\n";
         vtk << "POINT_DATA " << p.Nx * p.Ny * p.Nz << "\n";
         vtk << "SCALARS Temperature float 1\n";
         vtk << "LOOKUP_TABLE default\n";

         for (int k = 0; k < p.Nz; k++)
             for (int j = 0; j < p.Ny; j++)
                 for (int i = 0; i < p.Nx; i++)
                     vtk << T[idx(i, j, k)] << "\n";

         vtk.close();
         std::cout << "Vysledky ulozeny do " << filename << std::endl;
     }
};



int main() {
     // 1. Nastavení parametrů
     parametry params;

     // 2. Vytvoření řešiče
     HeatSolver solver(params);

     // 3. Spuštění výpočtu
     solver.solve(1000, 1e-6);

     // 4. Uložení dat
     solver.saveVTK("temperature.vtk");

     return 0;
}