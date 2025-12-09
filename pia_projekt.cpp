#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

struct parametry {

    // Geometrie
    double Lx = 0.1, Ly = 0.06, Lz = 0.01;
    int Nx = 100, Ny = 50, Nz = 25;

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
    // Rozdílné konvekční koeficienty reprezentují polohu tělesa v prostoru  

    // Konstruktor
    parametry() {
        dx = Lx / (Nx - 1);
        dy = Ly / (Ny - 1);
        dz = Lz / (Nz - 1);
    }
};

class HeatSolver {
private:
    parametry p;
    std::vector<double> T;         // teplota
    std::vector<double> Q;         // zdroj
    std::vector<double> residualHistory;

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

    // ===============================
    //         Hlavní výpočet
    // ===============================
    void solve(int maxIter, double tol) {
        std::cout << "Zacinam vypocet..." << std::endl;

        for (int iter = 0; iter < maxIter; iter++) {
            double maxChange = 0.0;

            // -----------------------------
            //   Aktualizace vnitřních bodů
            // -----------------------------
            for (int k = 1; k < p.Nz - 1; k++) {
                for (int j = 1; j < p.Ny - 1; j++) {
                    for (int i = 1; i < p.Nx - 1; i++) {

                        int id = idx(i, j, k);

                        double Tnew =
                            (T[idx(i + 1, j, k)] + T[idx(i - 1, j, k)]) / (p.dx * p.dx) +
                            (T[idx(i, j + 1, k)] + T[idx(i, j - 1, k)]) / (p.dy * p.dy) +
                            (T[idx(i, j, k + 1)] + T[idx(i, j, k - 1)]) / (p.dz * p.dz) +
                            Q[id] / p.lambda;

                        double denom = 2.0 / (p.dx * p.dx) +
                                       2.0 / (p.dy * p.dy) +
                                       2.0 / (p.dz * p.dz);

                        Tnew /= denom;

                        double diff = std::abs(Tnew - T[id]);
                        if (diff > maxChange)
                            maxChange = diff;

                        T[id] = Tnew;
                    }
                }
            }

            // Okraje
            okrajovapodm();

            // -----------------------------
            //         Výpočet rezidua
            // -----------------------------
            double maxRes = 0.0;

            for (int k = 1; k < p.Nz - 1; k++) {
                for (int j = 1; j < p.Ny - 1; j++) {
                    for (int i = 1; i < p.Nx - 1; i++) {
                        int id = idx(i, j, k);

                        double lap =
                            (T[idx(i + 1, j, k)] - 2*T[id] + T[idx(i - 1, j, k)]) / (p.dx * p.dx) +
                            (T[idx(i, j + 1, k)] - 2*T[id] + T[idx(i, j - 1, k)]) / (p.dy * p.dy) +
                            (T[idx(i, j, k + 1)] - 2*T[id] + T[idx(i, j, k - 1)]) / (p.dz * p.dz);

                        double R = p.lambda * lap + Q[id];

                        if (std::abs(R) > maxRes)
                            maxRes = std::abs(R);
                    }
                }
            }

            residualHistory.push_back(maxRes);

            // Výpis každých 100 iterací
            if (iter % 100 == 0) {
                std::cout << "Iterace " << iter
                          << " | maxChange = " << maxChange
                          << " | reziduum = " << maxRes << std::endl;
            }

                if (iter > 10 && maxChange < tol) {
                    std::cout << "Konvergence po " << iter << " iteracich." << std::endl;
                    break;
            }
        }
    }

    // ===============================
    //    Okrajové podmínky Robin
    // ===============================
    void okrajovapodm() {
        // X směr
        for (int j = 0; j < p.Ny; j++) {
            for (int k = 0; k < p.Nz; k++) {
                T[idx(0, j, k)] =
                    (p.lambda / p.dx * T[idx(1, j, k)] + p.h_x * p.Tinf) /
                    (p.lambda / p.dx + p.h_x);

                T[idx(p.Nx - 1, j, k)] =
                    (p.lambda / p.dx * T[idx(p.Nx - 2, j, k)] + p.h_x * p.Tinf) /
                    (p.lambda / p.dx + p.h_x);
            }
        }

        // Y směr
        for (int i = 0; i < p.Nx; i++) {
            for (int k = 0; k < p.Nz; k++) {
                T[idx(i, 0, k)] =
                    (p.lambda / p.dy * T[idx(i, 1, k)] + p.h_y * p.Tinf) /
                    (p.lambda / p.dy + p.h_y);

                T[idx(i, p.Ny - 1, k)] =
                    (p.lambda / p.dy * T[idx(i, p.Ny - 2, k)] + p.h_y * p.Tinf) /
                    (p.lambda / p.dy + p.h_y);
            }
        }

        // Z směr
        for (int i = 0; i < p.Nx; i++) {
            for (int j = 0; j < p.Ny; j++) {
                T[idx(i, j, 0)] =
                    (p.lambda / p.dz * T[idx(i, j, 1)] + p.h_z * p.Tinf) /
                    (p.lambda / p.dz + p.h_z);

                T[idx(i, j, p.Nz - 1)] =
                    (p.lambda / p.dz * T[idx(i, j, p.Nz - 2)] + p.h_z * p.Tinf) /
                    (p.lambda / p.dz + p.h_z);
            }
        }
    }

    // ===============================
    //    Uložení teploty do VTK
    // ===============================
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

    // ===============================
    //    Uložení reziduí do CSV
    // ===============================
    void saveResiduals(const std::string& filename) const {
        std::ofstream f(filename);
        for (size_t i = 0; i < residualHistory.size(); i++)
            f << i << "," << residualHistory[i] << "\n";

        std::cout << "Rezidua ulozena do " << filename << std::endl;
    }
};



// ==========================================
//                  MAIN
// ==========================================

int main() {
    parametry params;

    HeatSolver solver(params);

    solver.solve(10000, 1e-6);

    solver.saveVTK("temperature.vtk");
    solver.saveResiduals("residuals.csv");

    return 0;
}
