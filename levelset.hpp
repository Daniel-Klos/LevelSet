#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>

struct Vector2 {
    float x;
    float y;
};

struct LevelSet {
    int nX;
    int nY;
    int gridSize;
    float cellSpacing;
    float halfSpacing;

    std::vector<float> phi;

    std::vector<int> grid;

    std::vector<sf::Vector2f> dphi;

    int AIR = 0;
    int SOLID = 1;

    sf::RectangleShape rect;
    sf::RectangleShape tinyRect;
    float stepSize = 2.f;

    sf::RenderWindow &window;

    std::mt19937 gen;
    std::bernoulli_distribution distrib;

    sf::Text text;

    sf::CircleShape circle;

    float WIDTH;
    float HEIGHT;

    LevelSet(sf::Text text_, sf::RenderWindow &window_, int nX_, float WIDTH_, float HEIGHT_): WIDTH(WIDTH_), HEIGHT(HEIGHT_), text(text_), window(window_), nX(nX_), gen(std::random_device{}()), distrib(0.5) {
        cellSpacing = WIDTH / nX;
        halfSpacing = cellSpacing * 0.5;

        rect.setFillColor(sf::Color::Transparent);
        rect.setSize(sf::Vector2f{cellSpacing, cellSpacing});
        rect.setOutlineThickness(1.f);
        rect.setOutlineColor(sf::Color::White);

        tinyRect.setSize(sf::Vector2f{stepSize, stepSize});

        float circleRad = 5.f;
        circle.setFillColor(sf::Color::Transparent);
        circle.setOutlineColor(sf::Color::Red);
        circle.setOutlineThickness(5.f);
        circle.setRadius(circleRad);
        circle.setOrigin(circleRad * 0.5f, circleRad * 0.5f);
        
        nY = std::floor(HEIGHT / cellSpacing);

        gridSize = nX * nY;

        grid.resize(gridSize);
        phi.resize(gridSize);

        /*for (int i = 0; i < gridSize; ++i) {
            grid[i] = distrib(gen);
        }*/

        /*int center = (nX / 2) * nY + (nY / 2);
        grid[center] = 1;*/

        float R = 120;   // outer radius
        float r = 80;    // inner radius

        for (int ix = 0; ix < nX; ++ix) {
            for (int iy = 0; iy < nY; ++iy) {
                int idx = ix * nY + iy;
            
                float dx = ix - nX / 2;
                float dy = iy - nY / 2;
                float d2 = dx*dx + dy*dy;
            
                grid[idx] = (d2 <= R*R && d2 >= r*r) ? 1 : 0;
            }
        }

        FastSweep();

        //printPhiAndGrid();
    }

    sf::Vector2f gridCellToPos(int idx) {
        int i = idx % nY;
        int j = idx / nY;
        float x = j * cellSpacing;
        float y = i * cellSpacing;
        return sf::Vector2f{x, y};
    }

    void DrawSolidCellAt(int cellID) {
        sf::Vector2f cellPos = gridCellToPos(cellID);

        rect.setPosition(cellPos);
        window.draw(rect);

        /*text.setString(std::to_string(phi[cellID]));
        text.setPosition(cellPos);
        window.draw(text);*/
    }

    void DrawSolidCells() {
        for (int i = 0; i < nX; ++i) {
            for (int j = 0; j < nY; ++j) {
                int idx = i * nY + j;
                if (grid[idx] == 1) {
                    DrawSolidCellAt(idx);
                }
            }
        }
    }

    void FastSweep(float solidCells) {
        for (int i = 0; i < gridSize; ++i) {
            if (grid[i] == SOLID) {
                phi[i] = 0.0f;            // inside obstacle
            }
            else {
                phi[i] = 1e6f;            // large distance for air
            }
        }

        const int NSweeps = 4;
        // sweep directions { start, end, step }
        const int dirX[NSweeps][3] = { {0, nX - 1,  1}, {nX - 1, 0, -1}, {nX - 1, 0,  -1}, {0, nX - 1,   1} };
        const int dirY[NSweeps][3] = { {0, nY - 1, 1}, {0, nY - 1, 1}, {nY - 1, 0, -1}, {nY - 1, 0, -1} };

        double aa[2], eps = 1e-6;
        double d_new, a, b;
        int s, ix, iy, gridPos;
        const double h = 1.0, f = 1.0;

        for (s = 0; s < NSweeps; s++) {
            for (ix = dirX[s][0]; dirX[s][2] * ix <= dirX[s][1]; ix += dirX[s][2]) {
                for (iy = dirY[s][0]; dirY[s][2] * iy <= dirY[s][1]; iy += dirY[s][2]) {
                    gridPos = ix * nY + iy;
                    int left   = gridPos - nY;
                    int right  = gridPos + nY;
                    int top    = gridPos - 1;
                    int bottom = gridPos + 1;

                    if (!grid[gridPos]) {
                        if (iy == 0 || iy == (nY - 1)) {
                            if (iy == 0) {
                                aa[1] = phi[gridPos] < phi[bottom] ? phi[gridPos] : phi[bottom];
                            }
                                if (iy == (nY - 1)) {
                                aa[1] = phi[top] < phi[gridPos] ? phi[top] : phi[gridPos];
                            }
                        }
                        else {
                            aa[1] = phi[top] < phi[bottom] ? phi[top] : phi[bottom];
                        }
                    
                        if (ix == 0 || ix == (nX - 1)) {
                            if (ix == 0) {
                                aa[0] = phi[gridPos] < phi[right] ? phi[gridPos] : phi[right];
                            }
                            if (ix == (nX - 1)) {
                                aa[0] = phi[left] < phi[gridPos] ? phi[left] : phi[gridPos];
                            }
                        }
                        else {
                            aa[0] = phi[left] < phi[right] ? phi[left] : phi[right];
                        }
                    
                        a = aa[0]; b = aa[1];
                        d_new = (fabs(a - b) < f * h ? (a + b + sqrt(2.0 * f * f * h * h - (a - b) * (a - b))) * 0.5 : std::fminf(a, b) + f * h);
                    
                        phi[gridPos] = phi[gridPos] < d_new ? phi[gridPos] : d_new;
                    }
                }
            }
        }
    }

    void ComputeGradients() {
        for (int i = 0; i < nX; ++i) {
            for (int j = 0; j < nY; ++j) {
                int idx = i * nY + j;
                float x = i * cellSpacing;
                float y = j * cellSpacing;
                dphi[idx] = QueryLevelSet(sf::Vector2f{x, y});
            }
        }
    }

    sf::Vector2f QueryLevelSet(sf::Vector2f pos) {
        float x = pos.x;
        float y = pos.y;

        float d = samplePhi(pos);

        float eps = 1e-1;

        float ddx = (samplePhi(pos + sf::Vector2f{pos.x + eps, pos.y}) - samplePhi(pos + sf::Vector2f{pos.x - eps, pos.y})) / (2 * eps);

        float ddy = (samplePhi(pos + sf::Vector2f{pos.x, pos.y + eps}) - samplePhi(pos + sf::Vector2f{pos.x, pos.y - eps})) / (2 * eps);

        float div = sqrt(ddx * ddx + ddy * ddy);
        sf::Vector2f grad = sf::Vector2{ddx / div, ddy / div};

        return pos - d * grad;
    }

    float samplePhi(sf::Vector2f pos) {
        int x = pos.x;
        int y = pos.y;

        int i = std::clamp(std::floor(x / cellSpacing), 0.f, static_cast<float>(nX - 1));
        int j = std::clamp(std::floor(y / cellSpacing), 0.f, static_cast<float>(nY - 1));

        //return phi[i * nY + j];

        float fx = (x / cellSpacing - i);
        float fy = (y / cellSpacing - j);

        int idx = i * nY + j;
        
        float topLeft = phi[idx];
        float topRight = phi[idx + nY];
        float bottomLeft = phi[idx + 1];
        float bottomRight = phi[idx + nY + 1];

        float v0 = (1.f - fx) * topLeft + fx * topRight;
        float v1 = (1.f - fx) * bottomLeft + fx * bottomRight;

        return (1.f - fy) * v0 + fy * v1;
    }

    sf::Vector2f sampleGradient(sf::Vector2f pos) {
        float eps = cellSpacing * 0.5f; // small step in world units
        float ddx = (samplePhi(pos + sf::Vector2f{eps, 0}) -
                     samplePhi(pos - sf::Vector2f{eps, 0})) / (2 * eps);
        float ddy = (samplePhi(pos + sf::Vector2f{0, eps}) -
                     samplePhi(pos - sf::Vector2f{0, eps})) / (2 * eps);

        sf::Vector2f g{ddx, ddy};
        float len = std::sqrt(g.x * g.x + g.y * g.y);
        return (len > 1e-6f) ? g / len : sf::Vector2f{0,0};
    }

    sf::Vector2f ClosestSurfacePoint(sf::Vector2f pos) {
        float d = samplePhi(pos);
        sf::Vector2f grad = sampleGradient(pos);
        return pos - d * grad;
    }

    void DrawClosestSurfacePoint() {
        sf::Vector2f mousePos = static_cast<sf::Vector2f>(sf::Mouse::getPosition(window));
        //sf::Vector2f point = QueryLevelSet(mousePos);
        sf::Vector2f surfacePoint = ClosestSurfacePoint(mousePos);

        /*text.setPosition(mousePos.x, mousePos.y);
        std::string grad = "{" + std::to_string(point.x) + ", " + std::to_string(point.y) + "}";
        text.setString(grad);
        window.draw(text);*/

        circle.setPosition(surfacePoint);

        window.draw(circle);
    }

    void DrawSDF() {
        // Find min and max phi for normalization
        auto [minIt, maxIt] = std::minmax_element(phi.begin(), phi.end());
        float minPhi = *minIt;
        float maxPhi = *maxIt;

        // Avoid division by zero
        float range = (maxPhi - minPhi > 1e-6f) ? (maxPhi - minPhi) : 1.f;

        for (int i = 0; i < WIDTH / stepSize; ++i) {
            for (int j = 0; j < HEIGHT / stepSize; ++j) {
                float x = i * stepSize;
                float y = j * stepSize;
                sf::Vector2f pos = {x, y};

                float dist = samplePhi(pos);

                // Normalize dist → [0,1]
                float intensity = (dist - minPhi) / range;
                intensity = std::clamp(intensity, 0.f, 1.f) * 255.f;

                sf::Color color(intensity, intensity, intensity);

                tinyRect.setFillColor(color);
                tinyRect.setPosition(pos);
                window.draw(tinyRect);
            }
        }
    }

    void update() {
        DrawSDF();
        //DrawSolidCells();
        DrawClosestSurfacePoint();
    }

    void printPhiAndGrid() {
        std::cout << "phi:" << std::endl;
        std::cout << "[";
        for (int j = 0; j < nY; ++j) {
            for (int i = 0; i < nX; ++i) {
                int idx = i * nY + j;
                std::cout << phi[idx] << ", ";
            }
            if (j < nY - 1) {
                std::cout << std::endl;
            }
            else {
                std::cout << "]";
            }
        }

        std::cout << std::endl << std::endl;

        std::cout << "grid:" << std::endl;
        std::cout << "[";
        for (int j = 0; j < nY; ++j) {
            for (int i = 0; i < nX; ++i) {
                int idx = i * nY + j;
                std::cout << grid[idx] << ", ";
            }
            if (j < nY - 1) {
                std::cout << std::endl;
            }
            else {
                std::cout << "]";
            }
        }
    }
};