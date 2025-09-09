#pragma once

#include <SFML/Graphics.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <iostream>

#include <chrono>

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

    sf::VertexArray SDFva{sf::PrimitiveType::Quads};
    sf::Texture texture;
    sf::RenderStates states;
    sf::Vector2f texture_size;
    float stepSize = 2.f;
    int vaGridSize;

    sf::RenderWindow &window;

    std::mt19937 gen;
    std::bernoulli_distribution distrib;

    sf::Text text;

    sf::CircleShape circle;

    float WIDTH;
    float HEIGHT;

    int numStepIters = 2;
    int numInsideStepIters = 2;

    bool drawSolidCells = true;

    bool drawing = false;
    bool erasing = false;

    bool wave = true;

    int drawerRadius = 1;

    std::array<std::array<float, 3>, 2> SDFcolors = {{{255.f, 150.f, 0.f}, {0.f, 80.f, 255.f}}};

    std::array<std::array<int, 3>, 50> SDFgradient;

    LevelSet(sf::Text text_, sf::RenderWindow &window_, int nX_, float WIDTH_, float HEIGHT_): WIDTH(WIDTH_), HEIGHT(HEIGHT_), text(text_), window(window_), nX(nX_), gen(std::random_device{}()), distrib(0.5) {
        cellSpacing = WIDTH / nX;
        halfSpacing = cellSpacing * 0.5;

        rect.setFillColor(sf::Color::White);
        rect.setSize(sf::Vector2f{cellSpacing, cellSpacing});
        rect.setOutlineThickness(1.f);
        rect.setOutlineColor(sf::Color::White);

        vaGridSize = (WIDTH / stepSize) * (HEIGHT / stepSize);
        SDFva.resize(4 * vaGridSize);
        texture.loadFromFile("white_square.png");
        //texture.generateMipmap();
        texture_size = static_cast<sf::Vector2f>(texture.getSize());
        for (uint32_t i = 0; i < WIDTH / stepSize; ++i) {
            for (uint32_t j = 0; j < HEIGHT / stepSize; ++j) {
                int index = i * (HEIGHT / stepSize) + j;
                int vaIdx = 4 * index;

                float halfStepSize = stepSize * 0.5f;
                float px = i * stepSize + halfStepSize;
                float py = j * stepSize + halfStepSize;
                SDFva[vaIdx].position = {px - halfStepSize, py - halfStepSize};
                SDFva[vaIdx + 1].position = {px + halfStepSize, py - halfStepSize};
                SDFva[vaIdx + 2].position = {px + halfStepSize, py + halfStepSize};
                SDFva[vaIdx + 3].position = {px - halfStepSize, py + halfStepSize};

                SDFva[vaIdx].texCoords = {0.f, 0.f};
                SDFva[vaIdx + 1].texCoords = {texture_size.x, 0.f};
                SDFva[vaIdx + 2].texCoords = {texture_size.x, texture_size.y};
                SDFva[vaIdx + 3].texCoords = {0.f, texture_size.y};
            }
        }
        states.texture = &texture;

        float num_colors = SDFcolors.size() - 1; // number of colors - 1
        float num_steps = 1.f * SDFgradient.size() / num_colors; //num_steps = 50 * key_range
        int index = 0;
        for (int i = 0; i < num_colors; ++i) {  
            for (int x = 0; x < num_steps; ++x) {
                float t = 1.f * x / num_steps;
                int r = (int)(SDFcolors[i][0] * (1 - t) + SDFcolors[i + 1][0] * t);
                int g = (int)(SDFcolors[i][1] * (1 - t) + SDFcolors[i + 1][1] * t);
                int b = (int)(SDFcolors[i][2] * (1 - t) + SDFcolors[i + 1][2] * t);
                SDFgradient[index] = std::array<int, 3>{r, g, b};

                index++;
            }
        }

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

        for (int i = 0; i < gridSize; ++i) {
            grid[i] = distrib(gen);
        }

        int donutOuterRad = nX / 3.333333;
        int donutInnerRad = nX / 6.25;

        CreateDonut(donutOuterRad, donutInnerRad);

        int squareWidth = nX / 10;
        CreateCenterSquare(squareWidth);

        FastSweep();
    }

    void CreateDonut(int outerRadius, int innerRadius) {
        int centerX = nX / 2;
        int centerY = nY / 2;
        int padding = 3;
        for (int ix = centerX - (outerRadius + padding); ix < centerX + (outerRadius + padding); ++ix) {
            for (int iy = centerY - (outerRadius + padding); iy < centerY + (outerRadius + padding); ++iy) {
                int idx = ix * nY + iy;
            
                float dx = ix - nX / 2;
                float dy = iy - nY / 2;
                float d2 = dx*dx + dy*dy;
            
                grid[idx] = (d2 <= outerRadius*outerRadius && d2 >= innerRadius*innerRadius) ? SOLID : AIR;
            }
        }
    }

    void CreateCenterSquare(int width) {
        int halfWidth = width / 2;
        int centerX = nX / 2;
        int centerY = nY / 2;

        for (int ix = centerX - halfWidth; ix <= centerX + halfWidth; ++ix) {
            for (int iy = centerY - halfWidth; iy <= centerY + halfWidth; ++iy) {
                int idx = ix * nY + iy;
                grid[idx] = SOLID;
            }
        }
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

    void SetUpPhi() {
        for (int i = 0; i < gridSize; ++i) {
            if (grid[i] == SOLID) {
                phi[i] = 1e10f;
            }
            else {
                phi[i] = 0.f;
            }
        }
    }

    void FastSweep() {
        auto start = std::chrono::high_resolution_clock::now();

        SetUpPhi();

        const int NSweeps = 4;
        // sweep directions { start, end, step }
        const int dirX[NSweeps][3] = { {0, nX - 1,  1}, {nX - 1, 0, -1}, {nX - 1, 0,  -1}, {0, nX - 1,   1} };
        const int dirY[NSweeps][3] = { {0, nY - 1, 1}, {0, nY - 1, 1}, {nY - 1, 0, -1}, {nY - 1, 0, -1} };

        double aa[2], eps = 1e-6;
        double d_new, a, b;
        int s, ix, iy, gridPos;
        const double h = cellSpacing, f = 1.0;

        for (s = 0; s < NSweeps; s++) {
            for (ix = dirX[s][0]; dirX[s][2] * ix <= dirX[s][1]; ix += dirX[s][2]) {
                for (iy = dirY[s][0]; dirY[s][2] * iy <= dirY[s][1]; iy += dirY[s][2]) {
                    gridPos = ix * nY + iy;
                    int left   = gridPos - nY;
                    int right  = gridPos + nY;
                    int top    = gridPos - 1;
                    int bottom = gridPos + 1;

                    if (grid[gridPos]) {
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

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        //std::cout << duration.count() << " milliseconds\n";
    }

    void ComputeGradients() {
        for (int i = 0; i < nX; ++i) {
            for (int j = 0; j < nY; ++j) {
                int idx = i * nY + j;
                float x = i * cellSpacing;
                float y = j * cellSpacing;
                dphi[idx] = sampleGradient(sf::Vector2f{x, y});
            }
        }
    }

    float samplePhi(sf::Vector2f pos) {
        int x = pos.x;
        int y = pos.y;

        float gx = x / cellSpacing - 0.5f;
        float gy = y / cellSpacing - 0.5f;

        int i0 = std::clamp(int(std::floor(gx)), 0, nX - 1);
        int j0 = std::clamp(int(std::floor(gy)), 0, nY - 1);
        int i1 = std::min(i0 + 1, nX - 1);
        int j1 = std::min(j0 + 1, nY - 1);
        
        float fx = gx - i0;
        float fy = gy - j0;

        int topLeft     = i0 * nY + j0;
        int topRight    = i1 * nY + j0;
        int bottomLeft  = i0 * nY + j1;
        int bottomRight = i1 * nY + j1;

        float topLeftVal     = phi[topLeft];
        float topRightVal    = phi[topRight];
        float bottomLeftVal  = phi[bottomLeft];
        float bottomRightVal = phi[bottomRight];

        float v0 = (1 - fx) * topLeftVal    + fx * topRightVal;
        float v1 = (1 - fx) * bottomLeftVal + fx * bottomRightVal;
        return (1 - fy) * v0 + fy * v1;
    }

    sf::Vector2f sampleGradient(sf::Vector2f pos) {
        float eps = cellSpacing * 0.25f;
        
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

    sf::Vector2f ClosestSurfacePointIterative(sf::Vector2f pos) {
        auto p = pos;
        auto phi_p = samplePhi(p);
        if (phi_p == 0) {
            return pos;
        }
        auto d = sampleGradient(p);

        float eps = 1e-1;

        for (int i = 0; i < numStepIters; ++i) {
            float alpha = 1.f;

            for (int i = 0; i < numInsideStepIters; ++i) {
                auto q = p - alpha * phi_p * d;
                auto phi_q = samplePhi(q);
                if (abs(phi_q) < abs(phi_p)) {
                    p = q;
                    phi_p = phi_q;
                    d = sampleGradient(q);
                    if (abs(phi_p) < eps) {
                        return p;
                    }
                }
                else {
                    alpha = 0.7 * alpha;
                }
            }
        }

        return p;
    }

    void DrawClosestSurfacePoint() {
        sf::Vector2f mousePos = static_cast<sf::Vector2f>(sf::Mouse::getPosition(window));
        //sf::Vector2f surfacePoint = ClosestSurfacePoint(mousePos);
        sf::Vector2f surfacePoint = ClosestSurfacePointIterative(mousePos);

        circle.setPosition(surfacePoint);
        window.draw(circle);
    }

    void DrawSDF() {
        auto [minIt, maxIt] = std::minmax_element(phi.begin(), phi.end());
        float minPhi = *minIt;
        float maxPhi = *maxIt;

        float range = (maxPhi - minPhi > 1e-6f) ? (maxPhi - minPhi) : 1.f;
        for (uint32_t i = 0; i < WIDTH / stepSize; ++i) {
            for (uint32_t j = 0; j < HEIGHT / stepSize; ++j) {
                int index = i * (HEIGHT / stepSize) + j;
                int vaIdx = 4 * index;

                float px = i * stepSize + stepSize * 0.5f;
                float py = j * stepSize + stepSize * 0.5f;

                int gx = px / cellSpacing;
                int gy = py / cellSpacing;
                int gIdx = gx * nY + gy;

                sf::Vector2f pos = {px, py};

                float dist = samplePhi(pos);
                if (grid[gIdx] == AIR) {
                    sf::Color black(0, 0, 0);
                    SDFva[vaIdx].color = black;
                    SDFva[vaIdx + 1].color = black;
                    SDFva[vaIdx + 2].color = black;
                    SDFva[vaIdx + 3].color = black;
                    continue;
                }

                if (wave) {
                    dist = maxPhi * sin(dist / 2.f);
                }

                float intensity = (dist - minPhi) / range;
                intensity = static_cast<int>(std::clamp(intensity, 0.f, 1.f) * (SDFgradient.size() - 1));

                auto gradientColor = SDFgradient[intensity];

                sf::Color color(gradientColor[0], gradientColor[1], gradientColor[2]);

                SDFva[vaIdx].color = color;
                SDFva[vaIdx + 1].color = color;
                SDFva[vaIdx + 2].color = color;
                SDFva[vaIdx + 3].color = color;
            }
        }
        window.draw(SDFva, states);
    }

    void update() {
        if (drawing) {
            drawSolids();
        }
        if (erasing) {
            eraseSolids();
        }
        if (!drawSolidCells) {
            DrawSDF();
        }
        if (drawSolidCells) {
            DrawSolidCells();
        }
        DrawClosestSurfacePoint();

        displayNumStepIters();
        displayInsideNumStepIters();
    }

    void changeDrawSolidCells() {
        drawSolidCells = !drawSolidCells;
    }

    void increaseNumStepIters() {
        numStepIters++;
    }

    void decreaseNumStepIters() {
        if (numStepIters > 1) {
            numStepIters--;
        }
    }

    void increaseInsideNumStepIters() {
        numInsideStepIters++;
    }

    void decreaseInsideNumStepIters() {
        if (numInsideStepIters > 1) {
            numInsideStepIters--;
        }
    }

    void displayNumStepIters() {
        text.setPosition(WIDTH - 50, 20);
        text.setString(std::to_string(numStepIters));
        window.draw(text);
    }

    void displayInsideNumStepIters() {
        text.setPosition(WIDTH - 20, 20);
        text.setString(std::to_string(numInsideStepIters));
        window.draw(text);
    }

    void drawSolids() {
        sf::Vector2i pos = sf::Mouse::getPosition(window);

        int gx = pos.x / cellSpacing;
        int gy = pos.y / cellSpacing;

        for (int i = -drawerRadius; i <= drawerRadius; ++i) {
            for (int j = -drawerRadius; j <= drawerRadius; ++j) {
                int x = gx + i;
                int y = gy + j;

                if (x >= 0 && y >= 0 && x <= nX - 1 && y <= nY - 1) {

                    int idx = x * nY + y;

                    grid[idx] = SOLID;
                }
            }
        }

        FastSweep();
    }

    void eraseSolids() {
        sf::Vector2i pos = sf::Mouse::getPosition(window);

        int gx = pos.x / cellSpacing;
        int gy = pos.y / cellSpacing;

        for (int i = -drawerRadius; i <= drawerRadius; ++i) {
            for (int j = -drawerRadius; j <= drawerRadius; ++j) {
                int x = gx + i;
                int y = gy + j;

                if (x >= 0 && y >= 0 && x <= nX - 1 && y <= nY - 1) {

                    int idx = x * nY + y;

                    grid[idx] = AIR;
                }
            }
        }

        FastSweep();
    }

    void SetDrawer(bool set) {
        drawing = set;
    }

    void SetEraser(bool set) {
        erasing = set;
    }

    void AddToDrawerRadius(int add) {
        if (drawerRadius + add < 0) {
            drawerRadius = 0;
        }
        else {
            drawerRadius += add;
        }
    }

    void SwitchSDFRendering() {
        wave = !wave;
    }
};