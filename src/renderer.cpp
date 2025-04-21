#include "renderer.h"
#include <cmath>
#include <algorithm>

void Renderer::initialize() {
    window = std::make_unique<sf::RenderWindow>(
        sf::VideoMode(windowWidth, windowHeight),
        "Lid-Driven Cavity Flow"
    );
    calculateVorticity();
    findMinMaxValues();
}

Renderer::Renderer(const int width, const int height,
                   std::vector<std::vector<double> > &u,
                   std::vector<std::vector<double> > &v,
                   std::vector<double> &x,
                   std::vector<double> &y)
    : windowWidth(width), windowHeight(height),
      u(u), v(v), x(x), y(y),
      vorticity(u.size(), std::vector(u[0].size(), 0.0)), minVorticity(0), maxVorticity(0) {
}

bool Renderer::isWindowOpen() const {
    return window && window->isOpen();
}

void Renderer::handleEvents() const {
    if (!window) return;

    sf::Event event;
    while (window->pollEvent(event)) {
        if (event.type == sf::Event::Closed) {
            window->close();
        }
    }
}

void Renderer::render() const {
    if (!window) return;

    window->clear(sf::Color::White);
    drawColorMap();
    window->display();
}

void Renderer::updateData(const std::vector<std::vector<double> > &u,
                          const std::vector<std::vector<double> > &v) {
    this->u = u;
    this->v = v;
    calculateVorticity();
    findMinMaxValues();
}

void Renderer::drawColorMap() const {
    const int gridSize = u.size();
    const float cellSize = std::min(windowWidth, windowHeight) / static_cast<float>(gridSize);
    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            double normalized = (vorticity[i][j] - minVorticity) / (maxVorticity - minVorticity);
            normalized = std::max(0.0, std::min(1.0, normalized));

            sf::Color color;

            if (normalized < 0.15) {
                // Black to Deep Blue
                color.r = 0;
                color.g = 0;
                color.b = static_cast<sf::Uint8>(normalized / 0.15 * 255 * 0.5);  // Max 127
            } else if (normalized < 0.3) {
                // Blue to Purple
                color.r = static_cast<sf::Uint8>((normalized - 0.15) / 0.15 * 200);
                color.g = 0;
                color.b = 127 + static_cast<sf::Uint8>((normalized - 0.15) / 0.15 * 128);
            } else if (normalized < 0.5) {
                // Purple to Magenta to Red
                color.r = 200 + static_cast<sf::Uint8>((normalized - 0.3) / 0.2 * 55);
                color.g = 0;
                color.b = static_cast<sf::Uint8>((0.5 - normalized) / 0.2 * 255);
            } else if (normalized < 0.7) {
                // Red to Orange-Yellow
                color.r = 255;
                color.g = static_cast<sf::Uint8>((normalized - 0.5) / 0.2 * 255);
                color.b = 0;
            } else if (normalized < 0.9) {
                // Orange-Yellow to Yellow-White
                color.r = 255;
                color.g = 255;
                color.b = static_cast<sf::Uint8>((normalized - 0.7) / 0.2 * 200);
            } else {
                // Bright white flash for high vorticity
                sf::Uint8 intensity = static_cast<sf::Uint8>(200 + (normalized - 0.9) / 0.1 * 55);
                color.r = intensity;
                color.g = intensity;
                color.b = intensity;
            }

            sf::RectangleShape rect(sf::Vector2f(cellSize, cellSize));
            rect.setPosition(i * cellSize, j * cellSize);
            rect.setFillColor(color);
            window->draw(rect);
        }
    }
}

void Renderer::calculateVorticity() {
    const int gridSize = u.size();
#pragma omp parallel for
    for (int i = 1; i < gridSize - 1; ++i) {
        for (int j = 1; j < gridSize - 1; ++j) {
            vorticity[i][j] = (v[i][j + 1] - v[i][j - 1]) / (x[j + 1] - x[j - 1]) -
                              (u[i + 1][j] - u[i - 1][j]) / (y[i + 1] - y[i - 1]);
        }
    }
}

void Renderer::findMinMaxValues() {
    minVorticity = vorticity[1][1];
    maxVorticity = vorticity[1][1];

    const int gridSize = u.size();
#pragma omp parallel for reduction(min:minVorticity) reduction(max:maxVorticity)
    for (int i = 1; i < gridSize - 1; ++i) {
        for (int j = 1; j < gridSize - 1; ++j) {
            if (vorticity[i][j] < minVorticity) minVorticity = vorticity[i][j];
            if (vorticity[i][j] > maxVorticity) maxVorticity = vorticity[i][j];
        }
    }
}
