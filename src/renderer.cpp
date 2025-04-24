#include "renderer.h"
#include <cmath>
#include <algorithm>

Renderer::Renderer(int width, int height,
                   std::vector<std::vector<double> > &u,
                   std::vector<std::vector<double> > &v,
                   std::vector<double> &x,
                   std::vector<double> &y)
    : windowWidth(width), windowHeight(height),
      u(u), v(v), x(x), y(y),
      maxVelocity(0.0), gridSize(u.size()) {
}

void Renderer::initialize() {
    window = std::make_unique<sf::RenderWindow>(
        sf::VideoMode(windowWidth, windowHeight),
        "Lid Driven Cavity",
        sf::Style::Titlebar | sf::Style::Close
    );
    window->setVerticalSyncEnabled(true);
    // Calculate max velocity for normalization
    for (const auto &row: u) {
        for (double val: row) {
            maxVelocity = std::max(maxVelocity, std::abs(val));
        }
    }
    for (const auto &row: v) {
        for (double val: row) {
            maxVelocity = std::max(maxVelocity, std::abs(val));
        }
    }

    // Initialize render texture for anti-aliased rendering
    renderTexture.create(windowWidth, windowHeight);
}

bool Renderer::isWindowOpen() const {
    return window && window->isOpen();
}

void Renderer::handleEvents() {
    if (!window) return;

    sf::Event event;
    while (window->pollEvent(event)) {
        if (event.type == sf::Event::Closed) {
            window->close();
        }
    }
}

void Renderer::render() {
    if (!window) return;

    renderTexture.clear(sf::Color::White);
    drawColorMap();
    drawStreamlines();
    renderTexture.display();

    window->clear();
    sf::Sprite sprite(renderTexture.getTexture());
    // Single line to rotate contents upside down:
    sprite.setRotation(180); // Rotate 180 degrees
    sprite.setOrigin(sprite.getLocalBounds().width / 2, sprite.getLocalBounds().height / 2); // Rotate around center
    sprite.setPosition(windowWidth / 2, windowHeight / 2); // Center in window
    sprite.setScale(-1, 1);
    window->draw(sprite);
    window->display();
}

void Renderer::updateData(const std::vector<std::vector<double> > &u,
                          const std::vector<std::vector<double> > &v) {
    this->u = u;
    this->v = v;

    // Recalculate max velocity
    maxVelocity = 0.0;
    for (const auto &row: u) {
        for (double val: row) {
            maxVelocity = std::max(maxVelocity, std::abs(val));
        }
    }
    for (const auto &row: v) {
        for (double val: row) {
            maxVelocity = std::max(maxVelocity, std::abs(val));
        }
    }
}

void Renderer::drawColorMap() {
    const float cellWidth = static_cast<float>(windowWidth) / gridSize;
    const float cellHeight = static_cast<float>(windowHeight) / gridSize;

    sf::VertexArray vertices(sf::Quads, gridSize * gridSize * 4);

    for (int i = 0; i < gridSize; ++i) {
        for (int j = 0; j < gridSize; ++j) {
            double vel = std::sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]);
            sf::Color color = getVelocityColor(vel);

            int index = (i * gridSize + j) * 4;

            vertices[index].position = sf::Vector2f(j * cellWidth, i * cellHeight);
            vertices[index + 1].position = sf::Vector2f((j + 1) * cellWidth, i * cellHeight);
            vertices[index + 2].position = sf::Vector2f((j + 1) * cellWidth, (i + 1) * cellHeight);
            vertices[index + 3].position = sf::Vector2f(j * cellWidth, (i + 1) * cellHeight);

            vertices[index].color = color;
            vertices[index + 1].color = color;
            vertices[index + 2].color = color;
            vertices[index + 3].color = color;
        }
    }

    renderTexture.draw(vertices);
}

void Renderer::drawStreamlines() {
    const int numSeedsX = 15;
    const int numSeedsY = 15;
    const float minVelocity = maxVelocity * 0.05f;

    streamlines.clear();

    // Generate seed points
    for (int i = 0; i < numSeedsX; ++i) {
        for (int j = 0; j < numSeedsY; ++j) {
            float x = static_cast<float>(i + 0.5f) / numSeedsX * windowWidth;
            float y = static_cast<float>(j + 0.5f) / numSeedsY * windowHeight;

            // Check if velocity at this point is significant
            float u_val, v_val;
            bilinearInterpolate(x, y, u_val, v_val);
            if (std::sqrt(u_val * u_val + v_val * v_val) > minVelocity) {
                generateStreamline(x, y);
            }
        }
    }

    // Draw all streamlines
    for (const auto &streamline: streamlines) {
        renderTexture.draw(streamline);
    }
}

void Renderer::generateStreamline(float startX, float startY) {
    const float stepSize = 0.5f;
    const int maxSteps = 200;
    const float minVelocity = maxVelocity * 0.01f;

    sf::VertexArray forwardLine(sf::LineStrip);
    sf::VertexArray backwardLine(sf::LineStrip);

    // Trace forward
    float x = startX;
    float y = startY;
    for (int i = 0; i < maxSteps; ++i) {
        float u_val, v_val;
        bilinearInterpolate(x, y, u_val, v_val);

        float vel = std::sqrt(u_val * u_val + v_val * v_val);
        if (vel < minVelocity) break;

        forwardLine.append(sf::Vertex(sf::Vector2f(x, y),
                                      getVelocityColor(vel * 1.5f)));

        // RK4 integration
        float k1x = u_val;
        float k1y = v_val;

        float k2x, k2y;
        bilinearInterpolate(x + k1x * stepSize / 2, y + k1y * stepSize / 2, k2x, k2y);

        float k3x, k3y;
        bilinearInterpolate(x + k2x * stepSize / 2, y + k2y * stepSize / 2, k3x, k3y);

        float k4x, k4y;
        bilinearInterpolate(x + k3x * stepSize, y + k3y * stepSize, k4x, k4y);

        x += (k1x + 2 * k2x + 2 * k3x + k4x) * stepSize / 6;
        y += (k1y + 2 * k2y + 2 * k3y + k4y) * stepSize / 6;

        if (x < 0 || x >= windowWidth || y < 0 || y >= windowHeight) break;
    }

    // Trace backward
    x = startX;
    y = startY;
    for (int i = 0; i < maxSteps; ++i) {
        float u_val, v_val;
        bilinearInterpolate(x, y, u_val, v_val);

        float vel = std::sqrt(u_val * u_val + v_val * v_val);
        if (vel < minVelocity) break;

        backwardLine.append(sf::Vertex(sf::Vector2f(x, y),
                                       getVelocityColor(vel * 1.5f)));

        // RK4 integration (backward)
        float k1x = -u_val;
        float k1y = -v_val;

        float k2x, k2y;
        bilinearInterpolate(x + k1x * stepSize / 2, y + k1y * stepSize / 2, k2x, k2y);

        float k3x, k3y;
        bilinearInterpolate(x + k2x * stepSize / 2, y + k2y * stepSize / 2, k3x, k3y);

        float k4x, k4y;
        bilinearInterpolate(x + k3x * stepSize, y + k3y * stepSize, k4x, k4y);

        x += (k1x + 2 * k2x + 2 * k3x + k4x) * stepSize / 6;
        y += (k1y + 2 * k2y + 2 * k3y + k4y) * stepSize / 6;

        if (x < 0 || x >= windowWidth || y < 0 || y >= windowHeight) break;
    }

    // Combine forward and backward traces
    sf::VertexArray fullLine(sf::LineStrip);
    for (int i = backwardLine.getVertexCount() - 1; i >= 0; --i) {
        fullLine.append(backwardLine[i]);
    }
    for (unsigned int i = 0; i < forwardLine.getVertexCount(); ++i) {
        fullLine.append(forwardLine[i]);
    }

    if (fullLine.getVertexCount() > 1) {
        streamlines.push_back(fullLine);
    }
}

sf::Color Renderer::getVelocityColor(double magnitude) const {
    if (maxVelocity < 1e-6) return sf::Color::Black;

    double normalized = magnitude / maxVelocity;
    normalized = std::max(0.0, std::min(1.0, normalized));

    // Blue (low) -> Cyan -> Green -> Yellow -> Red (high)
    if (normalized < 0.25) {
        return sf::Color(
            0,
            static_cast<sf::Uint8>(normalized * 4 * 255),
            255
        );
    } else if (normalized < 0.5) {
        return sf::Color(
            0,
            255,
            static_cast<sf::Uint8>((1.0 - (normalized - 0.25) * 4) * 255)
        );
    } else if (normalized < 0.75) {
        return sf::Color(
            static_cast<sf::Uint8>((normalized - 0.5) * 4 * 255),
            255,
            0
        );
    } else {
        return sf::Color(
            255,
            static_cast<sf::Uint8>((1.0 - (normalized - 0.75) * 4) * 255),
            0
        );
    }
}

void Renderer::bilinearInterpolate(float x, float y, float &u_out, float &v_out) const {
    float gridX = x / windowWidth * (gridSize - 1);
    float gridY = y / windowHeight * (gridSize - 1);

    int x0 = static_cast<int>(gridX);
    int y0 = static_cast<int>(gridY);
    int x1 = std::min(x0 + 1, gridSize - 1);
    int y1 = std::min(y0 + 1, gridSize - 1);

    float sx = gridX - x0;
    float sy = gridY - y0;

    // Bilinear interpolation for u
    float u00 = u[y0][x0];
    float u01 = u[y0][x1];
    float u10 = u[y1][x0];
    float u11 = u[y1][x1];

    float u0 = u00 + sx * (u01 - u00);
    float u1 = u10 + sx * (u11 - u10);
    u_out = u0 + sy * (u1 - u0);

    // Bilinear interpolation for v
    float v00 = v[y0][x0];
    float v01 = v[y0][x1];
    float v10 = v[y1][x0];
    float v11 = v[y1][x1];

    float v0 = v00 + sx * (v01 - v00);
    float v1 = v10 + sx * (v11 - v10);
    v_out = v0 + sy * (v1 - v0);
}
