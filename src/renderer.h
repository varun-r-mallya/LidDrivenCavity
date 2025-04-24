#ifndef RENDERER_H
#define RENDERER_H

#include <SFML/Graphics.hpp>
#include <vector>
#include <memory>
#include <list>

class Renderer {
public:
    Renderer(int width, int height,
             std::vector<std::vector<double> > &u,
             std::vector<std::vector<double> > &v,
             std::vector<double> &x,
             std::vector<double> &y);

    void initialize();

    [[nodiscard]] bool isWindowOpen() const;

    void handleEvents();

    void render();

    void updateData(const std::vector<std::vector<double> > &u,
                    const std::vector<std::vector<double> > &v);

private:
    void drawColorMap();

    void drawStreamlines();

    void generateStreamline(float startX, float startY);

    sf::Color getVelocityColor(double magnitude) const;

    void bilinearInterpolate(float x, float y, float &u, float &v) const;

    std::unique_ptr<sf::RenderWindow> window;
    mutable sf::RenderTexture renderTexture; // Made mutable
    int windowWidth;
    int windowHeight;

    std::vector<std::vector<double> > &u;
    std::vector<std::vector<double> > &v;
    std::vector<double> &x;
    std::vector<double> &y;

    double maxVelocity;
    std::vector<sf::VertexArray> streamlines;
    int gridSize;
};

#endif // RENDERER_H
