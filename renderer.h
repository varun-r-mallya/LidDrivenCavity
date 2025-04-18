#ifndef RENDERER_H
#define RENDERER_H

#include <SFML/Graphics.hpp>
#include <vector>
#include <memory>

class Renderer {
public:
    Renderer(int width, int height,
             std::vector<std::vector<double> > &u,
             std::vector<std::vector<double> > &v,
             std::vector<double> &x,
             std::vector<double> &y);

    void initialize();

    bool isWindowOpen() const;

    void handleEvents() const;

    void render() const;

    void updateData(const std::vector<std::vector<double> > &u,
                    const std::vector<std::vector<double> > &v);

private:
    void drawVelocityField() const;

    void drawColorMap() const;

    void calculateVorticity();

    void findMinMaxValues();

    std::unique_ptr<sf::RenderWindow> window;
    int windowWidth;
    int windowHeight;

    std::vector<std::vector<double> > &u;
    std::vector<std::vector<double> > &v;
    std::vector<double> &x;
    std::vector<double> &y;

    std::vector<std::vector<double> > vorticity;
    double minVorticity;
    double maxVorticity;
};

#endif // RENDERER_H
