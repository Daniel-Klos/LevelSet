#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <array>
#include <sstream>
#include <iomanip>

#include "levelset.hpp"

int main()
{
    // Adjust the sim to a good size for you. Note that if you change these numbers, you will have to change the settings
    int WIDTH = 1500; // 2500 
    int HEIGHT = 1200; // 1300

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "level set");

    sf::Font font;
    font.loadFromFile("C:\\Users\\dklos\\vogue\\Vogue.ttf");

    sf::Text text;
    text.setFont(font);
    text.setPosition(10, 10);
    text.setFillColor(sf::Color::Red);

    sf::Clock deltaClock;

    window.setFramerateLimit(120); // 120

    int frame = 0;
    int fps = 0;

    int nX = 60;
    LevelSet SDF = LevelSet(text, window, nX, WIDTH, HEIGHT);

    while (window.isOpen())
    {
        sf::Time deltaTime = deltaClock.restart();
        double setDT = 1.0 / 1000.f;
        double trueDT = deltaTime.asSeconds();

        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            else if (event.key.code == sf::Keyboard::Q) {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::D) {
                    SDF.changeDrawSolidCells();
                }
            }
        }

        window.clear();

        SDF.update();

        frame++;
        if (frame == 20) {
            fps = (int)(1.f / trueDT);
            frame = 0;
        }

        window.display();
       
    }

    return 0;
}