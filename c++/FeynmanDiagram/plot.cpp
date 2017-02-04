#include "FeynmanDiagram.h"

void FeynmanDiagram::drawPhonon (sf::RenderWindow *window, int start, int end) {
  double radius = 0.5*(end - start - this->thickness);
  double middle = 0.5*this->height;

  sf::CircleShape circle(radius, (int) 2*radius);
  circle.setOrigin(-0.5*this->thickness, radius);
  circle.setPosition(start, middle);

  circle.setFillColor(sf::Color(255, 255, 255, 0));
  circle.setOutlineThickness(this->thickness);
  circle.setOutlineColor(sf::Color(255, 0, 0));

  window->draw(circle);

  sf::RectangleShape rectangle(sf::Vector2f(2*(radius + this->thickness), radius + this->thickness));
  rectangle.setFillColor(sf::Color(255, 255, 255));
  rectangle.setPosition(start - 0.5*this->thickness, middle);

  window->draw(rectangle);
}

void FeynmanDiagram::drawVertex (sf::RenderWindow *window, int center) {
  sf::CircleShape circle(this->thickness, (int) 2*thickness);
  circle.setOrigin(this->thickness, this->thickness);
  circle.setPosition(center, 0.5*this->height);
  circle.setFillColor(sf::Color(0, 0, 0));

  window->draw(circle);
}

void FeynmanDiagram::drawElectron (sf::RenderWindow *window, int start, int end) {
  sf::RectangleShape rectangle(sf::Vector2f(end - start, thickness));
  rectangle.setFillColor(sf::Color(0, 255, 255));
  rectangle.setPosition(start, 0.5*this->height - 0.5*thickness);

  window->draw(rectangle);
}

void FeynmanDiagram::plot () {
  this->width = 2000;
  this->height = 1000;
  this->thickness = 10;

  // antialiasing
  sf::ContextSettings settings;
  settings.antialiasingLevel = 8;

  sf::RenderWindow window(sf::VideoMode(this->width, this->height), "SFML works!", sf::Style::Default, settings);
  window.setFramerateLimit(10);

  int rendered = false;
  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) {
        window.close();
      }
    }

    if (!rendered) {
      rendered = true;

      window.clear(sf::Color(255, 255, 255));


      this->drawPhonon(&window, 0, 970);
      this->drawPhonon(&window, 900, 1100);

      this->drawVertex(&window, 900);

      this->drawElectron(&window, 900, 970);


      sf::Vertex line[] = {
          sf::Vertex(sf::Vector2f(0, 500), sf::Color(255, 0, 255)),
          sf::Vertex(sf::Vector2f(2000, 500), sf::Color(255, 0, 255))
      };
      window.draw(line, 2, sf::Lines);

      window.display();
    }
  }
}