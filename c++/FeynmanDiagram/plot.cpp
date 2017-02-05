#include "FeynmanDiagram.h"

void FeynmanDiagram::drawPhonon (sf::RenderWindow *window, double yPos, double start, double end) {
  double radius = 0.5*(end - start - this->thickness);

  sf::CircleShape circle(radius, (int) 2*radius);
  circle.setOrigin(-0.5*this->thickness, radius);
  circle.setPosition(this->horizontalMargin + start, this->verticalMargin + yPos);

  circle.setFillColor(sf::Color(255, 255, 255, 0));
  circle.setOutlineThickness(this->thickness);
  circle.setOutlineColor(sf::Color(255, 0, 0));

  window->draw(circle);

  sf::RectangleShape rectangle(sf::Vector2f(2*(radius + this->thickness), radius + this->thickness));
  rectangle.setFillColor(sf::Color(255, 255, 255));
  rectangle.setPosition(this->horizontalMargin + start - 0.5*this->thickness, this->verticalMargin + yPos);

  window->draw(rectangle);
}

void FeynmanDiagram::drawVertex (sf::RenderWindow *window, double yPos, double center, bool solid) {
  sf::CircleShape circle(this->thickness, (int) 2*thickness);
  circle.setOrigin(this->thickness, this->thickness);
  circle.setPosition(this->horizontalMargin + center, this->verticalMargin + yPos);

  // gray or black vertex
  if (solid) {
    circle.setFillColor(sf::Color(0, 0, 0));
  } else {
    circle.setFillColor(sf::Color(0.5*255, 0.5*255, 0.5*255));
  }

  window->draw(circle);
}

void FeynmanDiagram::drawElectron (sf::RenderWindow *window, double yPos, double start, double end) {
  sf::RectangleShape rectangle(sf::Vector2f(end - start, thickness));
  rectangle.setFillColor(sf::Color(0, 255, 255));
  rectangle.setPosition(this->horizontalMargin + start, this->verticalMargin + yPos - 0.5*thickness);

  window->draw(rectangle);
}

int FeynmanDiagram::longestPhonon () {
  int counter, longestPhonon = 0;

  for (list<Vertex>::iterator i = this->Vs.begin(); i != this->Vs.end(); ++i) {
    if (i->D[1]) {

      counter = 0;
      for (list<Vertex>::iterator j = next(i, 1); j != this->Vs.end(); ++j) {
        counter++;

        if (i->D[1]->end == &(*j) && counter > longestPhonon) {
          longestPhonon = counter;
          break;
        }
      }
    }
  }

  return longestPhonon;
}

void FeynmanDiagram::plot () {
  this->windowWidth = 2000;
  this->windowHeight = 1000;
  this->horizontalMargin = 20;
  this->verticalMargin = 20;
  this->thickness = 10;

  // antialiasing
  sf::ContextSettings settings;
  settings.antialiasingLevel = 8;

  sf::RenderWindow window(
    sf::VideoMode(this->windowWidth + 2*this->horizontalMargin, this->windowHeight + 2*this->verticalMargin),
    "Feynman Diagram",
    sf::Style::Default,
    settings
  );

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


      int numVertices, longestPhonon;
      numVertices = this->Vs.size();
      longestPhonon = this->longestPhonon();

      double diagramRatio, windowRation;
      diagramRatio = (numVertices - 1)/(0.5*this->longestPhonon());
      windowRation = this->windowWidth/this->windowHeight;

      double gLength, left, right, top, bottom;
      if (diagramRatio >= windowRation) {
        // propagator lenth according to window length
        gLength = this->windowWidth/(numVertices - 1);
        left = 0;
        right = this->windowWidth;
        top = 0.5*this->windowHeight - 0.25*longestPhonon*gLength;
        bottom = 0.5*this->windowHeight + 0.25*longestPhonon*gLength;
      } else {
        // propagator lenth according to window length
        gLength = this->windowHeight/(0.5*longestPhonon);
        left = 0.5*this->windowWidth - 0.5*(numVertices - 1)*gLength;
        right = this->windowWidth + 0.5*(numVertices - 1)*gLength;
        top = 0;
        bottom = this->windowHeight;
      }

      cout << gLength << endl;
      cout << left << " "
           << right << " "
           << top << " "
           << bottom << endl;

      window.clear(sf::Color(255, 255, 255));


      // this->drawPhonon(&window, bottom, 0, 970);
      // this->drawPhonon(&window, bottom, 900, 1100);

      // this->drawElectron(&window, bottom, 900, 970);

      // this->drawVertex(&window, bottom, 900);
      // this->drawVertex(&window, bottom, 970);


      //
      // phonons first, electrons second and vertices last
      //
      int start = 0, end;
      for (list<Vertex>::iterator i = this->Vs.begin(); i != this->Vs.end(); ++i) {
        if (i->D[1]) {

          end = start + 1;
          for (list<Vertex>::iterator j = next(i, 1); j != this->Vs.end(); ++j) {
            if (i->D[1]->end == &(*j)) {
              this->drawPhonon(&window, bottom, start*gLength, end*gLength);
              break;
            }

            end++;
          }
        }

        start++;
      }

      int count = 0;
      for (list<Vertex>::iterator i = this->Vs.begin(); i != prev(this->Vs.end(), 1); ++i) {
        this->drawElectron(&window, bottom, count*gLength, (count + 1)*gLength);

        this->drawVertex(&window, bottom, count*gLength, (bool) count);
        count++;
      }
      this->drawVertex(&window, bottom, count*gLength, false);



      //
      // temp outline
      //
      sf::Vertex corners[] = {
          sf::Vertex(sf::Vector2f(this->horizontalMargin + left, this->verticalMargin + top), sf::Color(255, 0, 255)),
          sf::Vertex(sf::Vector2f(this->horizontalMargin + right, this->verticalMargin + top), sf::Color(255, 0, 255)),
          sf::Vertex(sf::Vector2f(this->horizontalMargin + right, this->verticalMargin + bottom), sf::Color(255, 0, 255)),
          sf::Vertex(sf::Vector2f(this->horizontalMargin + left, this->verticalMargin + bottom), sf::Color(255, 0, 255)),
          sf::Vertex(sf::Vector2f(this->horizontalMargin + left, this->verticalMargin + top), sf::Color(255, 0, 255))
      };
      window.draw(corners, 5, sf::LinesStrip);

      window.display();
    }
  }
}