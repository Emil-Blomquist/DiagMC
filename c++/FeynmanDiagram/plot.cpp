// #include "FeynmanDiagram.h"

// void FeynmanDiagram::drawElectron (sf::RenderWindow *window, double yPos, double start, double end, sf::Color color) {
//   sf::RectangleShape rectangle(sf::Vector2f(end - start, thickness));
//   rectangle.setFillColor(color);
//   rectangle.setPosition(this->horizontalMargin + start, this->verticalMargin + yPos - 0.5*thickness);

//   window->draw(rectangle);
// }

// void FeynmanDiagram::drawPhonon (sf::RenderWindow *window, double yPos, double start, double end, sf::Color color) {
//   double radius = 0.5*(end - start - this->thickness);

//   sf::CircleShape circle(radius, (int) 2*radius);
//   circle.setOrigin(-0.5*this->thickness, radius);
//   circle.setPosition(this->horizontalMargin + start, this->verticalMargin + yPos);

//   circle.setFillColor(sf::Color(255, 255, 255, 0));
//   circle.setOutlineThickness(this->thickness);
//   circle.setOutlineColor(color);

//   window->draw(circle);

//   sf::RectangleShape rectangle(sf::Vector2f(2*(radius + this->thickness), radius + this->thickness));
//   rectangle.setFillColor(sf::Color(255, 255, 255));
//   rectangle.setPosition(this->horizontalMargin + start - 0.5*this->thickness, this->verticalMargin + yPos);

//   window->draw(rectangle);
// }

// void FeynmanDiagram::drawVertex (sf::RenderWindow *window, double yPos, double center, bool solid) {
//   sf::CircleShape circle(this->thickness, (int) 2*thickness);
//   circle.setPosition(this->horizontalMargin + center, this->verticalMargin + yPos);

//   if (solid) {
//     circle.setFillColor(sf::Color(0, 0, 0));
//     circle.setOrigin(this->thickness, this->thickness);
//   } else {
//     circle.setFillColor(sf::Color(255, 255, 255));
//     circle.setRadius(this->thickness - 2);
//     circle.setOutlineThickness(2);
//     circle.setOutlineColor(sf::Color(0, 0, 0));
//     circle.setOrigin(this->thickness - 2, this->thickness - 2);
//   }

//   window->draw(circle);
// }

// void FeynmanDiagram::drawText (sf::RenderWindow *window, sf::Font* font, double toStr, double xPos, double yPos, double aboveOrBelow) {
//   int fontSize = 30;

//   // round to three decimals
//   toStr = round(toStr*1000)/1000;

//   // convert to string
//   ostringstream strs;
//   strs << toStr;
//   string str = strs.str();

//   sf::Text text;
//   text.setFont(*font);
//   text.setCharacterSize(fontSize);
//   text.setFillColor(sf::Color(0, 0, 0));
//   text.setString(str);

//   // centering
//   sf::FloatRect textRect = text.getLocalBounds();
//   text.setOrigin(textRect.width/2,textRect.height/2 + fontSize*0.3);
//   text.setPosition(this->horizontalMargin + xPos, this->verticalMargin + yPos + (0.5*fontSize + this->thickness)*aboveOrBelow);

//   window->draw(text);
// }

// int FeynmanDiagram::longestPhonon () {
//   int counter, longestPhonon = 0;

//   for (list<Vertex>::iterator i = this->Vs.begin(); i != this->Vs.end(); ++i) {
//     if (i->D[1]) {

//       counter = 0;
//       for (list<Vertex>::iterator j = next(i, 1); j != this->Vs.end(); ++j) {
//         counter++;

//         if (i->D[1]->end == &(*j) && counter > longestPhonon) {
//           longestPhonon = counter;
//           break;
//         }
//       }
//     }
//   }

//   return longestPhonon;
// }

// sf::Color FeynmanDiagram::colorCode (double pMin, double pMax, Vector3d P) {
//   if (pMin == pMax) {
//     return sf::Color(0, 0, 255);
//   } else {
//     double r = (P.norm() - pMin)/(pMax - pMin);
//     return sf::Color(255*r, 0, 255*(1 - r));
//   }
// }

// void FeynmanDiagram::plot () {
//   this->windowWidth = 2000;
//   this->windowHeight = 1000;
//   this->horizontalMargin = 100;
//   this->verticalMargin = 100;
//   this->thickness = 10;

//   // antialiasing
//   sf::ContextSettings settings;
//   settings.antialiasingLevel = 8;

//   sf::RenderWindow window(
//     sf::VideoMode(this->windowWidth + 2*this->horizontalMargin, this->windowHeight + 2*this->verticalMargin),
//     "Feynman Diagram",
//     sf::Style::Default,
//     settings
//   );

//   window.setFramerateLimit(10);

//   //
//   // prepare for drawing Feynman diagram
//   //
//   int numVertices, longestPhonon;
//   numVertices = this->Vs.size();
//   longestPhonon = this->longestPhonon();

//   double diagramRatio, windowRation;
//   diagramRatio = (numVertices - 1)/(0.5*this->longestPhonon());
//   windowRation = this->windowWidth/this->windowHeight;

//   double gLength, left, right, top, bottom;
//   if (diagramRatio >= windowRation) {
//     // propagator lenth according to window length
//     gLength = this->windowWidth/(numVertices - 1);
//     left = 0;
//     right = this->windowWidth;
//     top = 0.5*this->windowHeight - 0.25*longestPhonon*gLength;
//     bottom = 0.5*this->windowHeight + 0.25*longestPhonon*gLength;
//   } else {
//     // propagator lenth according to window length
//     gLength = this->windowHeight/(0.5*longestPhonon);
//     left = 0.5*this->windowWidth - 0.5*(numVertices - 1)*gLength;
//     right = this->windowWidth + 0.5*(numVertices - 1)*gLength;
//     top = 0;
//     bottom = this->windowHeight;
//   }

//   // 
//   // declare font
//   //
//   sf::Font font;
//   font.loadFromFile("OpenSans-Regular.ttf");

//   // 
//   // find max and min momentum
//   //
//   double p = this->Vs.begin()->G[1]->momentum.norm(), pMin = p, pMax = p;
//   for (list<Vertex>::iterator i = this->Vs.begin(); i != prev(this->Vs.end(), 1); ++i) {
//     p = i->G[1]->momentum.norm();
//     pMin = min(pMin, p);
//     pMax = max(pMax, p);

//     if (i->D[1]) {
//       p = i->D[1]->momentum.norm();
//       pMin = min(pMin, p);
//       pMax = max(pMax, p);
//     }
//   }

//   //
//   // rendering loop
//   //
//   int rendered = false;
//   while (window.isOpen()) {
//     sf::Event event;
//     while (window.pollEvent(event)) {
//       if (event.type == sf::Event::Closed) {
//         window.close();
//       }
//     }

//     if (!rendered) {
//       rendered = true;

//       window.clear(sf::Color(255, 255, 255));

//       //
//       // phonons first, electrons second and vertices last
//       //
//       int start = 0, end;
//       for (list<Vertex>::iterator i = this->Vs.begin(); i != this->Vs.end(); ++i) {
//         if (i->D[1]) {

//           end = start + 1;
//           for (list<Vertex>::iterator j = next(i, 1); j != this->Vs.end(); ++j) {
//             if (i->D[1]->end == &(*j)) {
//               this->drawPhonon(&window, bottom, start*gLength, end*gLength, this->colorCode(pMin, pMax, i->D[1]->momentum));
//               this->drawText(&window, &font, i->D[1]->momentum.norm(), 0.5*(start + end)*gLength, bottom - 0.5*(end - start)*gLength, -1);
//               break;
//             }

//             end++;
//           }
//         }

//         start++;
//       }

//       int count = 0;
//       for (list<Vertex>::iterator i = this->Vs.begin(); i != prev(this->Vs.end(), 1); ++i) {
//         this->drawElectron(&window, bottom, count*gLength, (count + 1)*gLength, this->colorCode(pMin, pMax, i->G[1]->momentum));
//         this->drawText(&window, &font, i->G[1]->momentum.norm(), (count + 0.5)*gLength, bottom, -1);

//         count++;
//       }

//       count = 0;
//       bool solidVertex;
//       for (list<Vertex>::iterator i = this->Vs.begin(); i != this->Vs.end(); ++i) {
        
//         solidVertex = true;
//         if (i == this->Vs.begin() || i == prev(this->Vs.begin(), 2)) {
//           solidVertex = false;
//         }

//         this->drawVertex(&window, bottom, count*gLength, solidVertex);
//         this->drawText(&window, &font, i->position, count*gLength, bottom, 1);

//         count++;
//       }

//       window.display();
//     }
//   }
// }