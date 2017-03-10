// #include "Display.h"

// Display::Display (FeynmanDiagram *FD) {
//   this->FD = FD;

//   // settings
//   this->windowWidth = 2000;
//   this->windowHeight = 1000;
//   this->horizontalMargin = 100;
//   this->verticalMargin = 100;
//   this->thickness = 10;
// }

// void Display::render () {
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
//   int numElectrons, longestPhonon;
//   numElectrons = this->FD->Gs.size();
//   longestPhonon = this->longestPhonon();

//   double diagramRatio, windowRation;
//   diagramRatio = numElectrons/(0.5*this->longestPhonon());
//   windowRation = this->windowWidth/this->windowHeight;

//   double gLength, left, right, top, bottom;
//   if (diagramRatio >= windowRation) {
//     // propagator lenth according to window length
//     gLength = this->windowWidth/numElectrons;
//     left = 0;
//     right = this->windowWidth;
//     top = 0.5*this->windowHeight - 0.25*longestPhonon*gLength;
//     bottom = 0.5*this->windowHeight + 0.25*longestPhonon*gLength;
//   } else {
//     // propagator lenth according to window length
//     gLength = this->windowHeight/(0.5*longestPhonon);
//     left = 0.5*this->windowWidth - 0.5*numElectrons*gLength;
//     right = this->windowWidth + 0.5*numElectrons*gLength;
//     top = 0;
//     bottom = this->windowHeight;
//   }

//   // 
//   // declare font
//   //
//   sf::Font font;
//   font.loadFromFile("Display/OpenSans-Regular.ttf");

//   // 
//   // find max and min momentum
//   //
//   double p = this->FD->Gs[0]->momentum.norm(), pMin = p, pMax = p;
//   for (auto g : this->FD->Gs) {
//     p = g->momentum.norm();
//     pMin = min(pMin, p);
//     pMax = max(pMax, p);

//     if (g->start->D[1]) {
//       p = g->start->D[1]->momentum.norm();
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

//       shared_ptr<Phonon> d;
//       for (int i = 0; i < this->FD->Gs.size(); i++) {
//         d = this->FD->Gs[i]->start->D[1];
//         if (d) {

//           for (int j = i; j < this->FD->Gs.size(); j++) {
//             if (d->end == this->FD->Gs[j]->end) {
//               this->drawPhonon(&window, bottom, i*gLength, (j + 1)*gLength, this->colorCode(pMin, pMax, d->momentum));
//               this->drawText(&window, &font, d->momentum.norm(), 0.5*(i + j + 1)*gLength, bottom - 0.5*(j - i + 1)*gLength, -1);
//               break;
//             }
//           }

//         }
//       }

//       shared_ptr<Electron> g;
//       for (int i = 0; i < this->FD->Gs.size(); i++) {
//         g = this->FD->Gs[i];

//         this->drawElectron(&window, bottom, i*gLength, (i + 1)*gLength, this->colorCode(pMin, pMax, g->momentum));
//         this->drawText(&window, &font, g->momentum.norm(), (i + 0.5)*gLength, bottom, -1);
//       }

//       // 
//       // validate vertices
//       //
//       shared_ptr<Vertex> v;
//       Vector3d dP;
//       bool validVertex;
//       for (int i = 1; i < this->FD->Gs.size(); i++) {
//         validVertex = true;
//         v = this->FD->Gs[i]->start;

//         // three lines connected
//         if ( ! (v->G[0] && v->G[1] && ((v->D[0] != NULL) ^ (v->D[1] != NULL)))) {
//           cout << i << ": DIAGRAM ERROR: Not exactly three lines connected to a vertex" << endl;
//           v->print();
//           validVertex = false;
//         }

//         // momentum conservation
//         dP = v->G[1]->momentum - v->G[0]->momentum;
//         if (v->D[1]) {
//           dP += v->D[1]->momentum;
//         } else if (v->D[0]) {
//           dP -= v->D[0]->momentum;
//         }
//         if (dP.norm() > pow(10.0, -10.0)) {
//           cout << i << ": DIAGRAM ERROR: Momentum nonconservation at vertex" << endl;
//           v->print();
//           validVertex = false;
//         }

//         // increasing vertex times
//         if ( ! (v->G[0]->start->position < v->position && v->position < v->G[1]->end->position)) {
//           cout << i << ": DIAGRAM ERROR: Vertices not chronologically ordered" << endl;
//           v->print();
//           validVertex = false;
//         }

//         // plot invalid vertex
//         if ( ! validVertex) {
//           this->markVertex(&window, i*gLength, bottom);
//         }
//       }

//       for (int i = 0; i < this->FD->Gs.size(); i++) {
//         v = this->FD->Gs[i]->end;

//         this->drawVertex(&window, (i + 1)*gLength, bottom, i < this->FD->Gs.size() - 1);
//         this->drawText(&window, &font, v->position, (i + 1)*gLength, bottom, 1);
//       }
//       this->drawVertex(&window, 0, bottom, false);
//       this->drawText(&window, &font, 0, 0, bottom, 1);

//       //
//       // validate phonon propagators
//       //
//       double validPhonon;
//       for (auto d : this->FD->Ds) {
//         validPhonon = true;

//         if (d->end->position <= d->start->position) {
//           cout << "DIAGRAM ERROR: Phonon propagating backwards" << endl;
//           d->print();
//           validPhonon = false;
//         }
//       }

//       window.display();
//     }
//   }
// }




















// //
// // methods used in Display::render
// //
// void Display::drawElectron (sf::RenderWindow *window, double yPos, double start, double end, sf::Color color) {
//   sf::RectangleShape rectangle(sf::Vector2f(end - start, thickness));
//   rectangle.setFillColor(color);
//   rectangle.setPosition(this->horizontalMargin + start, this->verticalMargin + yPos - 0.5*thickness);

//   window->draw(rectangle);
// }

// void Display::drawPhonon (sf::RenderWindow *window, double yPos, double start, double end, sf::Color color) {
//   double radius = abs(0.5*(end - start - this->thickness));

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

// void Display::drawVertex (sf::RenderWindow *window, double xPos, double yPos, bool solid) {
//   sf::CircleShape circle(this->thickness, (int) 2*thickness);
//   circle.setPosition(this->horizontalMargin + xPos, this->verticalMargin + yPos);

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

// void Display::markVertex (sf::RenderWindow *window, double xPos, double yPos) {
//   sf::CircleShape circle(this->thickness, (int) 2*thickness);
//   circle.setOrigin(this->thickness, this->thickness);
//   circle.setPosition(this->horizontalMargin + xPos, this->verticalMargin + yPos);

//   circle.setOutlineThickness(3*this->thickness);
//   circle.setFillColor(sf::Color(255, 255, 255, 0));
//   circle.setOutlineColor(sf::Color(42, 252, 92));

//   window->draw(circle);
// }

// void Display::drawText (sf::RenderWindow *window, sf::Font* font, double toStr, double xPos, double yPos, double aboveOrBelow) {
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

// int Display::longestPhonon () {
//   int longestPhonon = 0;

//   for (int i = 0; i < this->FD->Gs.size(); i++) {
//     if (this->FD->Gs[i]->start->D[1]) {

//       for (int j = i; j < this->FD->Gs.size(); j++) {
//         if (
//           this->FD->Gs[i]->start->D[1]->end == this->FD->Gs[j]->end &&
//           j - i + 1 > longestPhonon
//         ) {
//           longestPhonon = j - i + 1;
//           break;
//         }
//       }

//     }
//   }

//   return longestPhonon;
// }

// sf::Color Display::colorCode (double pMin, double pMax, Vector3d P) {
//   if (pMin == pMax) {
//     return sf::Color(0, 0, 255);
//   } else {
//     double r = (P.norm() - pMin)/(pMax - pMin);
//     return sf::Color(255*r, 0, 255*(1 - r));
//   }
// }