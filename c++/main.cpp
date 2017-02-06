#include "FeynmanDiagram/FeynmanDiagram.h"
// #include "Vertex/Vertex.h"
// #include "Propagator/Electron.h"
// #include "Propagator/Phonon.h"

int main() {
  FeynmanDiagram FD(Vector3d(0,0,0), 1, 2, 3);

  auto v1 = FD.insertVertex(0, 0.1);
  auto v2 = FD.insertVertex(1, 0.1);
  auto v3 = FD.insertVertex(2, 0.1);

  FD.addInternalPhonon(v1, v3, Vector3d(1,2,1), 1, 2);
  // FD.addInternalPhonon(FD.Vs.begin(), prev(FD.Vs.end(), 1), Vector3d(3,4,5), 1, 2);

  FD.print();
}

// #include <iostream>
// #include <memory>
// #include <vector>
// using namespace std;

// class Class1 : public enable_shared_from_this<Class1> {
//   public:
//     int param;
//     Class1 (int p) {this->param = p;};

//     shared_ptr<Class1> getptr() {
//         return shared_from_this();
//     }

//     void print () {
//       cout << this << endl;
//       cout << shared_from_this() << endl;
//     }
// };

// int main() {
//   vector<shared_ptr<Class1>> vec;

//   vec.emplace_back(new Class1 {0});
//   vec.emplace_back(new Class1 {1});


//   shared_ptr<Class1> test = vec[0]->getptr();

//   test->param = 3;

//   // vec[1].reset();
//   // // vec[0].reset();

//   // cout << vec[0] << endl;
//   // cout << vec[1] << endl;
//   // cout << test << endl;


//   vec[0]->print();


//   // cout << vec[0].use_count() << endl;
//   // cout << vec[1].use_count() << endl;
//   // cout << test.use_count() << endl;
// };
