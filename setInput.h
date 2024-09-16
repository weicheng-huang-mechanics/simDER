#ifndef SETINPUT_H
#define SETINPUT_H

#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "Option.h"
#include "eigenIncludes.h"

class setInput
{
public:

  typedef std::map<std::string, Option> OptionMap;
  OptionMap m_options;

  setInput();
  ~setInput();

  template <typename T>
	int AddOption(const std::string& name, const std::string& desc, const T& def);

  Option* GetOption(const std::string& name);
  bool& GetBoolOpt(const std::string& name);
  int& GetIntOpt(const std::string& name);
  double& GetScalarOpt(const std::string& name);
  Vector3d& GetVecOpt(const std::string& name);
  string& GetStringOpt(const std::string& name);

  int LoadOptions(const char* filename);
  int LoadOptions(const std::string& filename)
  {
    return LoadOptions(filename.c_str());
  }
  int LoadOptions(int argc, char** argv);

private:
	double RodLength;
	double rodRadius;
	double helixradius;
	double helixpitch;
	int numVertices;
	double youngM;
	double Poisson;
	double deltaTime;
	double totalTime;
	double tol;
	double stol;
	int maxIter; // maximum number of iterations
	double density;
	Vector3d gVector;
	Vector3d bVector;
	double viscosity;
	bool render;
	bool saveData;
};

#include "setInput.tcc"

#endif // PROBLEMBASE_H
