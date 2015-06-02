#include "io-utils.h"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using std::ifstream;
using std::runtime_error;
using std::string;
using std::vector;

void readPoints(const char *src, int dims, vector<double>* points) {
  ifstream inp(src);
  if (!inp.is_open()) {
    fprintf(stderr, "Cannot open file %s\n", src);
    throw runtime_error("Cannot open file");
  }
  string line;
  while (!inp.eof()) {
    getline(inp, line);
    if (line.length() == 0) continue;
    const char *sline = line.c_str();
    char *tok = strtok((char*) sline, " \t\r\n");
    if (tok && tok[0] == '#') continue;
    int ndims = 0;
    vector<double> readPts;
    while (tok && ndims < dims) {
      readPts.push_back(atof(tok));
      points->push_back(atof(tok));
      tok = strtok(NULL, " \t\r\n");
      ndims++;
    }
    if (ndims < dims) {
      fprintf(stderr, "Error reading input file: num.dims=%d, expecting %d\n",
          ndims, dims);
      throw runtime_error("Error reading input file");
    }
  }
  inp.close();
}

