# ngl
Neighbor Graph Library

NGL is a lightweight library to compute a variety of geometric neighborhood
graphs in arbitrary dimensions, including the Relative Neighbor graph, the
Gabriel graph, and the beta skeleton.

It is a re-write of the NGL library available in ngraph.org, and implements
the algorithms described here: 
Carlos D. Correa and Peter Lindstorm, "Towards Robust Topology of Sparsely
Sampled Data". IEEE Transactions on Visualization and Computer Graphics
(Proceedings Visualization / Information Visualization 2011), vol. 17, no. 12,
Dec. 2011.

## Building

```
cd ngl/src
make all
cd ../test
make all
```

## Example

```
#include "ngl.h"

void computeRelativeNeighborGraph(const vector<double*>& points, int dims) {
  // Create a vector space of double precision
  ngl::DoubleVectorSpace space(dims);

  // Create a relative neighbor graph builder
  ngl::RelativeNeighborGraphBuilderD builder(space);
  // Add points to the builder
  builder.addPoints(points);

  ngl::NeighborGraphImpl neighborGraph;
  // Compute neighborhood graph
  builder.computeNeighborGraph(&neighborGraph);

  // Write edges of the graph to stdout
  for (int i = 0; i < neighborGraph.size(); ++i) {
    int src;
    int dst;
    neighborGraph.getEdge(i, &src, &dst);
    printf("%d %d\n", src, dst);
  }
}
```