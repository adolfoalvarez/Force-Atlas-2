Force-Atlas-2
=============

R implementation of the Force Atlas 2 graph layout designed for Gephi. The algorithm is detailed in:

[Jacomy M, Venturini T, Heymann S, Bastian M (2014) ForceAtlas2, a Continuous Graph Layout Algorithm for Handy Network Visualization Designed for the Gephi Software. PLoS ONE 9(6): e98679. doi:10.1371/journal.pone.0098679](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098679)


# Installation

This is not a package yet, so you should download the "ForceAtlas2.R" file and source() it.

# Usage
 ```R
    library(gephi)
    g <- graph.ring(100)
    layout.forceatlas2(g, iterations=10000, plotstep=500)
  ```

# Authors

- Bazyli Klockiewicz (bazyli.klockiewicz@analyx.com)
- Adolfo Alvarez (adolfo.alvarez@analyx.com)

