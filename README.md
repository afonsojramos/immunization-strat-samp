## Immunization Strategies and Sampling - A Complex and Social Networks Project

### Requirements

This project will explore several strategies for immunization in real and synthetic (big) networks, but for which we do not (or can not) have full knowledge, hence must combine immunisation with some crawling strategy to sample part of the network. The problem is as follows: a disease spreads in a network according to the SIS model. We are allowed to vaccinate a limited number of "nodes" -- what this means is that these vaccinated nodes will never contract the disease. A typical technique is to direct vaccination to those nodes of highest degree. This project will explore other immunization strategies such as directing vaccination to those nodes of highest pagerank or betweenness centrality.
In addition, since we do not have full knowledge of the graph we will crawl part of the graph and determine which nodes to vaccinate during the crawl (using some reasonable heuristic). Then, we will see which crawling strategy leads to better results in terms of total fraction of infecteds after a finite amount of time.

In particular, I believe it is interesting to check whether the crawling strategy of "expansion sampling" [1] beats other strategies, such as a random walks (Random node or random edge sampling).

The result of this project should be an empirical study of what strategies are better suited for different network models, networks with different degree distribution, etc.
(An alternative to sampling is to use one of available web crawls: http://webdatacommons.org/hyperlinkgraph/ .
Some immunization strategies that you can try out and make cost-benefit analysis can be found in
https://snap.stanford.edu/class/cs224w-readings/pastorsatorras02immunization.pdf )

[1] Maiya, A. S. and Berger-Wolf, T. Y. (2010). Sampling community structure. In Proceedings of the 19th International Conference on World Wide Web, WWW’10, pages 701–710, New York, NY, USA. ACM.
