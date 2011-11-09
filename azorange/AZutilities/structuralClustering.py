""" Module for the integration of the TUM structural clustering algorithm
	author: Tobias Girschick; tobias.girschick@in.tum.de
			TUM - I12 (wwwkramer.in.tum.de/girschic)
	dependencies: gSpan
	
	Please cite the following article if you use the structural clustering procedure or results produced with it in any publication:
	
@InProceedings{ seeland2010,
	title = "Online Structural Graph Clustering Using Frequent Subgraph Mining",
	booktitle = "Proceedings of the ECML/PKDD'10",
	series = "Lecture Notes in Computer Science",
	author = "M. Seeland and T. Girschick and F. Buchwald and S. Kramer",
	editor = "J.L. Balc{\'a}zar and F. Bonchi and A. Gionis and M. Sebag",
	pages = "213--228",
	volume = "6323",
	year = "2010",
	booktitle1 = "Machine Learning and Knowledge Discovery in Databases, European Conference, {ECML} {PKDD} 2010, Barcelona, Spain, September 20-24, 2010, Proceedings, Part {III}"
}
"""

import os
import sys
import time

import subprocess
from subprocess import Popen, PIPE
import tempfile
