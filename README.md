# cbstools
CBS High-Res Brain Processing Tools

This repository contains the core elements of the CBS Tools software developed by Pierre-Louis Bazin and his team 
at the Max Planck Institute for Human Cognitive and Brain Sciences in Leipzig, Germany, starting with version 3.1. 
Other contributors may add modules in separate repositories.


Repository overview

The CBS Tools include software modules for the JIST and MIPAV environment as well as python classes.
Source Java code is located under de/mpg/cbs/ and organized as described below.
Python wrappers and python compiling scripts are located under python/.

The CBS Tools software comes with atlas data under atlases/, and standard processing pipeline layouts for JIST under layouts/
Sample data sets for testing purposes are provided under data/


Organization of the CBS Tools code base

This details general organizing rules for the code within the CBS Tools.
Please follow these rules if you want to contribute.

The Java code for the CBS Tools is organized under de/mpg/cbs as follows:
- core: contains the core algorithms for data processing, to be called either as a JIST module or an encapsulated python class.
- jist: contains all the JIST module definitions, including algorithms not encapsulated in pyhon.

Both core and jist have sub-sections corresponding to algorithm categories. 
New sub-sections may only be added at that level, and require justification.
Both core and jist have a rigid i/o structure to allow proper interfacing with JIST and python.

- libraries: contains simple image and data processing libraries. 
Static classes (i.e.without internal data structures) are encouraged.
Any application-specific algorithm or function should not be included here.
No sub-sections.

- methods: contains the methods and classes used in algorithms.
There are no restrictions on the structuring of classes within this section.
No sub-sections.

- structures: contains generic data structures for processing.
Application-specific structure types should not be added here.
Structures defined here should be lightweight and self-contained.
No sub-sections.

- utilities: contains convenience classes for i/o, interfacing with JIST and MIPAV, simple numerical operations.
Static classes (i.e.without internal data structures) are encouraged.
Any application-specific algorithm or function should not be included here.
No sub-sections.

