ManBo
=====

ab initio Molecular Dynamics Program that uses Many Body Expansion

Description
===========

ManBo makes use of the Electrostatically Embedded Many-Body Expansion
in order to reduce the computational cost of ab initio Molecular
Dynamics. The aim is to perform that while keeping an acceptable
accuracy on the results.

Installing ManBo
================

In order to install ManBo, clone this git repository and change to
the directory you cloned to. Then, change to the directory that
contains ManBo's source code:

```
cd src
```

There are different compilation options. Each option is present in
a different make file. Choose the make file you want and copy to
Makefile:


```
cp Makefile-[COMPILATION OPTION] Makefile
```

Then execute:

```
make install
```

The executables will be placed at the ../bin/ directory.

Authors and Contributors
========================

The following people have contributed to ManBo's coding development:

* Vinicius Wilian D. Cruzeiro (principal developer) | vwcruzeiro@ufl.edu
* Herbert C. Georg


Some history
============

ManBo started as an undergraduate project at the University of Goiás, Brazil.
