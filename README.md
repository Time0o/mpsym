[![GitHub Releases](https://img.shields.io/github/release/Time0o/mpsym.svg)](https://github.com/Time0o/mpsym/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-doxygen-blue.svg)](https://Time0o.github.io/mpsym)
[![Build Status](https://travis-ci.com/Time0o/mpsym.svg?branch=master)](https://travis-ci.com/Time0o/mpsym)
[![Coverage
status](https://codecov.io/gh/Time0o/mpsym/branch/master/graph/badge.svg)](https://codecov.io/gh/Time0o/mpsym?branch=master)

<p align="center">
  <img src="image.png">
</p>

# MPsym - Multiprocessor System  Symmetry Reduction

## Contents

- [Introduction](#introduction)
- [A Motivating Example](#a-motivating-example)
- [How It Works](#how-it-works)
- [Installation](#installation)
- [Examples](#references)
  - [Defining Architecture Graphs](#defining-architecture-graphs)
  - [Hierarchical Architecture Graphs](#hierarchical-architecture-graphs)
  - [Persisting Architecture Graphs](#persisting-architecture-graphs)
  - [Initializing Architecture Graphs](#initializing-architecture-graphs)
  - [Orbits and Representatives](#orbits-and-representatives)
  - [Automorphism Groups](#automorphism-groups)
- [Limitations](#limitations)
- [Developer Notes](#developer-notes)
  - [Main Project Structure](#main-project-structure)
  - [Profiling](#profiling)
  - [Continuous Integration](#continuous-integration)
- [References](#references)

## Introduction

_MPsym_ is a C++/Python library that makes it possible to determine whether
mappings of computational tasks to multiprocessor systems are _equivalent by
symmetry_. This is useful when trying to find an optimal/good (with respect to
runtime/energy consumption etc.) mapping or set of mappings to a given system.
MPsym is able to implicitly partition the search space of all possible mappings
into _equivalence classes_ of mappings that have almost identical runtime
properties due to architectural symmetries. A search algorithm can then work
with  _representatives_ of these equivalence classes, effectively reducing the
size of the search space ifself.

MPsym uses algorithms and data structures from _computational group theory_ as
presented in e.g. [[1]](#1).  These are implemented from scratch to avoid
reliance on existing computational algebra systems like _GAP_ [[2]](#2). As a
result, MPsym can be more easily integrated with other C++ and Python
applications and could be released under the liberal [MIT License](LICENSE.txt).

Although MPsym has been developed in the context of multiprocessor systems, it
is potentially applicable to other problems involving symmetries (more
formally: automorphism groups) of graphs.

## A Motivating Example

As an introductory example, consider an abstract _Multiprocessor System-on-Chip
(MPSoC)_ architecture consisting of four _processing elements_ connected
circularly by four bidirectional _communication channels_ (which could
represent direct wired/wireless links, shared memory etc.). We can represent
this architecture by the following _architecture graph_:

<p align="center">
  <img src="./.img/svg/example_1.svg" />
</p>

Every vertex corresponds to a processing element and each edge corresponds to a
communication channel. In general, such architecture graphs might also be
_vertex-_ and or _edge-colored_ to reflect non-identical processing elements
and communication channels.

Say we want to map two tasks _T_red_ and _T_green_ to this architecture.
Assuming that we don't wish to map both of them to the same processing element,
there are twelve distinct ways of doing so:

<p align="center">
  <img src="./.img/svg/example_2.svg" />
</p>

We refer to this set of all possible mappings for a given set of tasks as the
_full mapping space_. MPsym is able to (implicitly) partition the full mapping
space into sets of mappings equivalent by symmetry:

<p align="center">
  <img src="./.img/svg/example_3.svg" />
</p>

We refer to such sets of mappings equivalent by symmetry as _orbits_. MPsym
_collapses_ the full mapping space by reducing each orbit to one arbitrary
mapping contained in it. If we represent mappings by tuples of processor
indices, it is natural to choose the lexicographically smallest such mapping
(as shown above) which we call the _canonical_ representative.

Since explicitly enumerating all orbits for a given mapping space is often
prohibitively expensive, MPsym offers functions for directly determining the
canonical representative for any given mapping.

## How It Works

For simple architectures, MPsym performs the following steps in order to
determine canonical representatives for a number of mappings:

1. Parse an architecture configuration file.
2. Construct a totally colored architecture graph.
3. Determine that architecture graph's _automorphism group_ using _nauty_ [[3]](#3).
4. Construct a _base and strong generating set_ representation for this group.
5. Find the canonical representative for a given mapping by:
    * Enumerating the orbit.
    * Enumerating the automorphism group.
    * Using local search.

The automorphism group of an architecture graph captures all of the graphs
inherent symmetries in the form of permutations.  For certain classes of
[_hierarchical_ architectures](#hierarchical-architecture-graphs), MPsym uses
the methods presented in [[4]](#4) to speed up this process, making it viable
for large but highly symmetrical architectures.

## Installation

This section explains how to install MPsym to your system. Note that MPsym
currently only runs on Linux. It should in principle be possible to build it on
non-Linux systems but this is currently either untested or simply not  yet
supported by the build system.

### Via pip

If you only plan on calling MPsym from Python, then the easiest way to install
it is via `pip install pympsym`. This requires `Python >= 3.6` and `pip >= 19.3`.

### From Source

If you plan to install MPsym from source you will need the following installed
on your system:

* `CMake >= 3.6`
* `Boost >= 1.72.0`
* `Lua >= 5.2.0`
* `LuaRocks`

If you only plan on calling MPsym from Python simply run `pip install .` or
`pip install --user .`.

If you want to call MPsym from C++, you need to directly build MPsym using
CMake. Run the following from the root of this repository:

```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
make -j $(nproc)
```
Afterwards, run `make install` to install the C++ header files and shared
objects as well as the `mpsym` Lua rock to your system. If you do not want to
install the Lua rock you can pass `-DLUA_EMBED=ON` to CMake to embed the
`mpsym` Lua module into a shared object. If you do not want to use Lua
[architecture graph configuration files](#defining-architecture-files) at all,
you can instead pass `DLUA_NO_ROCK=ON` to CMake.

You can also pass `-DPYTHON_BINDINGS=ON` to CMake to additionally install the
Python bindings without separately invoking `pip`.

## Examples

The following brief examples showcase how to use the Python interface of MPsym.
For more examples, in both Python and C++, take a look at the unit tests under
`test/tests`.

### Defining Architecture Graphs

MPsym can parse architecture graph descriptions given as Lua scripts or JSON
files (the latter mostly for [serialization
purposes](#persisting-architecture-graphs)).  It is also possible to construct
architecture graphs programmatically.

The aforementioned Lua scripts must return a table describing an architecture
graph. This table is constructed using the `mpsym` Lua module.  This module
defines several convenience functions and tables which can be used to quickly
construct complex architecture graphs via the following steps:

1. Define a set of processing elements.
2. Connect them with communication channels.
3. Construct an `ArchGraph` table.
4. (Optional) Repeat steps 1-4 to construct additional `ArchGraph` tables.
5. (Optional) Compose the constructed `ArchGraph` tables using `ArchGraphCluster` and `ArchUniformSuperGraph`.

Here, `ArchGraphCluster` represents a set of architecture graphs "isolated"
from each other, such that mappings within an orbit never map the same task to
different architecture graphs within the cluster. In contrast ,
`ArchUniformSuperGraph` represents a "super graph", i.e. a set of identical
architecture graphs connected among each other.

As a simple example, let's consider again the architecture graph from our
[introductory example](#a-motivating-example).  We can construct it as follows:

```lua
local mpsym = require 'mpsym'

return mpsym.ArchGraph:create{
  directed = false,
  processors = {
    {0, 'P'},
    {1, 'P'},
    {2, 'P'},
    {3, 'P'}
  },
  channels = {
    {0, 1, 'C'},
    {1, 2, 'C'},
    {2, 3, 'C'},
    {3, 0, 'C'}
  }
}
```

Here, the integers are processor indices and `P` / `C` are arbitrary processing
element and communication channel-type labels. The exact form of these labels
is not important, but it is cruccial, that identical processing
elements/communication channels receive identical labels.

Explicitly specifying all processing elements and communication channels can be
tedious and error prone for more complex architecture graphs. While it is often
possible to more succinctly construct the respective tables using Lua language
features, the `mpsym` module also offers several convenience functions for this
purpose.  For instance,  the architecture graph from the previous example can
be more easily constructed as such:

```lua
local mpsym = require 'mpsym'

local processors = mpsym.identical_processors(4, 'P')
local channels = mpsym.grid_channels(processors, 'C')

return mpsym.ArchGraph:create{
  directed = false,
  processors = processors,
  channels = channels
}
```

We can parse a Lua configuration file in Python as follows:

```python
import pympsym

ag = pympsym.ArchGraphSystem.from_lua_file('arch_graph.lua')
```

We can also explicitly construct architecture graphs, e.g.:

```python
import pympsym

ag = pympsym.ArchGraph()

ag.add_processors(4, 'P')

for i in range(4):
    ag.add_channel(i, (i + 1) % 4, 'C')
```

### Hierarchical Architecture Graphs

MPsym can determine representatives especially efficiently when working with
certain hierarchical graph. For example, Kalray's MPPA3 Coolidge processor
consists or sixteen identical clusters of processing elements which are fully
connected internally (via shared memory) and connected in a grid fashion among
each other. We can use `ArchUniformSuperGraph` to model this (here "proto"
refers to the architecture of each cluster and "super" to the interconnections
between them):

```lua
local mpsym = require 'mpsym'

local super_graph_clusters = mpsym.identical_clusters(16, 'SoC')
local super_graph_channels = mpsym.grid_channels(super_graph_clusters, 'C')

local proto_processors = mpsym.identical_processors(16, 'P')
local proto_channels = mpsym.fully_connected_channels(proto_processors, 'shared memory')

return mpsym.ArchUniformSuperGraph:create{
  super_graph = mpsym.ArchGraph:create{
    directed = false,
    clusters = super_graph_clusters,
    channels = super_graph_channels
  },
  proto = mpsym.ArchGraph:create{
    directed = false,
    processors = proto_processors,
    channels = proto_channels
  }
}
```

This also works in Python:

```python
import pympsym

ag_super = pympsym.ArchGraph()
# ... build super graph

ag_proto = pympsym.ArchGraph()
# ... build proto graph

ag = pympsym.ArchUniformSuperGraph(ag_super, ag_proto)
```

Both `ArchGraph` and `ArchUniformSuperGraph` are subclasses of the abstract
`ArchGraphSystem` class which defines methods for [determining orbits,
representatives etc](#orbits-and-representatives). There is also
`ArchGraphCluster` which combines an arbitrary number of different and
unconnected architecture graphs.

### Persisting Architecture Graphs

Once constructed, it's possible to convert `ArchGraphSystem` objects to and
from JSON. This is especially useful because
[initializing](#initializing-architecture-graphs) an architecture graph after
its construction can require a non-trivial amount of computing time which can
be avoided on subsequent runs by loading an already initialized architecture
graph from a JSON file.

```python
>>> ag.to_json()
'{"automorphisms": [4,[1, 0],["(0, 1)(2, 3)", "(0, 3)", "(1, 2)"]]}'
>>> ag = pympsym.ArchGraphSystem.from_json(...)
```

### Initializing Architecture Graphs

Before we can perform any useful operations on an `ArchGraphSystem` object, we
need to _initialize_ it. This is a separate step that must be performed exactly
once after construction (and again if the architecture graph changes):

```python
>>> ag.initialize()
```

Note that it is not necessary to call `ArchGraphSystem.initialize` explicitly.
Operations that require it will call this method implicitly. However, since it
might take a significant amount of time to complete, it can be more sensible to
separate initialization and use of an architecture graph. We can also pass a
`timeout` argument, such that an exception will be raised if initialization does
not complete in the given number of seconds:

```python
>>> ag.initialize(timeout=2.5) # raise exception after 2.5 seconds of runtime
```

Several other `ArchGraphSystem` methods also take a `timeout` parameter that
works the same way.

### Orbits and Representatives

Given an architecture graph, we can easily determine the orbit of an arbitrary
mapping:

```python
>>> list(ag.orbit((0,1)))
[(0, 1), (0, 3), (1, 0), (1, 2), (2, 1), (2, 3), (3, 0), (3, 2)]
>>> list(ag.orbit((0,2)))
[(0, 2), (1, 3), (2, 0), (3, 1)]
```

Orbits are constructed lazily, i.e. the orbit elements are determined
incrementally while iterating through the object returned by
`ArchGraphSystem.orbit`. The lexicographically smallest mapping in each
orbit is its representative. We can also directly determine this representative
(possibly much more efficiently) as follows:

```python
>>> ag.representative((1,0))
(0, 1)
```

The `method` argument controls how the representative is determined. `iterate`
and `orbit` always produce the correct representative and which one is faster
depends on the given architecture graph and mapping. `local_search_bfs` and
`local_search_dfs` are very fast, but the returned representative is not
guaranteed to be correct (the likelihood of an incorrect result again varies
with architecture graphs and mappings):

```python
>>> ag.representative((1,0), method='auto') # default
(0, 1)
>>> ag.representative((1,0), method='iterate') # iterate through automorphism group
(0, 1)
>>> ag.representative((1,0), method='orbit') # enumerate orbit
(0, 1)
>>> ag.representative((1,0), method='local_search_bfs') # BFS local search
(0, 1)
>>> ag.representative((1,0), method='local_search_dfs') # DFS local search
(0, 1)
```

`ArchGraphSystem.representative` also takes an optional parameter of type
`Representatives` which conveniently stores all determined representatives and
causes `ArchGraphSystem.representative` to return a boolean flag and an integer
in addition to the determined representative. The boolean flag indicates
whether or not the representative has not been encountered before and the
integer uniquely identifies the orbit which the representative belongs to.

```python
>>> representatives = pympsym.Representatives
>>> ag.representative((0, 1), representatives)
((0, 1), True, 0)
>>> ag.representative((0, 2), representatives)
((0, 2), True, 1)
>>> ag.representative((0, 3), representatives)
((0, 1), False, 0)
```

This makes it possible to pass a number of mappings to
`ArchGraphSystem.representative` in sequence and to immediately decide whether
or not to e.g. skip a computationally expensive simulation step if the current
mapping is equivalent by symmetry to a previously simulated one:

```python
mappings = [...]

simulation_results = {}

for mapping in mappings:
    _, new, index = ag.representative(mapping)

    if new:
        simulation_results[index] = simulate(mapping)

    print('simulation results: {}'.format(simulation_results[index]))
```

### Automorphism Groups

We can directly retrieve the automorphism group of an `ArchGraphSystem` object:

```python
>>> pg = ag.automorphisms()
```

The returned object is of type `PermGroup` and acts like a set of `Perm` objects:

```python
>>> pg.degree()
4
>>> len(pg)
8
>>> next(iter(pg))
(0, 1)(2, 3)
>>> '(1,2)' in pg
True
>>> '(1,3)' in pg
False
```

The "degree" of a permutation group is the largest element it acts on + 1, for
an architecture graph's automorphism group this corresponds to the number of
processing elements. Each `Perm` object in turn represents a permutation of the
architecture graphs vertixes:

```python
>>> p.degree()
4
>>> p
(0, 1)(2, 3)
>>> p[0]
1
>>> p[2]
3
>>> p[4]
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
IndexError: not in domain
```

We can also directly construct permutation groups and use them in place of
architecture graphs via the `ArchGraphAutomorphisms` class. This can be useful
when the automorphism group of an architecture is known ahead of time and there
is thus no need to let MPsym determine it:

```python
>>> pg = pympsym.PermGroup(5, ['(0, 1)', '(3, 4)'])
>>> ag = pympsym.ArchGraphAutomorphisms(pg)
```

A number of convenience methods are available to construct and combine common
permutation groups:

```python
>>> pg1 = pympsym.PermGroup.symmetric(5)            # S_5
>>> pg2 = pympsym.PermGroup.cyclic(10)              # C_10
>>> pg = pympsym.PermGroup.direct_product(pg1, pg2) # S_5 x C_10
```

## Limitations

* Architecture graph initialization could be significantly sped up by employing
  more advanced BSGS construction algorithms.
* Heuristic methods could improve accuracy of local search but this would
  require significant experimentation and fine-tuning work.
* MPsym contains code for automatically decomposing hierarchical architecture
  graphs but it is as of now mostly untested and not explicitly useable from
  either the C++ or Python interface.
* MPsym also contains code for dealing with _partial symmetries_ derived from
  [[5]](#5), however, this is currently broken due to both technical and
  theoretical problems.

## Developer Notes

### Main Project Structure

The `include` and `source` directories contain the main C++ code which can be
compiled into a shared object (or static library if the `LINK_STATIC` CMake
flag is set). The C++ code has four dependencies:

* `Boost`
* `Lua`
* `nauty`
* `nlohmann/json`

The former two must already be installed on your system. The latter two are
automatically downloaded during the CMake configuration step. `nauty` is
compiled into a separate shared object (or static library), see
`nauty/CMakeLists.txt`.

The `lua` directory contains the `mpsym.lua` Lua module which can be used to
construct Lua architecture graph description files. When trying to parse these
files from MPsym, `mpsym.lua` must thus be made available via the `LUA_PATH`
environment variable, i.e. by setting:

```bash
export LUA_PATH=$LUA_PATH;$(readlink -f mpsym/lua/?.lua)
```

Alternatively, you can set the `LUA_EMBED` CMake flag to embed `mpsym.lua` into
the MPsym shared object/static library (which is arguably a weird thing to do
but quite handy in practice).

The `test/tests` directory contains C++ unit tests. these are built by CMake
if `-DCMAKE_BUILD_TYPE=Debug` is specified. The tests use the `Googletest`
framework which is automatically downloaded during the CMake configuration step.

The `python` directory contains Python binding code and tests. It has the
following structure:

```
python/
├── setup.py
├── mpsym
│   ├── __init__.py
│   └── _mpsym_tests.py
└── source
    ├── CMakeLists.txt
    └── _mpsym.cpp
```

`python/mpsym` is the binding module directory. `python/_mpsym.cpp` contains
[pybind11](https://github.com/pybind/pybind11) wrapper code for MPsym's public
C++ interface. When the `PYTHON_BINDINGS` CMake flag is set, a corresponding
shared object is created under `python/mpsym/_mpsym.*.so`. `_mpsym_tests.py`
contains a number of unit tests. The `mpsym` binding module loads both the
module created via pybind11 and `_mpsym_tests.py` in its `__init__.py`. The
latter can be run by invoking `mpsym.test(verbosity=...)`, a return value of
`0` indicates success.

If you don't care about the C++ interface, you can directly install the binding
module to your system via `python setup.py install --user`.

### Profiling

Besides `Release` and `Debug`, a third build mode, `Profile`, is also
supported.  In this mode, the programs under `profile/source` are compiled.
They can be used to profile the runtime of the Schreier-Sims algorithm as
implemented by MPsym, as well as the various canonical representative
algorithms. These programs implement `--help` flags that should more or less
explain how to use them. Some related example architecture graphs and scripts
can be found [here](https://github.com/Time0o/mpsym_experiments).

### Deploying

Running `deploy.sh` will create test coverage data and Doxygen documentation
and will upload these to Codecov and GitHub pages respectively. Of course you
shouldn't be able to do this unless you're me :o).

### Continuous Integration

Previously, MPsym used Travis for CI. Since Travis is unfortunately no longer
free for FOSS projects, GitHub Actions are now used instead. Take a look at
`.github/workflows/workflow.yml` for details. On every commit and pull request
to the master branch, MPsym is built and tested inside a Ubuntu/macOS image
using recent versions of gcc/clang.

Previously, Travis also took care of deployment. To save on GitHub Action
credits, coverage and documentation must now be deployed manually. PyPi wheels
are still built automatically (for Linux and macOS), but now in a [separate
repository](https://github.com/Time0o/mpsym_wheels) that uses
[multibuild](https://github.com/matthew-brett/multibuild) (currently a work in
progress).

## References

<a id="1">[1]</a>
Holt, D. F. (2005).
_Handbook of Computational Group Theory_.
CRC Press

<a id="2">[2]</a>
Groups, T. G. (2020).
_GAP - Groups, Algorithms, and Programming, Version 4.11.0_.

<a id="3">[3]</a>
McKay, B. D. and Piperno, A. (2014).
_Practical graph isomorphism, ii._
Journal of Symbolic Computation, 60:94–112.

<a id="4">[4]</a>
Donaldson, A. F. and Miller, A. (2009).
_On the constructive orbit problem_.
Annals of Mathematics and Artificial Intelligence 57:1-35.

<a id="5">[5]</a>
East, J., Egri-Nagy, A., Mitchell, J., and Péresse, Y. (2019).
_Computing finite semigroups_.
Journal of Symbolic Computation, 92:110–155.
