# Point in polygon problem
Implementation of the *ray-crossing algorithm* to decide whether a point is inside or outside a polygon.

## 1 - Objective

The aim of this project is to give an implementation of the ray-crossing algorithm to decide whether a point is inside or outside a polygon and is purely pedagogical in purpose. The algorithm is implemented in both Euclidean dn spherical spaces.

## 2 - Repo organisation

**`point_in_polygon/`: The point-in-polygon modules**
It contains the point-in-polygon modules. See also Module architecture below.

**`notebooks/:` Notebooks demonstrating the modules**
- `PointInPolygon.ipynb`: Notebook discussing the basics of the the ray-crossing algorithm and showing the main fetures of the modules

## 3 - Module architecture

Description of the module `molecular_dynamics` architecture.

- `point_in_polygon/__init__.py`
  - Initialises the module.
  - Imports the `PolygonE` class that handles polygons on Euclidean space.
  - Imports the `PolygonS` class that handles polygons on spherical space.

- `molecular_dynamics/polygon_euclidean.py`: defines the `PolygonE` class with the methods
  -  `creates`,
  -  `internal_point`.
- `molecular_dynamics/polygon_spherical.py`: defines the `PolygonS` class with the same methods as above.

## 4 - Features

- The `PolygonE` class:
  - generates a polygon in Euclidean space;
  - the polygons can be regular, irregular, concave or convex;
  - the irregular polygons are generated randomly;
  - the vertices can be ordered in clockwise or anticlockwise directions;
  - the generated polygons are mainly to test the ray-corssing algorithm for general enough polygon shapes. 
  - It has the following methods:
    - `create`: generates the vertices of the polygon from a chosen anchor point;
    - `internal_point`, check whether a given point is inside or outside the polygon using the ray-crossing algorithm.

- The `PolygonS` class:
  - generates a polygon on a 2-sphere;
  - the polygons can be regular, irregular, concave or convex;
  - the irregular polygons are generated randomly;
  - the vertices can be ordered in clockwise or anticlockwise directions;
  - the generated polygons are mainly to test the ray-corssing algorithm for general enough polygon shapes. 
  - It has the same methods as the `PolygonE` class.

## 5 - Results

TO DO

## 6 - Bibliography

TO DO
