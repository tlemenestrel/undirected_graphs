# Undirected Graphs

## Introduction

This repository contains a collection of **C++ softwares** based on the implementation of an _Undirected Graph_ Class. This allows for the usage of **ODE** and **PDE** solvers, the implementation of a _Breadth-First Search (BFS) algorithm_, a solver of systems of _sparse matrices_ and more, all using this class.
    
## Table of Contents

- [Graph Class](#Graph-Class)
- [Breadth-First Search (BFS) Algorithm](#BFS)
- [Mass Spring](#Mass-Spring)
- [Conjugate Gradient and Sparse Matrices](#Conjugate-Gradient-and-Sparse-Matrices)
- [Parallel Computing using OpenMP](#Parallel-Computing)

## Graph Class

<table>
<tr>
<td>
  
The Graph class is implemented using two sub-classes for Nodes and Edges. A Graph is made of two Hash Maps (called unorderded maps in the C++ STL), which contain the different Nodes and Edges of the Graph. The Nodes are implemented using a **Design Pattern** called a _Proxy Pattern_. Instead of storing all the information of a Node in the Node class, we define a Proxy Class which points to the underlying data of the Node. This allows for much lighter data structures and better encapsulation.

All three classes (Node, Edge and Graph) are thoroughly tested using GoogleTest, a testing framework for C++ code by Google.

<p align="center">
<img src="https://github.com/tlemenestrel/undirected_graphs/blob/main/data/53919192-8211-4DC7-8EF5-B69487EA4384.png" width="700">
</p>

</td>
</tr>
</table>

## Breadth-First Search (BFS) Algorithm

<table>
<tr>
<td>

After implementing the Graph Class, a viewer is implemented using SFML combined with OpenGL for high-performance graphics. It is able to read sets of Nodes and Egdes files that contain their respective 3D positions.

Thanks to appropriate usage of C++ STL data structures, the code is able to read through 10,000 Nodes and 300,000 Egdes within less than 1 second and automatically output the result. After this, random 3D coordinates are chosen. Starting from the closest Node to those coordinates, a BFS algorithm is implemented to parse through all Nodes of the Graph and compute the shortest path in terms of the number of Edges to go through to get to the original Node. Using this, each Node is coloured according to how close or far it is from the original Node.
<p align="center">
<img src="https://github.com/tlemenestrel/undirected_graphs/blob/main/data/A6760F25-AD77-4D39-BEE5-18B1B57BD93B.jpeg" width="700">
</p>

</td>
</tr>
</table>
