# Undirected Graphs

## Introduction

This repository contains a collection of **C++ softwares** based on the implementation of an _Undirected Graph_ Class. This allows for the usage of **ODE** and **PDE** solvers, the implementation of a _Breadth-First Search (BFS) algorithm_, a solver of systems of _sparse matrices_ and more, all using this class.
    
## Table of Contents

- [Graph Class](#Graph-Class)
- [Breadth-First Search (BFS) Algorithm](#Breadth-First-Search-(BFS)-Algorithm)
- [Mass Spring](#Mass-Spring)
- [Conjugate Gradient and Sparse Matrices](#Conjugate-Gradient-and-Sparse-Matrices)
- [Parallel Computing using OpenMP](#Parallel-Computing-using-OpenMP)

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

Thanks to appropriate usage of C++ STL data structures, the code is able to read through 10,000 Nodes and 300,000 Egdes within less than 1 second and automatically output the result. After this, a set of random 3D coordinates are chosen. Starting from the closest Node to those coordinates, a BFS algorithm is implemented to parse through all Nodes of the Graph and compute the shortest path in terms of the number of Edges to go through to get to the original Node. Using this, each Node is coloured according to how close or far it is from the original Node.
<p align="center">
<img src="https://github.com/tlemenestrel/undirected_graphs/blob/main/data/A6760F25-AD77-4D39-BEE5-18B1B57BD93B.jpeg" width="700">
</p>

</td>
</tr>
</table>

## Mass Spring

<table>
<tr>
<td>

We use the Graph class to implement a solver for Spring Mass systems. Each node has zero as initial velocity and a mass equal to 1 / N, N being the number of nodes in the graph. We have a spring constant set to 100 and a rest-length L set uniformly across all Edges of the Graph. The mass spring simulation works by solving a differential equation using an iterative approach.

To solve the numerical system, we update the position based on the velocity, compute the force using the new position and finally update the velocity.

We use an OOP framework to represent forces by building a parent Force class and having 3 subclasses inheriting from it for gravity, damping and mass spring.

<p align="center">
<img src="https://github.com/tlemenestrel/undirected_graphs/blob/main/data/9E269D24-D723-4089-BF5B-766B7B2C3FFE.jpeg" width="700">
</p>

</td>
</tr>
</table>

## Conjugate Gradient and Sparse Matrices

<table>
<tr>
<td>

We implement the Conjugate Gradient (CG) algorithm for solving systems of equations with symmetric matrices and represent those using the Graph class. The CG solver is implemented to work efficiently for sparse matrix-vector products. Then, we implement the boundary conditions and solve the system using the Conjugate Gradient Solver. 

Compared to the mass spring approach where we approximated derivatives in time, here we approximate them in space. We can approximate the Laplacian operator by a discrete Laplacian matrix, which we use to solve the system Ax = b.

Finally, we visualize the convergence of the solution by plotting the x and y coordinates of each node and the scalar solution to the Poisson problem of each node using the SFML viewer.

<p align="center">
<img src="https://github.com/tlemenestrel/undirected_graphs/blob/main/data/CB7867AE-F7CE-4560-88F0-8346EF1DBFC8.png" width="700">
</p>

</td>
</tr>
</table>

# Parallel Computing using OpenMP



<table>
<tr>
<td>

We re-implement the Mass Spring algorithm but in a parallel manner by using OpenMP and Thrust. This allows for much faster processing as each velocity and position of the Nodes can be computed at the same time.  

This allows to speed up the mass spring computation by a factor of 40 on a large grid (1M Nodes and 4M Edges).

<p align="center">
<img src="https://github.com/tlemenestrel/undirected_graphs/blob/main/data/BFDB2FA3-24D7-48EC-AF7C-C893DD7925AE.png" width="700">
</p>

</td>
</tr>
</table>
