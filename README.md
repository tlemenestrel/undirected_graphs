# Undirected Graphs

## Introduction

This repository contains a collection of **C++** softwares based on the implementation of an undirected graph Class. This allows for the usage of ODE and PDE solver, the implementation of a Breadth-First Search (BFS) algorithm, a solver of systems of sparse matrices and more, all using this class.
    
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

**Truss Linear Equations** is a Python software to define a Truss class for loading and analyzing a 2D truss using the method of joints.

<p align="center">
<img src="https://github.com/tlemenestrel/swe_scientific_projects/blob/master/truss_linear_equations/truss.png" width="700">
</p>

</td>
</tr>
</table>

### Terminal commands and outputs

1. First, cd into the truss_linear_equations folder:

    ```
    $ cd truss_linear_equations
    ```

2. Then, run the following command to process the data of the beams and truss files:

    ```
    $ python main truss1/joints.dat truss1/beams.dat
    ```
    
3. You should get the following output:
   ```
   Beam       Force
   ----------------
      1       0.000
      2      -1.000
      3       0.000
      4      -1.000
      5       0.000
      6       0.000
      7      -0.000
      8      -1.414
   ```
   

## Maze solver

<table>
<tr>
<td>

Maze Solver is a c++ software to solve a maze, given as an input txt file as a 2D array. It uses the right-hand wall following algorithm. Based on the current location and direction, the algorithm changes direction and/or moves to a new location to find the way out.

The code checks if the command line arguments are provided and reads them, reads the input file, first by reading the size of the array and then the maze. Then, it populates the array with 0s everywhere and 1s where there is no wall. After, it finds the entry of the maze by looking at the first row and runs a while loop that while iterate until the exit is found. At each step, the algorithm will attempt to go in different directions and store the positions at each step. Finally, the positions are given in the output file and make up the solution to the maze.

Afterwards, if all the arguments are provided properly, it reads the maze file and stores the dimensions of the maze. Then, it buils a CSR matrix to store the maze (useful because the maze is mostly made of 0s). After, it reads the solution file, checks that the maze was properly entered on the first row, that each position change is valid and that the code reaches the exit of the maze on the last row. Finally, it prints feedback about whether the solution is valid or invalid.

<p align="center">
<img src="https://github.com/tlemenestrel/swe_scientific_projects/blob/master/maze_solver/maze.jpg" width="500">
</p>

</td>
</tr>
</table>

## Image Processing

<table>
<tr>
<td>

Image processing is an important category of computations utilized in many science and engineering fields. Some important computations in image processing include focus detection and various filtering operations to smooth or sharpen an image, detect edges, etc.

The goal of this project is to developp a C++ image class that can read and write JPEG files, and has methods to compute the sharpness of the image and smooth the image using a box blur kernel of a specified size.

<p align="center">
<img src="https://github.com/tlemenestrel/swe_scientific_projects/blob/master/image_processing/boxblur.png" width="700">
</p>

</td>
</tr>
</table>

### Code description

The code first implements an image class with a constructor that accepts a string
containing the name of the JPEG file to be read. The constructor reads the file 
and stores the image data as a data attribute. It calls a function called
ReadGrayscaleJPEG(), which is pre-defined.The class also has a Save() method that writes the current version of the image 
to a JPEG file. It takes a string containing the name of the JPEG file to be 
written. However, if the string is empty the original filename is used instead.

Afterwards, the class has a Convolution() method. It has specific requirements
for both the kernel and the input and output. First, the input and output should
 be of the same size. Second, the method only supports square kernels of odd 
 sizes and at least size 3. Then, it has a BoxBlur() method that takes an argument specifying the kernel 
size and calls the Convolution() methd to smooth the image. Also, it uses a Sharpness() method that returns the sharpness of the image. 

Finally, the main file loads the image stanford.jpg and computes and outputs 
the sharpness for kernel sizes of 3, 7, . . . , 23, 27 and then compute and 
output the sharpness of the resulting images.
