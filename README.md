# RouteEditor
ISOBUS Guidance Path editor

* Author: Juha Backman
* Contact: juha.backman@luke.fi
* Org: [Natural Resources Institute Finland (Luke)](https://www.luke.fi/en)

## Overview
ISOBUS Task file is commonly used in agriculture for precision farming. Most of the softwares that creates Task files is made to generate prescription maps for agricultural implements. However, the Task file can contain also guidance paths and other usefull information for robots to be used. See more: [Building a Robot Tractor Using Commercial Components and Widely Used Standards (Backman et.al, 2022)](https://doi.org/10.1016/j.ifacol.2022.11.106)

This program is a tool for creating and editing guidance paths that can be used either for precision farming, or for example robot paths.

This program also include "headland"-function that generates a path from start to goal state that respects constraints on the 0th, 1st and 2nd derivative of curvature and velocity. The algorithm is originally published in
[Smooth turning path generation for agricultural vehicles in headlands (backman et.al, 2015)](https://doi.org/10.1016/j.biosystemseng.2015.08.005). In this program it is used to generate a headland turning or connection between other driving paths.

## Dependencies

The program is made using Matlab and [Matlab App Designer](https://se.mathworks.com/help/matlab/ref/appdesigner.html). 

headland.m is a standalone function that can be used in Matlab without Route Editor. In addition, Matlab's m-scripts can be converted to Python quite easily. If you need only the code that generate a path between two arbitrary positions, you can also use the Python implementation by Rayne Milner:  https://github.com/rmm286/backman_smooth_turning_path_generation

## Usage

The usage of the program is explained with two different scenarios. In addition to the two scenarios described below, the guidance path can be created and edited point by point.

### Creating a path using obstacles

The first scenario is the creation of a path for a robot that moves in a strawberry tunnel. The locations of the strawberry table legs are used to create a path at a constant distance from the legs. Between the rows, a turning is created using parameters obtained from the robot (max turning angle, lock-to-lock steering time etc.) and the headland.m function.

![Screenshot](/images/Strawberry%20tunnel%201.png)

First, the locations of the strawberry table legs are loaded. It is done by loading a existing Task file that contains the legs locations described as obstacles.

The file is loaded using `File->Load`

In the popup window, select only `Obstacles`

![Screenshot](/images/Strawberry%20tunnel%202.png)

The editor shows the content that was loaded from the task file. The map can be zoomed using the tools found in the upper right corner of the map.

![Screenshot](/images/Strawberry%20tunnel%203.png)

The path is create by selecting `Copy` command using the button on the right. Then the first obstacle is selected using left mouse button. The selected obstacle is shown as red cross.

![Screenshot](/images/Strawberry%20tunnel%204.png)

The copied locations can be selected one by one, or it is possible to select the entire row by selecting any other obstacle in the row with the right mouse button.

![Screenshot](/images/Strawberry%20tunnel%205.png)

When selection of the copied locations is done, the `Copy` function is ended by pushing a `Ready` button. After that, in the popup window fill other information.  

![Screenshot](/images/Strawberry%20tunnel%206.png)

If there are existing paths in the list, the editor asks if the new path will be connected to the previous path automatically.

![Screenshot](/images/Strawberry%20tunnel%207.png)

If the path is connected automatically, fill the information of the tractor or robot. Also turning type can be selected manually, or the fastest can be chosen automatically. The method generates a path that respects the given parameters.

![Screenshot](/images/Strawberry%20tunnel%208.png)

The path now consists of three path segments. The last position of the path is marked with blue cross.

![Screenshot](/images/Strawberry%20tunnel%209.png)

By selecting `Modify path` command, the existing paths can be edited. The paths to be modified can be selected by selecting those onthe map with the mouse, or after pressing `Ready` button, from the list.

![Screenshot](/images/Strawberry%20tunnel%2010.png)

The available modify functions are:

* Merge paths: the selected paths will be merged together
* Insert zero to direction change: If the driving speed goes from negative to positive or vice versa, a position with zero speed is added.
* Cut paths when speed is zero: the path is divided to two different paths if the driving speed is zero
* Normalize point distances: the waypoints are set at a constant distance from each other 

One ore more functions can be selected at the same time. The functions are executed in the top-down order.

In addition to these functions, the path points can be deleted, inserted or moved one by one using the command buttons found from right.


![Screenshot](/images/Strawberry%20tunnel%2011.png)

If the `Normalize point distances` was selected, the desired distance between waypoints is asked. The current waypoints are marked with the red cross in the map.

![Screenshot](/images/Strawberry%20tunnel%2012.png)

After the execution of the functions, the result looks like this.

![Screenshot](/images/Strawberry%20tunnel%200.png)

The final map for the whole strawberry tunnel. The area where work is done (work state is 1), is drawn green.

### Creating a path using field boundaries

The second scenario is the creation of a path for a tractor driving in field.

![Screenshot](/images/Field%200.png)

First, the field boundaries is selected like in the first scenario.

![Screenshot](/images/Field%202.png)

The first path is created by selecting `Copy` command and one point from the field boundary with left mouse button and then continuing selection with right mouse button. The straight lines are selected automatically. When desired partion of the field boundary is selected, create the path by pressing `Ready` button as in first scenario.

![Screenshot](/images/Field%204.png)

After the first path, select again `Copy` command and the previous path with right mouse button. In the popup window, select `Direction: -1`. 

![Screenshot](/images/Field%203.png)

The path is created at a constant distance from the previous path. There might be extra points in the the corners that can be deleted manually. Also path starting or ending points can be edited manually. 

![Screenshot](/images/Field%201.png)

The`Copy` function can be continued until the desired area of the field is covered.
