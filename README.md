# HelloEmbree

<p align="center">
  <img width="600" height="600" src="https://github.com/pavlovdenis/HelloEmbree/blob/master/pics/pt.png">
</p>

This project is a template used for our Advanced Computer Graphics course at [MSU Graphics Lab](http://graphics.cs.msu.ru/en), focused mostly on path tracing and more advanced global illumination algorithms.


Basically, it's a very simple ray-tracer that writes the color of the first hit to the image and provides some useful utilities like image i/o, .obj model parser, e.t.c. The template contains blank functions and TODO instructions for students to follow.

## Installing Dependencies - Intel Embree

Go to [the official Intel Embree Page](https://github.com/embree/embree). Download the [tar.gz archive](https://github.com/embree/embree#linux-targz-files) and extract it. Open CmakeLists.txt and replace the line:

```
set(EMBREE_INSTALL_DIR "/home/denispavlov/Software/embree-3.2.4.x86_64.linux")
```
with the path of your extracted Intel Embree library.

## Building the Project

First, install `CMake` with: `sudo apt-get install cmake`.

Go to the root of the project:

```
mkdir build
cd build
cmake ..
make
./HelloEmbree
```

Or you could just open the project with CLion.

If everything's done correctly, the program should output something like this:

<p align="center">
  <img width="600" height="600" src="https://github.com/pavlovdenis/HelloEmbree/blob/master/pics/template.png">
</p>

If you've got any troubles with the project, feel free to email me: `denis.pavlov at graphics.cs.msu.ru`
