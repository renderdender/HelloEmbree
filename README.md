# HelloEmbree
This is a template used for our Advanced Computer Graphics course at [MSU Graphics Lab](http://graphics.cs.msu.ru/en).

Basically, it's a very simple ray-tracer that writes the color of the first hit to the image and provides some useful utilities like image i/o, .obj model parser, e.t.c.
## Installing Dependencies - Intel Embree

Go to [the official Intel Embree Page](https://github.com/embree/embree). Download the [tar.gz archive](https://github.com/embree/embree#linux-targz-files) and extract it. Open CmakeLists.txt and replace the line:
> set(EMBREE_INSTALL_DIR      "/home/denispavlov/Software/embree-3.2.4.x86_64.linux")

with the path of your extracted Intel Embree library.

## Building the Project

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

![Template Image](https://github.com/pavlovdenis/HelloEmbree/pics/template.png)
