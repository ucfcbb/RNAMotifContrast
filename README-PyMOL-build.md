### Steps to build PyMOL from open-source code (Python 3.6+ required)

#### 5.1 Installing required libraries

The following commands can be used to install all the required libraries (must be run as root):

```
# Debian/Ubuntu/Mint
apt-get install git build-essential python3-dev libglew-dev \
  libpng-dev libfreetype6-dev libxml2-dev libmsgpack-dev \
  python3-pyqt5.qtopengl libglm-dev libnetcdf-dev freeglut3-dev

# CentOS
yum install gcc gcc-c++ kernel-devel python-devel tkinter python-pmw glew-devel \
  freeglut-devel libpng-devel freetype-devel libxml2-devel glm-devel \
  msgpack-devel netcdf-devel

# Fedora
dnf install gcc gcc-c++ kernel-devel python3-devel glew-devel PyQt5 msgpack-devel \
  freeglut-devel libpng-devel freetype-devel libxml2-devel glm-devel

# Anaconda (Both Linux and MacOS)
conda install -c menpo glew
conda install -c conda-forge glm
conda install -c anaconda netcdf4
```

In MacOS, it is also possible to install the dependecies using different package managers (e.g., Homebrew, MacPorts, Fink, etc.).

#### 5.2 Collect PyMOL open-source code and other required repositories from Git

```
git clone https://github.com/schrodinger/pymol-open-source.git
git clone https://github.com/rcsb/mmtf-cpp.git
mv mmtf-cpp/include/mmtf* pymol-open-source/include/
```

#### 5.3 Build source

```
cd pymol-open-source
prefix=$HOME/pymol-py3
python3 setup.py build install --home=$prefix --glut
```

Sources of these instructions to build PyMOL from source:
* https://pymolwiki.org/index.php/Linux_Install
* https://pymolwiki.org/index.php/MAC_Install