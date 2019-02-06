# CaLib

CaLib calibration database system

### Installation

#### Dependencies
* ROOT 5.34 / ROOT 6.10 or newer (with MySQL or/and SQLite support)
* ncurses
* CMake 2.8

#### Installation
```
cd /some/directory
git clone https://github.com/werthm/CaLib.git
cd CaLib
mkdir build
cd build
cmake ..
make -j
export CALIB="/some/directory/CaLib"
export LD_LIBRARY_PATH="$CALIB/build/lib:$LD_LIBRARY_PATH"
export PATH="$CALIB/build/bin:$PATH"
```
It is recommended to set all environment variables in your shell configuration file.

#### Upgrade from 0.3.x to 0.4.x
* The database has to be updated to version 5 using

```
root -b $CALIB/macros/Upgrade_5.C
```

#### Upgrade from 0.2.x to 0.3.x
* The database has to be updated to version 4 using

```
root -b $CALIB/macros/Upgrade_4.C
```

* Exports to ROOT files created with CaLib < 0.3.0 cannot be imported by Calib > 0.3.0!

#### Upgrade from 0.1.11 to 0.2.x
* The database has to be updated to version 3 using

```
root -b $CALIB/macros/Upgrade_3.C
```

* Exports to ROOT files created with CaLib < 0.2.0 cannot be imported by Calib > 0.2.0!

### Configuration

All the configuration is done in config/config.cfg.  
config/example.cfg contains comments and basic settings that should
help to understand the configuration. Rename this example file to config.cfg
and modify it according to your setup.

### Documentation

The ROOT html documentation can be built by setting the cmake variable -DDOCS=ON.
Located in the directory htmldoc it gives an overview of the CaLib library and 
all of its classes.
Further information and examples can be found in the macros directory.

### Changelog

#### 0.4.0
February 6, 2019
* changed number of maximum tagger channels to 408 (requires DB upgrade to version 5)
* small improvements and bugfixes

#### 0.3.0
March 29, 2018
* added SQLite support
* use CMake building
* support for ROOT 6
* support for Pizza detector
* extended calib_manager
* generic peak calibration module
* improved support for bad scaler reads
* graphics/fitting improvements (marker line, refitting, ignoring elements, convergence factors, etc.)
* added data type for the beam polarization
* bugfixes

#### 0.2.0
January 7, 2014
* added support for bad scaler reads
* added calibration cloning
* added run range setting for calibration sets
* removed livetime data type

#### 0.1.11
August 16, 2013
* added support for Mk2 format and xz-compressed files
* bugfixes

#### 0.1.10
April 26, 2012
* several bugfixes

#### 0.1.9
November 14, 2011
* added calibration cloning function
* bugfixes

#### 0.1.8
October 10, 2011
* added TAPS CFD and veto LED to AcquRoot config file writer
* bugfixes

#### 0.1.7
September 11, 2011
* updated PID energy calibration

#### 0.1.6
August 31, 2011
* added TAPS CFD and veto LED calibration

#### 0.1.5
August 24, 2011
* initial public release

