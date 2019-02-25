# WM_S288C: whole-cell model for S. cerevisiae
This package provides the S. cerevisiae whole-cell knowledge base and simulation, the actively developed repo is in Matlab at https://github.com/fmmelab/WM_S288C.

Required software

(1) MATLAB R2009b (7.9) or newer is required for the whole-cell model software.

(2) An XAMPP (5.6.40) installation package, which brings together the Apache web server with PHP, Perl, and MariaDB, can be downloaded from https://www.apachefriends.org/download.html.

Installation guide

(1) Building of local mysql database:

After install the XAMPP software, s288c.sql, which contains the whole data used for construction of S. cerevisiae whole-cell model, should be imported.

(2) Configing the whole-cell toolbox:

First， the whole-cell toolbox should be added to the path of MATLAB,

Then, in MATLAB, change the path to the root directory of this toolbox,

Last, run install.m and follow the instructions on screen.


Demo


Instructions for use
