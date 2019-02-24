The SDT format
**************
The SDT format (for *SOFT Descriptive Text* or *Self-Descriptive Text* or
*Semi-Descriptive Text*) was developed in order to import magnetic
field data from systems with no HDF5 or MATLAB support. It is a very simple
text-based format without fancy features and with little safety. It is
recommended that users stick to HDF5 or MATLAB files whenever possible.

Basic structure
---------------
Just as MATLAB files (and in a sense, HD5 files), SDT files contains a set of
variables. Each variable consists of a *header* and a *body*. The header is
always one line and specifies the name and dimensions of the variable. The body
(which comes on the very next line) is at least one line and contains the data
of the variable, in ASCII format. Variables are separated by empty lines.

The header always consists of two integers and a string, all separated by
spaces. The first integer specifies the number of rows in the data. The second
integer specifies the number of columns in the data (or number of characters, in
the case of strings). The name is an ASCII string of arbitrary length (but
without any whitespace characters).

Note that strings and numeric variables are, technically, encoded differently.
Data should always be readable by a human in a text-editor, meaning that numeric
values are converted to their ASCII equivalent, while strings are stored
directly without converting every single character to an ASCII digit.
There is no indication in the header about which type of data a variable
contains and so it is up to the user to read variables using the correct
decoder.

Matrices are stored by converting each element to its ASCII equivalent. Elements
in the same row are separated by single spaces, while rows are separated by
(Unix) newlines (i.e. just one ``\n`` or ``0xA`` character).

Example SDT file
----------------
An example SDT file containing three variables is shown below::

   1 2 maxis
   1.688510 0.048496

   3 3 someMatrix
   1.1 2.2 3.3
   4.4 5.5 6.6
   7.7 8.8 9.9

   1 29 someString
   This is an SDT example string

