

## Plotting extra data on PSSECT plots

`pssect` is a program that makes plots of a calculation done by `vertex`.
It has a control file called `perplex_plot_option.dat` that selects
specific options for a plot.

If the plot option file contains the option, `plot_extra_data` and has the
value `T` (true), then `pssect` reads a file that is used to annotate the plot.
The annotation file name is read by `pssect` from the standard input and
processed for further additions to the basic plot.

`pssect` has a built-in color table whose colors are selected by an integer
number in many of the following commands.  There is also a set of line styles
that may be used for line drawing.

| Code | Color | Code | Color      || Code | Line            | Code| Line |
| :--: | ----- | :--: | -----      |-| :--: | :--:            | :--:| :--: |
| 0    | black |  7   | brown      ||  0   | none            |  7  | very short dash|
| 1    | white |  8   | orange     ||  1   | solid           |  8  | sparse dots |
| 2    | red   |  9   | dark blue  ||  2   | short dash      |  9  | heavy dots |
| 3    | green |  10  | dark red   ||  3   | uneven dash     | 10  | sparse dash |
| 4    | blue  |  11  | dark green ||  4   | very sparse dots|   | |
| 5    | purple|  12  | dark yellow||  5   | long dash | | |
| 6    | yellow|      |            ||  6   | dashed | | |

The types of annotations are:
* comments (not processed by `pssect` but for human readability/understanding)
* lines (with or without symbols at each point)
* points (with or without associated error bars)
* text line or block of text, the individual lines of which are justified

### Comments

Any line with `*` in the first column is a comment, and any
text on any line following a `#` or a `|` is ignored.

### Lines

A line is drawn on the resulting plot.  Lines are introduced by

`> L` *options*\
...\
`>`
   
Zero or more (x,y) pairs appear in the ... part, in the units of the plot
(e.g. P in bars, T in K, composition in mole fraction).
The points trace out a line drawn on the plot.

*options* indicates a group of keyword=value pairs that govern the drawing of
the line.  Presently recognized keywords are:
* col=*number* - set the line color (see above)
* lty=*number* - set the line style (see above)
* lwd=*number* - set the line width (1 is "normal")
* sym=*number* - draw symbol *number* at each line point (see below)
* ssiz=*number* - line symbol size (1 is "normal")
* sfil=*number* - line symbol fill color (see color table)

### Points

Points are described with a location, symbol type, size, and color
fill.  Points may have error bars plotted along each coordinate.
The general form of a point is:
```
x y symb size color
x y ex ey symb size color
```
The (`x`,`y`) values give the coordinates of the point.  (`ex`,`ey`) give the
length of the error bars for the symbol.
`size` is the size of the symbol (real number) - 1 is "normal."
`color` is a color for the symbol outline (see color table above) or fill.

`symb` is a symbol code from this table:
n| symbol       |n | symbol                  | n| symbol
-:|  :--:       |-:|  :--:                   |-:| :--:
0| square       |10|O + overlaid             |20|outlined small circle
1| circle       |11|up/down triangle overlaid|21|outlined small & large circle
2| up triangle  |12|+ square overlaid        |22|outlined small & large square
3| plus (+)     |13|O x overlaid             |23|outlined small & large diamond
4| times (x)    |14|square up triangle overlaid|24|outlined up triangle
5| diamond      |15|solid square             |25|outlined down triangle
6| down triangle|16|solid down triangle      |26|outlined circle
7| square       |17|solid up triangle        |28|outlined square
8| + x overlaid |18|solid diamond            |29|outlined diamond
9|+ diamond overlaid|19|solid circle         |  |

0-14 are outlined symbols; 15-19 are a solid color; and the rest are outlined
in black and then filled with a color.

### Text

Text is drawn on the resulting plot.  Text lines (or blocks)  are introduced by

`> T x y` *options*\
...\
`>`
   
`x y` are the coordinates of the text, in units of the plot (the same as for
points and lines).
Zero or more lines of text make up the ... part.
The text may be drawn at an angle, as well as being placed relative to a point
(above, below, beside, etc; see *adj*, below).
If there is more than one line, the lines are treated as a block of text and
may be justified (left- or right-aligned or centered).

*options* indicates a group of keyword=value pairs that govern the display of
the text.  Presently recognized keywords are:
* col=*number* or col=*number*,*number*- set the text color
(foreground and, optionally, background; see above)
* siz=*number* - set the text size (1 is "normal")
* sym=*number* - draw symbol *number* (see above) at `x y` along with the text
* ssiz=*number* - symbol size (1 is "normal")
* sfil=*number* - symbol fill (see color table above)
* ang=*number* - angle at which to draw the text
* adj=*number*,*number* - position of `x y` relative to the text, in units of
the text width and height.  For example, 0,0 is bottom left, 1,1 is top right,
-0.05,0.5 is left-justified and centered on the text height,
0.5,1.1 is centered text below the point
* just=L, just=R, just=C - justification of a block of text (multiple lines):
left, right or centered
* box=*number*,*number* - offset the box by (*number*,*number*) from `x y` but
connect the text to the point with a pointer; the box is outlined with the
foreground color and filled with the background color; any adj= and ang= are ignored if specified

Examples:
```
> T 50e4 4500 sym=26 ssiz=1 sfil=1 adj=-0.15,0.5 ang=45 col=2 siz=0.5
OPX
>
```
This plots a symbol (26) at (50e4,4500) of normal size (1) and color filled (1).
The text (OPX) is centered vertically and offset to the right of the point, and
written in color (2) at an angle of 45 degrees, in a font half of the normal size (0.5).

```
> T 0.5 1500 sym=1 box=0.3,-100 just=C col=2,1
Reversal
bracket
>
```
This plots a symbol (1) at (0.5,1500) and puts the text centered (C)
line-by-line in a box offset by (0.3,-100) connected to the point with an arrow.
The text is written in the foreground color (2) against the
background color (1) that fills the box.

```
> T 132e4 2700 adj=0.5,-0.1 
(liquid)
>
```
This puts a label (liquid) centered above the point (132e4,2700).
