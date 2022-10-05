Forked from http://mgl.scripps.edu/people/goodsell/research/bend/index.html on 2008.10.25:

<pre>
<h2>BEND and BEND_TRI:
<i>Programs for the Prediction of Bending and Curvature for Arbitrary Sequences</i></h2>

Version 2.0, January 1995

These simple programs calculate the magnitude of the local bending
and macroscopic curvature at each point along an arbitrary B-DNA sequence,
using any desired bending model that specifies values of twist, roll,
and tilt as a function of sequence.

The programs are described in full in:
Goodsell, D.S. &amp; Dickerson, R.E. (1994) "Bending and Curvature Calculations
in B-DNA" Nucl. Acids. Res 22, 5497-5503.

This page contains FORTRAN source for the program, and information on how
to get started immediately. Full descriptions of the input files are included
in the headers to the FORTRAN source code.

Mensur Dlakic at the University of Nevada, Reno, has developed versions
of BEND for IBM PCs and Macintoshes. Please contact him directly at
mensi@scs.unr.edu for information on availability.

---FILES INCLUDED IN THIS RELEASE-----------------------------------

Six files are included in this release:

   index.html - this file
   <a href="BEND.FOR">BEND.FOR</a> - source code for program using dimer models
   <a href="BEND_TRI.FOR">BEND_TRI.FOR </a>- source code for program using trimer models
   <a href="SEQUENCE.DAT">SEQUENCE.DAT</a> - sequence file containing an intrinsically bent
                  base sequence from kinetoplast DNA
   <a href="ANGLES.DAT">ANGLES.DAT</a> - sample dimer angle data
   <a href="ANGLES_TRI.DAT">ANGLES_TRI.DAT</a> - sample trimer angle data

----HOW TO GET STARTED----------------------------------------------

Grab the five files listed above.
Be sure to give them the same file names -- I have hard-wired them
into the code!
Compile the source code -- either BEND.FOR or BEND_TRI.FOR. 
Run the program - it will read angles from ANGLES.DAT (or ANGLES_TRI.DAT)
and the base sequence from SEQUENCE.DAT. It will write diagnostic
material to the terminal, and produce two files: PROFILE.PS is a
postscript profile and PROFILE.LOG contains values for angles.

Our best results were obtained using the program BEND_TRI and the
angles included in ANGLES_TRI.DAT. 

----REFERENCE-----------------------------------------------------------------

Please reference the report listed above in any publication resulting
from use of these programs.
</pre>