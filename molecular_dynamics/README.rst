What is this code?
==================

This code illustrates the *Verlet Algorithm* used in **Molecular Dynamics** simulations.

What is the model?
==================

For the purpose of this demo code, I used a very simple Lennard-Jones pair potential.

How to use the code?
====================

Just launch the script::

    python molecular_dynamics.py

And a window will pop up, from which you can use your mouse to select the positions of as many atom as you'd like for the MD run. Then close the window and the MD run will start. At the very end, a plot of the different terms in the total energy of the system will be displayed. The key feature to look for here is the constancy of the total energy.
