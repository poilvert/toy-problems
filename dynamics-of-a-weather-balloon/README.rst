Motivation for the model
========================

A group of bright and DIY-type people from the Rowland Institute
at Harvard (actually friends of mine) decided to participate in
the `Hackerspaces in space <http://hackerspacesinspace.com/>`_
competition which consisted in the launch of a *high altitude
weather balloon* equipped with a camera from a smartphone and
some other GPS recording devices. The whole experiment (balloon +
Helium + equipment) was less than $200 !

Given the fact that the team recorded a full battery of physical
variables like time and altitude, I was keen on knowing whether
a simple dynamical model of the balloon flight could describe
with some accuracy the altitude versus time function.

The model
=========

The model that I used was mostly based on my intuition of what
could be the main forces at stake. Also, since the balloon was
to flight up to very high altitudes, I could not disregard the
importance of pressure and temperature with altitude. I describe
below the two important pieces that constitute my model, that is
to say :
    - a model of the atmosphere that the balloon is flying
      through
    - the forces at stake and the equation of motion

Atmospheric model
-----------------

Searching on the web for a **model** of the earth's atmosphere, I
quite quickly found my way to the Wikipedia page of the US standard
atmospheric model_. This model basically assumes a given empirical
temperature profile of the atmosphere up until something like
70,000 meters (that's almost 240,000 feet!). Given that temperature
profile, the model assumes a simple hydrostatic equilibrium and
the perfect gas equation of state to find out what the pressure and
air density profiles are.

I then integrated the model and implemented it inside a module called
``atmosphere.py``. This module contains all the information about
the US standard atmospheric model that I chose and I encourage anyone
to change it / tweak it / play with it. As a sanity check, I display
the temperature and pressure profiles of the model using matplotlib_
if one calls the program without any argument::

    >>> python atmosphere.py

The result is shown below:

.. image:: https://github.com/poilvert/toy-problems/blob/master/dynamics-of-a-weather-balloon/atmosphere.png

Dynamical model
---------------

In a simple model, one should stick with the essentials. In my
dynamical model, I only used what I considered to be the major
forces having an impact on the balloon. Namely I chose to describe
the dynamics of the balloon through:

    1. The gravitational force of the earth (so total mass times
       the acceleration of gravity). Given the highest altitude
       reached which peaked at around 23-24 km, I neglected the
       decrease of the earth gravitational field with altitude.
    2. The buoyancy, which I considered to be exerted only on the
       balloon itself (spherical with a large radius between 7 and
       12 feet) and neglected the payload. This force is pushing
       upward if enough Helium has been injected in the balloon.
       The intensity of the force is the acceleration of gravity
       times the mass of air displaced by the balloon.
    3. At last, given the relatively high Reynolds number for
       the air flow around the balloon at typical velocities
       reached by the balloon when ascending, I had to take into
       account the drag force. I took the standard formula used
       in fluid dynamics for objects moving into a fluid. You can
       find more details here_. The force is basically proportional
       to the area of a section of the balloon (pi times radius
       squared), the density of the air, the square of the balloon
       velocity and a geometric coefficiant called the
       *drag coefficiant*.

The assumptions made by the dynamical model are as follows:

    1. motion only along z (lateral forces completely
       decoupled from the vertical motion)
    2. No mass loss (i.e. no loss of Helium through
       the balloon membrane)
    3. Spherical approximation for the shape of the
       balloon at every altitude.
    4. Two **thermal** models are used. One in which the temperature
       of Helium inside the balloon does not have time to adjust with
       the outside temperature of the air. We could call this the
       **adiabatic** model. A second model in which the ascent is
       sufficently slow so as to make the temperature of Helium to
       be equal to the temperature of the outside air. We could call
       that latter model a **quasi-static** model.

Running the program
===================

In order to run the simulation, just type::

    >>> python dynamics.py

What the code does is this:

    1. It parses the experimental data and extracts the time and
       altitude values for which we have some information.
       This part is taken care of by ``parse_raw_data.py``.
    2. It loads the atmospheric model from ``atmosphere.py``.
    3. It optimizes some free parameters in the dynamical model
       to find a best fit to the experimental data. This last
       part is done by ``dynamics.py``.

The program then prints the values for the *optimal* free parameters
and plots the theoretical solution and the experimental one on the
same graph.

An example of the kind of plot that one obtains is given below:

.. image:: https://github.com/poilvert/toy-problems/blob/master/dynamics-of-a-weather-balloon/balloon_flight_dynamics.png

Discussion
==========

As one can readily see on the previous plot, this small model is pretty
decent in predicting the dynamics of the balloon. Nevertheless, it seems
that something strange happens at about 11,000 meters, where the theory
starts to deviate from experiment.
Two possibilities came to my mind. First of all, the atmospheric model
has an abrupt change in temperature profile at exactly 11,000 meters. This
could potentially explain the discrepancy above that altitude.
A second explanation, that is actually related to the first one, is the
onset of the **jet streams** layer_, which could have rendered the
"upward motion only" hypothesis completely wrong (because of the strong
turbulence brought about by the jet streams impacting the motion of the
balloon).

.. external links

.. _model: http://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
.. _matplotlib: http://matplotlib.sourceforge.net
.. _here: http://en.wikipedia.org/wiki/Drag_coefficient
.. _layer: http://en.wikipedia.org/wiki/Jet_stream
