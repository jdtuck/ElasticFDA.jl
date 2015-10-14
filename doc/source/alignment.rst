Functional Alignment
====================

The main functions deal with the alignment of functional data using the
*square-root slope (srsf)* framework. Where an input into a function is
expecting an array of functions the shape is assumed to be ``(M,N)`` with ``M``
being the number of sample points and ``N`` being the number of functions.

SRSF Functions
-------------
.. function:: params(d)

    Return a tuple of parameters.

    **Note:** Let ``d`` be a distribution of type ``D``, then
    ``D(params(d)...)`` will construct exactly the same distribution as
     ``d``.

.. function:: succprob(d)

    Get the probability of success.

Alignment
--------

