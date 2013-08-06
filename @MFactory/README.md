MFactory
========

Represents a Semiempirical Model, via a set of mixers that it contains and manages. These mixers establish the "policies" (how are matrix elements modified) and the parameters that do the adjusting.

Given an MSet, is able to produce a fitme corresponding to that data set.

The setPolicies method of MFactory is where the decisions on the nature of the models are made, and we give these names. The goal is to add to this, but to never change the previous models, since we want to be able to say that "hybridslater1" worked well on ethane. If you want to change the nature of the decisions, call it hybridslater2.

Policies represent rules for attaching mixers to a model. Decisions are of the form:

* What matrix elements does a given mixer modify:
    * Operator (KE, EN(Z), E2).
    * Diag or off-diagonal matrix elements.
    * Element types: 'i' for diagonal and 'j' for off-diagonal.
* What functional form is used for the modifications:
    * Multiply by constant, add a constant, or interpolate between narrow and diffuse.
