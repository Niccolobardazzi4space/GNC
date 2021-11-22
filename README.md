# GNC

## Mars-Earth Δv minimization
It has been studied through a patched conic method how a local optimizer works in the minimization of the Δv in a 2-impulse manoeuvre Earth to Mars.

## Best bielliptic trajectory Earth-L2(EL2EM)
It has been minimized with MATLAB fmincon the Δv of the Earth-L2(L2EM) in a 4-body propagation approach (spacecraft, Earth, Moon, Sun). As clear visible in results, the manoeuvering point is close to L2(L2ES).

## Low thrust manoeuvre, minimum time problem
It has been designed the minimum time low thrust manoeuvre through a shooting method where the lagrangian multipliers have been guessed at the initial point and then iterated to find the minimizing vectorial variable.

### Installation
The program uses the basic MATLAB Toolbox and CSPICE (NASA) Toolkit for MATLAB.
