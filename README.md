# Guidance

## Mars-Earth Δv minimization
It has been studied through a patched conic method how a local optimizer works in the minimization of the Δv in a 2-impulse manoeuvre Earth to Mars.

## Best bielliptic trajectory Earth-L2(EL2EM)
It has been minimized with MATLAB fmincon the Δv of the Earth-L2(L2EM) in a 4-body propagation approach (spacecraft, Earth, Moon, Sun). As clearly visible in results, the manoeuvering point is close to L2(L2ES).

## Low thrust manoeuvre, minimum time problem
It has been designed the minimum time low thrust manoeuvre through a shooting method where the lagrangian multipliers have been guessed at the initial point and then iterated to find the minimizing vectorial variable.

# Navigation

## Uncertainties propagation
Three models of uncertainties propagation have been applied and compared each others: LinCov, Unscented Transform and Montecarlo simulation (100) points.

## Batch filters 
From data directly coming from a TLE measurements, the navigation problem has been solved with and without a priori information.

## UKF 
Through a one week monitoring campaign from 2 stations of Exomars mission, all the data acquired have been processed to be used for an Unscented Kalman Filter whose estimation provided a very accurate measurement about the state of the satellite.

### Installation
The program uses the basic MATLAB Toolbox and CSPICE (NASA) Toolkit for MATLAB.
