# SolarCalculator
Solar Calculator that predicts time-varying values of incident radiation and temperature worldwide. It  can be coupled to photoactivated kinetic models (SODIS in the example) to estimate its potential

The solar positioning uses the NREL SPA algorithm (spa.c and spa.h), which must be compiled and linked alltogether with the main code file, SolarCalculator.cpp

The program includes a kinetic model of SODIS (Solar Disinfection), based on Hom's law. It can be modified to study the potential impact of any process requiring instant values of temperature and/or incident radiation. To ease this modification, a #CHANGEKINETICS flag has been added in several locations of the code.
