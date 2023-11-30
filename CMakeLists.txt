cmake_minimum_required(VERSION 3.26)
project(ProjectFiniteElements C)

set(CMAKE_C_STANDARD 11)

add_executable(ProjectFiniteElements
        errorFunction.m
        FEMDiNe.m
        FEMDiNeQuadratura.m
        FEMDiNeQuadratura2.m
        FEMDiNoHom.m
        FEMParabolicEquations.m
        integralQuadratura.m
        main.m
        nodes_weights.m
        P2.m
        Triangolator.m)
