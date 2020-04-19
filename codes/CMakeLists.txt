cmake_minimum_required(VERSION 3.6)
project(Samaritan)

set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES
        header/analyse.h
        header/dominance.h
        header/global.h
        header/indicators.h
        header/initialization.h
        header/memory.h
        header/metaheuristics.h
        header/population.h
        header/print.h
        header/problems.h
        header/rand.h
        header/rank_sort.h
        header/reproduction.h
        header/selection.h
        header/utility.h
        header/vector.h

        core/dominance.c
        core/initialization.c
        core/list.c
        core/memory.c
        core/population.c
        core/print.c
        core/rand.c
        core/rank_sort.c
        core/utility.c
        core/vector.c
        core/g_evaluate.c

        metaheuristics/gnsga2.c
        metaheuristics/r2nsga2.c
        metaheuristics/pbea.c
        metaheuristics/pbeaGeneral.c
        metaheuristics/R-MEAD2.c
        metaheuristics/nsga2.c
        metaheuristics/moead.c
        metaheuristics/ibea.c
        metaheuristics/nsga3.c
        metaheuristics/rnsga2.c
        metaheuristics/rnsga2General.c

        selection/ibea/fitness.c
        selection/ibea/ibea_selection.c
        selection/nsga2/crowddist.c
        selection/nsga2/fillnds.c
        selection/nsga2/qsort.c
        selection/moead/update.c
        selection/moead/utility.c
        selection/moead/fitness.c

        reproduction/crossover/crossover_moead.c
        reproduction/crossover/crossover_nsga2.c
        reproduction/crossover/crossover_spea2.c
        reproduction/crossover/crossover_steadystate.c
        reproduction/crossover/differential_evolution.c
        reproduction/crossover/sbx.c
        reproduction/mutation/mutation_real.c
        reproduction/mutation/polymut.c
        reproduction/mating/moead_mating.c
        reproduction/mating/tournament_selection.c

        analyse/analyse.c
        analyse/plot.c

        indicators/gd.c
        indicators/hv.c
        indicators/igd.c

        externals/WFG_1.15/read.c
        externals/WFG_1.15/wfg.c
        externals/WFG_1.15/wfg.h
        externals/IWFG/iwfg.c
        externals/IWFG/iwfg.h
        externals/IWFG/read.c

        testproblem/evaluation.c
        testproblem/stock/stock.c
        testproblem/ZDT/ZDT1.c
        testproblem/ZDT/ZDT2.c
        testproblem/ZDT/ZDT3.c
        testproblem/ZDT/ZDT4.c
        testproblem/ZDT/ZDT6.c
        testproblem/DTLZ/DTLZ1.c
        testproblem/DTLZ/DTLZ2.c
        testproblem/DTLZ/DTLZ3.c
        testproblem/DTLZ/DTLZ4.c
        testproblem/DTLZ/DTLZ5.c
        testproblem/DTLZ/DTLZ6.c
        testproblem/DTLZ/DTLZ7.c
        testproblem/mDTLZ/mDTLZ1.c
        testproblem/mDTLZ/mDTLZ2.c
        testproblem/mDTLZ/mDTLZ3.c
        testproblem/mDTLZ/mDTLZ4.c
        testproblem/minusDTLZ/minusDTLZ1.c
        testproblem/minusDTLZ/minusDTLZ2.c
        testproblem/minusDTLZ/minusDTLZ3.c
        testproblem/minusDTLZ/minusDTLZ4.c
        testproblem/UF/UF1.c
        testproblem/UF/UF2.c
        testproblem/UF/UF3.c
        testproblem/UF/UF4.c
        testproblem/UF/UF5.c
        testproblem/UF/UF6.c
        testproblem/UF/UF7.c
        testproblem/UF/UF8.c
        testproblem/UF/UF9.c
        testproblem/UF/UF10.c
        testproblem/CDTLZ/C1DTLZ1.c
        testproblem/CDTLZ/C2DTLZ2.c
        testproblem/CDTLZ/C1DTLZ3.c
        testproblem/CDTLZ/C3DTLZ1.c
        testproblem/CDTLZ/C3DTLZ4.c

        testproblem/WFG/toolkit.h
        testproblem/WFG/toolkit.c
        testproblem/WFG/WFG1.c
        testproblem/WFG/WFG2.c
        testproblem/WFG/WFG3.c
        testproblem/WFG/WFG4.c
        testproblem/WFG/WFG5.c
        testproblem/WFG/WFG6.c
        testproblem/WFG/WFG7.c
        testproblem/WFG/WFG8.c
        testproblem/WFG/WFG9.c
        testproblem/WFG/WFG4x.c
        config.txt

        main.c
        )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(Samaritan ${SOURCE_FILES})
set(CMAKE_C_FLAGS  ${CMAKE_C_FLAGS} -pthread)
target_link_libraries(Samaritan m)
