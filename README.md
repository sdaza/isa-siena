## Agent-Based Models for Assessing Complex Statistical Models: An Example Evaluating Selection and Social Influence Estimates from SIENA

Code and data to reproduce results.

- Authors: Sebastian Daza and L. Kurt Kreuger
- ISA: Influence Selection Agent-based Model
- SIENA: [Simulation Investigation for Empirical Network Analysis](https://www.stats.ox.ac.uk/~snijders/siena/)

### Anylogic

 - [Anylogic simulation file](https://github.com/sdaza/isa-siena/tree/master/Anylogic) (.alp)
 - All models can be run using the [**Personal Learning Edition (PLE) Anylogic Version**](https://www.anylogic.com/downloads/).
 - Data can be generated running the experiments defined in the Anylogic file. They used the same random seed to create replicates.

### R

 - [R Code](https://github.com/sdaza/isa-siena/tree/master/R) used to process the data generated in [**Anylogic**](https://www.anylogic.com/downloads/).
 - [RData files](https://github.com/sdaza/isa-siena/tree/master/data) are [data.table](https://github.com/Rdatatable/data.table/wiki) objects containing the data generated using Anylogic. The data.table are:
    - `b`: behavior data
    - `d`: distance between agents
    - `n`: network data
    - `p`: experimental data
