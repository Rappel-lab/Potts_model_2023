# Potts_model_2023
This is a Cellular Potts Model (CPM) that simulates the development and migration pattern of cancer cells in dense extra-cellular matrix (ECM). In addition to cell-cell adhesion/cell-medium interfacial tension/cell size constraint in conventional CPM, we've built additional modules to model cell growth and division, cell motility and cell-ECM interaction which incluces local and globlal remodeling.

`CellularPottsModel_CancerDevMigration.m` is the main code which simulates the development of cancer cell from a single cell to multi-cellular structure. By adjusting the modes of cell-ECM remodeling (`q.rho`, `q.KD1`) and cell polarization (`q.Kv`, `q.Kp`), we observe the two distinct modes of cell morphologies (spheroids, network) and migration patterns (rotaional,s invasive).

`CellularPottsModel_CancerDevMigration_Parameter2DGridSearch.m` is the code to conduct 2D parameter grid search over a given pair of parameters. `CPM_Simulations_Analysis.m` takes the outputs from the simulation to compute statistics of cell speed, circularity, protrusion dynamics etc., of the simulated cells.
