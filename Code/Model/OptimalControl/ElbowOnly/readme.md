## Optimisation of elbow flexion-extension

First, make sure you add the folders **TechForParalysis\Code\Model** and **TechForParalysis\Code\Model\OptimalControl** to the Matlab path.

To choose which muscles to include in the optimisation, open **include_subset_muscles.m** and edit line 11 to select the elements to include. You must then run this file to update the structure before running the model. You might want to run **get_muscle_names.m** first, so you can see which indices correspond to which muscles.

You run the optimisation in **runopt_elbow.m**. In this file you can set:

- the starting number of nodes
- maximum number of nodes
- initial guess
- weights for the kinematic and energy consumption terms
- target movement. You can choose a movement or a static posture (for which # nodes = 1).

At the end of each IPOPT run, it prints out something like this:

IPOPT info status: 1

Objfun (unweighted):   0.23037 (fit) +   0.14236 (effort)

Objfun (weighted):  23.17969 =  23.03733 (fit) +   0.14236 (effort)

Norm(ceq):   0.00056

The first line shows the reason for termination: 0 means "solved", 1 means "solved to acceptable level", -1 means "maximum number of iterations exceeded". For the rest, see [the Matlab Interface for IPOPT](https://ethz.ch/content/dam/ethz/special-interest/mavt/dynamic-systems-n-control/idsc-dam/Research_Onder/Downloads/IPOPT/IPOPT_MatlabInterface_V0p1.pdf), section 4.

**Objfun** is the value of the objective function, as a sum of the two terms, first unweighted, and then weighted. In the example above, the weights are 100 for the kinematic fit, and 1 for the effort.

**Norm(ceq)** is the magnitude of the linear equality constraints.

You can run **print_results_obj.m** to make graphs of the output elbow angle and muscle activations.

If you want to call the model (das3.mexw64) to get, for example, muscle moment arms or lengths for additional plots, see [das3.m in the DAS3 repository](https://github.com/dasproject/DAS3/blob/master/model/das3.m).
