function add_parameters_dofs(model_struct)

prev_model = load(model_struct);
new_model = prev_model;

% Add stiffness K1, K2 and damping to each DOF
for idof = 1:new_model.model.nDofs
    new_model.model.dofs{idof}.stiffness_K1 = 1;
    new_model.model.dofs{idof}.stiffness_K2 = 5000;
    new_model.model.dofs{idof}.B = 1;
end

% Set the linear elastic moment to zero for elbow flexion/extension
new_model.model.dofs{13}.stiffness_K1 = 0;

model = new_model.model;
save(model_struct,'model');