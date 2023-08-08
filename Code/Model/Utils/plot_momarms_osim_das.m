function plot_momarms_osim_das(imus)

% Plots the original Opensim moment arms and the real-time polynomial
% approximations

modelparams = load('simplified_model_struct.mat'); 
model = modelparams.model;

% Get moment arms and lengths out of the .mat files
% These are created using code in the DAS3 repository
% (https://github.com/dasproject/DAS3)
musfilename = ['..\..\DAS3\model\build\simplified_unscaled_model\path_',model.muscles{imus}.name,'.mat'];
ma = load(musfilename);
mus = model.muscles{1,imus};
ndofs = length(mus.dof_indeces); % number of dofs spanned by this muscle

musdof_indeces = zeros(ndofs,1);
for idof = 1:ndofs
    imusdof = mus.dof_indeces(idof);
    musdof_indeces(idof) = imusdof;
end
ang = (ma.alljnts(:,musdof_indeces) + 1e-6);	% protect against angle = 0.0

examine_momarms(mus, mus.dof_names, ma.allmomarms, ang);

end

function examine_momarms(musmodel, dof_names, moment_arms, angles)
% plot momentarm-angle data

% choose a subset of "angles" that contains only 100 points
a = 1;
b = size(angles,1);
r = a + (b-a).*rand(100,1);
indeces = ceil(sort(r));
sangles = angles(indeces,:);

% calculate moment arms from polynomial
pmoment_arms = zeros(length(indeces),length(dof_names));
for iframe = 1:length(indeces)
    for i=1:musmodel.lparam_count

        % add this term's contribution to the muscle length 
        term = musmodel.lcoefs(i);

        for j=1:length(dof_names)
            for k=1:musmodel.lparams(i,j)
                term = term * sangles(iframe,j); % this creates lcoeff(i) * product of all angles to the power lparams(i,j) 
            end
        end

        % first derivatives of length with respect to all q's
        for  k=1:length(dof_names)
            % derivative with respect to q_k is zero unless exponent is 1 or higher and q is not zero
            if ((musmodel.lparams(i,k) > 0) && (sangles(iframe,k)))	
                dterm = musmodel.lparams(i,k)*term/sangles(iframe,k);
                pmoment_arms(iframe,k) = pmoment_arms(iframe,k) + dterm;
            end
        end
    end
end

figure;
for idof=1:length(dof_names)
    subplot(length(dof_names),1,idof);
    plot(moment_arms(indeces,idof),'bx-'); hold on; plot(-pmoment_arms(:,idof),'ro-'); 	
    title([dof_names{idof}, ' momentarms for ',musmodel.name],'Interpreter', 'none'); 
end
legend('osim','poly');
end