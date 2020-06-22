function [Pdf] = DistributionMRImages(Sequences,Parameters,Ydens)


if nargin < 3, error('Not enought input arguments'); end

switch length(size(Sequences))
    case 4
        [s1,s2,t,slices] = size(Sequences);
    case 3
        [s1,s2,t]   = size(Sequences);
        slices      = 1;
    case 2
        [s1,t]      = size(Sequences);
        s2          = 1;
        Sequences   = reshape(Sequences, s1,s2,t);
        slices      = 1;
        if ~isempty(References) && length(size(References))==2
            References = reshape(References, s1,s2,size(References,2));
        end
    otherwise
        error('Invalid Sequences argument size')
end


for s = 1:slices
    
    Obs = reshape(Sequences(:,:,:,s),s1*s2,t);
    [m,~,~,alpha] = gllim_inverse_map_and_kurt(Obs', Parameters.theta, 0);
    
    for i = 1:size(Obs,1)
        [~,psi]         = gllim_inverse_dens(Obs(i,:)', Parameters.theta, ones(size(Obs(i,:)')), [], 0);
        
        for c = 1:size(psi.mu,1)
            gm          = gmdistribution(psi.mu(c,alpha(i,:)>0)', psi.S(c,c,alpha(i,:)>0), alpha(i,alpha(i,:)>0));
            p(i,c,:)	= pdf(gm, Ydens{c}');
        end    
    end
    
    Pdf(:,:,:,:,s) = reshape(p, s1,s2,size(p,2),size(p,3));
end




