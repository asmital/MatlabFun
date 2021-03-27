
expectation_r=[];
expectation_rsq=[];
expectation_1_by_r=[];
expectation_1_by_rsq=[];
reduced_mass=1836/1837;
for n=1:3
    for l=0:n-1
        k=reduced_mass/n;
        expectation_r=[expectation_r;integral(@(r) r.^3.*wavefunction(r,k,l,n).*wavefunction(r,k,l,n), 0, inf)];
    end
end
for n=1:3
    for l=0:n-1
        k=reduced_mass/n;
        expectation_rsq=[expectation_rsq;integral(@(r) r.^4.*wavefunction(r,k,l,n).*wavefunction(r,k,l,n), 0, inf)];
    end
end

for n=1:3
    for l=0:n-1
        k=reduced_mass/n;
        expectation_1_by_r=[expectation_1_by_r;integral(@(r) r.*wavefunction(r,k,l,n).*wavefunction(r,k,l,n), 0, inf)];
    end
end

for n=1:3
    for l=0:n-1
        k=reduced_mass/n;
        expectation_1_by_rsq=[expectation_1_by_rsq;integral(@(r) wavefunction(r,k,l,n).*wavefunction(r,k,l,n), 0, inf)];
    end
end
Row_headings=[1; 2; 3;4;5;6];
n_values=[1;2;2;3;3;3];
l_values=[0;0;1;0;1;2];
T=table(Row_headings,n_values, l_values,expectation_r, expectation_rsq,expectation_1_by_r,expectation_1_by_rsq);
T.Properties.VariableNames ={'Case no','n values','l values','E[r]','E[r^2]', 'E[1/r]', 'E[1/(r^2)]'};
disp(T)


function psi=wavefunction(r,k,l,n)
c_R=normalisation(r,k,l,n);
R=radialwavefunction(r,k,l,n);
psi=c_R.*R;
end
function R=radialwavefunction(r,k,l,n)
R=r.^l.*exp(-k.*r);
series=0;
j=n-l-1;
for i=0:j
    series=series+(recurrence(l,i,n)).*(r.^i).*(k.^i);
end
R=R.*series;
end
function c_R=normalisation(r,k,l,n)
temp=integral(@(r)(radialwavefunction(r,k,l,n).*r).^2, 0, inf);
c_R=1/(sqrt(temp));
end

function cj= recurrence(l,j,n)
    if(j==0)
        cj=1;
        return
    else
        j=j-1;
        temp=(2*(j+l+1)-2*n)/((j+2*l+2)*(j+1));
        cj=temp*recurrence(l,j,n);
    end
    return
end







       